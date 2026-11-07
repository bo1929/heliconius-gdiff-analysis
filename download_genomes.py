#!/usr/bin/env python3
"""
Download every Heliconius assembly listed in genomes.tsv, with annotations,
metadata, assembly stats and assembly reports.

Per genome this fetches, from the NCBI FTP mirror:

  <CID>.<acc>.fna.gz                  genomic FASTA
  <CID>.<acc>.gff.gz                  gene annotation  (preferred)
  <CID>.<acc>.gbff.gz                 ...or GenBank flat file when no GFF exists
  <CID>.<acc>.assembly_stats.txt      assembly statistics (N50s, counts, GC)
  <CID>.<acc>.assembly_report.txt     full assembly report (sequence-level table)

Note on annotations: most Heliconius assemblies here are GenBank scaffold
submissions that carry **no GFF**, only a GBFF. Asking for "annotations if
exists" while only requesting GFF would therefore leave 14 of 16 genomes with
nothing, so the script prefers GFF and falls back to GBFF. The choice made for
each genome is recorded in ``genomes/metadata.tsv`` (``annotation`` column).

It also writes:

  genomes/manifest.tsv    every file downloaded, with sizes and source URLs
  genomes/metadata.tsv    one row per genome: organism, accession, level,
                          submitter, bioproject, WGS, assembly method,
                          total length, GC%, scaffold/contig N50 and counts

Every file is verified against the assembly's published ``md5checksums.txt``
before being accepted, and re-running the script skips files already present and
correct.

Requirements: python3 (stdlib only) and curl.

Examples
--------
    ./download_genomes.py --dry-run            # show the plan, download nothing
    ./download_genomes.py                      # everything
    ./download_genomes.py --only Heliconius_sara,Heliconius_erato
    ./download_genomes.py --no-annotation      # genomes, stats and reports only
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
USER_AGENT = "gdiff-heliconius-analysis/1.0 (genome download; contact: local user)"
ACCESSION_RE = re.compile(r"^(GC[AF])_(\d{9,})\.(\d+)$")

_print_lock = threading.Lock()


def log(msg: str) -> None:
    with _print_lock:
        print(msg, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------- #
# input
# --------------------------------------------------------------------------- #


@dataclass
class Genome:
    cid: str
    taxonomy: str
    accession: str
    level: str = ""
    link: str = ""
    asm_name: str = ""
    asm_dir: str = ""


def read_genomes(tsv_path: str) -> list[Genome]:
    """Parse genomes.tsv (CID / Taxonomy / ASSEMBLY_ACCESSION / ...)."""
    genomes: list[Genome] = []
    malformed: list[str] = []
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"CID", "ASSEMBLY_ACCESSION"}
        if not required.issubset(reader.fieldnames or []):
            sys.exit(f"error: {tsv_path} must have columns: {sorted(required)}")
        for lineno, row in enumerate(reader, start=2):
            cid = (row.get("CID") or "").strip()
            acc = (row.get("ASSEMBLY_ACCESSION") or "").strip()
            if not cid or not acc:
                continue
            if not ACCESSION_RE.match(acc):
                malformed.append(f"  line {lineno}: {cid}: bad accession {acc!r}")
                continue
            genomes.append(
                Genome(
                    cid=cid,
                    taxonomy=(row.get("Taxonomy") or "").strip(),
                    accession=acc,
                    level=(row.get("LEVEL") or "").strip(),
                    link=(row.get("LINK") or "").strip(),
                    asm_name=(row.get("ASSEMBLY_NAME") or "").strip(),
                )
            )
    if malformed:
        sys.exit("error: unparseable accession(s):\n" + "\n".join(malformed))
    if not genomes:
        sys.exit(f"error: no genomes found in {tsv_path}")
    return genomes


# --------------------------------------------------------------------------- #
# NCBI helpers
# --------------------------------------------------------------------------- #


def _sleep_for(attempt: int, retry_after: str | None) -> None:
    if retry_after:
        try:
            time.sleep(min(float(retry_after), 60.0))
            return
        except ValueError:
            pass
    time.sleep(min(2.0 * (2**attempt), 30.0))


def http_get(url: str, tries: int = 6) -> bytes:
    """GET a URL, retrying the transient errors NCBI throws under load."""
    last: Exception | None = None
    for attempt in range(tries):
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code == 404:      # a real answer, not a transient failure
                raise
            _sleep_for(attempt, exc.headers.get("Retry-After"))
        except Exception as exc:
            last = exc
            _sleep_for(attempt, None)
    raise RuntimeError(f"GET failed after {tries} tries: {url} ({last})")


def accession_path(accession: str) -> str:
    """GCA_056956725.1 -> GCA/056/956/725 (NCBI's 3-digit fan-out)."""
    prefix, rest = accession.split("_", 1)
    digits = rest.split(".", 1)[0]
    groups = [digits[i : i + 3] for i in range(0, len(digits), 3)]
    return "/".join([prefix] + groups)


def resolve_assembly_dir(genome: Genome) -> tuple[str, str]:
    """Find the assembly directory name, e.g. ('GCA_..._ASM5695672v1', url)."""
    parent = f"{FTP_ROOT}/{accession_path(genome.accession)}/"
    html = http_get(parent).decode("utf-8", "replace")
    match = re.search(r'href="(' + re.escape(genome.accession) + r'_[^"/]+)/"', html)
    if not match:
        raise RuntimeError(f"no assembly directory for {genome.accession} under {parent}")
    asm_dir = match.group(1)
    return asm_dir, f"{parent}{asm_dir}/"


def list_files(dir_url: str) -> set[str]:
    """File names offered in an assembly directory (fallback inventory)."""
    html = http_get(dir_url).decode("utf-8", "replace")
    return set(re.findall(r'href="([^"/]+\.(?:gz|txt))"', html))


def fetch_md5s(dir_url: str) -> dict[str, str]:
    """Return {filename: md5} from the assembly's md5checksums.txt."""
    text = http_get(f"{dir_url}md5checksums.txt").decode("utf-8", "replace")
    checksums: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2 and re.fullmatch(r"[0-9a-f]{32}", parts[0]):
            checksums[os.path.basename(parts[1])] = parts[0]
    return checksums


# --------------------------------------------------------------------------- #
# downloading
# --------------------------------------------------------------------------- #


def md5_of(path: str, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def curl_download(url: str, dest: str) -> None:
    """Download url -> dest (atomically), resuming a partial .part file."""
    part = dest + ".part"
    cmd = [
        "curl", "--fail", "--location", "--silent", "--show-error",
        "--retry", "5", "--retry-delay", "5", "--retry-all-errors",
        "--connect-timeout", "30",
        "--speed-limit", "1024", "--speed-time", "120",
        "--user-agent", USER_AGENT,
    ]
    if os.path.exists(part) and os.path.getsize(part) > 0:
        cmd += ["--continue-at", "-"]
    cmd += ["--output", part, url]
    if subprocess.run(cmd).returncode != 0:
        raise RuntimeError(f"curl failed for {url}")
    os.replace(part, dest)


@dataclass
class Result:
    genome: Genome
    status: str = "skipped"
    files: list[tuple[str, int]] = field(default_factory=list)
    error: str = ""
    url: str = ""
    fetched: bool = False
    annotation: str = "none"


def wanted_files(genome: Genome, inventory: set[str], asm_dir: str,
                 args: argparse.Namespace) -> list[tuple[str, str]]:
    """[(remote file name, local file name)] for one genome."""
    stem = f"{genome.cid}.{genome.accession}"
    items = [(f"{asm_dir}_genomic.fna.gz", f"{stem}.fna.gz"),
             (f"{asm_dir}_assembly_stats.txt", f"{stem}.assembly_stats.txt"),
             (f"{asm_dir}_assembly_report.txt", f"{stem}.assembly_report.txt")]
    if not args.no_annotation:
        # Prefer a real GFF; fall back to the GenBank flat file, which is the
        # only annotation these scaffold assemblies usually ship.
        for suffix, ext in (("_genomic.gff.gz", "gff.gz"), ("_genomic.gbff.gz", "gbff.gz")):
            remote = f"{asm_dir}{suffix}"
            if remote in inventory:
                items.append((remote, f"{stem}.{ext}"))
                break
    return items


def process(genome: Genome, args: argparse.Namespace) -> Result:
    res = Result(genome=genome)
    try:
        asm_dir, dir_url = resolve_assembly_dir(genome)
        genome.asm_dir = asm_dir
        genome.asm_name = genome.asm_name or asm_dir[len(genome.accession) + 1 :]
        res.url = dir_url
        log(f"[{genome.cid}] {genome.accession} -> {asm_dir}")

        try:
            checksums = fetch_md5s(dir_url)
        except Exception as exc:
            checksums = {}
            log(f"[{genome.cid}] warning: no md5checksums.txt ({exc})")
        inventory = set(checksums)
        if not inventory:
            try:
                inventory = list_files(dir_url)
                log(f"[{genome.cid}] note: using directory listing as inventory "
                    "(no checksums available)")
            except Exception as exc:
                raise RuntimeError(f"cannot determine available files: {exc}")

        for remote_name, local_name in wanted_files(genome, inventory, asm_dir, args):
            dest = os.path.join(args.out, local_name)
            expected = checksums.get(remote_name)
            if remote_name not in inventory:
                log(f"[{genome.cid}] note: {remote_name} not offered by NCBI; skipping")
                continue
            if local_name.endswith((".gff.gz", ".gbff.gz")):
                res.annotation = "gff" if local_name.endswith(".gff.gz") else "gbff"

            if args.dry_run:
                log(f"[{genome.cid}] plan: {remote_name} -> {local_name}")
                res.files.append((local_name, 0))
                continue

            if os.path.exists(dest) and not args.force:
                if args.no_verify or expected is None:
                    log(f"[{genome.cid}] have {local_name} (unverified)")
                    res.files.append((local_name, os.path.getsize(dest)))
                    continue
                if md5_of(dest) == expected:
                    log(f"[{genome.cid}] ok   {local_name} (already present, md5 verified)")
                    res.files.append((local_name, os.path.getsize(dest)))
                    continue
                log(f"[{genome.cid}] stale/corrupt {local_name}; re-downloading")

            log(f"[{genome.cid}] get  {remote_name}")
            curl_download(f"{dir_url}{remote_name}", dest)
            res.fetched = True
            if expected is not None:
                actual = md5_of(dest)
                if actual != expected:
                    os.unlink(dest)
                    raise RuntimeError(
                        f"md5 mismatch for {local_name}: expected {expected}, got {actual}"
                    )
            size = os.path.getsize(dest)
            res.files.append((local_name, size))
            log(f"[{genome.cid}] done {local_name} ({size / 1e6:.1f} MB)")

        res.status = "planned" if args.dry_run else (
            "downloaded" if res.fetched else "skipped"
        )
    except Exception as exc:  # one bad genome must not kill the run
        res.status = "failed"
        res.error = str(exc)
        log(f"[{genome.cid}] FAILED: {exc}")
    return res


# --------------------------------------------------------------------------- #
# metadata
# --------------------------------------------------------------------------- #

METADATA_COLUMNS = [
    "CID", "taxonomy", "assembly_accession", "assembly_name", "organism", "isolate",
    "taxid", "biosample", "bioproject", "submitter", "release_date",
    "assembly_type", "assembly_level", "genome_representation", "wgs_project",
    "assembly_method", "genome_coverage", "sequencing_technology",
    "total_length", "gc_percent", "scaffold_count", "scaffold_N50", "scaffold_L50",
    "contig_count", "contig_N50", "contig_L50", "molecule_count", "annotation",
    "source_url",
]

# "# Key:  value" lines in assembly_stats.txt -> metadata column
STATS_HEADER_KEYS = {
    "Assembly name": "assembly_name",
    "Organism name": "organism",
    "Isolate": "isolate",
    "Taxid": "taxid",
    "BioSample": "biosample",
    "BioProject": "bioproject",
    "Submitter": "submitter",
    "Date": "release_date",
    "Assembly type": "assembly_type",
    "Assembly level": "assembly_level",
    "Genome representation": "genome_representation",
    "WGS project": "wgs_project",
    "Assembly method": "assembly_method",
    "Genome coverage": "genome_coverage",
    "Sequencing technology": "sequencing_technology",
}

# assembly-level ("all") statistics -> metadata column
STATS_VALUE_KEYS = {
    "total-length": "total_length",
    "gc-perc": "gc_percent",
    "scaffold-count": "scaffold_count",
    "scaffold-N50": "scaffold_N50",
    "scaffold-L50": "scaffold_L50",
    "contig-count": "contig_count",
    "contig-N50": "contig_N50",
    "contig-L50": "contig_L50",
    "molecule-count": "molecule_count",
}


def parse_assembly_stats(path: str) -> dict[str, str]:
    """Read the descriptive header and the assembly-wide ('all') statistics."""
    out: dict[str, str] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") and ":" in line:
                key, _, val = line.lstrip("# ").partition(":")
                col = STATS_HEADER_KEYS.get(key.strip())
                if col:
                    out[col] = val.strip()
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) != 6:
                continue
            unit, _mol, _type, role, stat, value = fields
            # The "all" rows are the whole-assembly summary.
            if unit == "all" and role == "all" and stat in STATS_VALUE_KEYS:
                out.setdefault(STATS_VALUE_KEYS[stat], value)
    return out


def build_metadata(out_dir: str, results: list[Result]) -> list[dict]:
    rows: list[dict] = []
    for res in results:
        g = res.genome
        stem = f"{g.cid}.{g.accession}"
        row = {
            "CID": g.cid,
            "taxonomy": g.taxonomy,
            "assembly_accession": g.accession,
            "assembly_name": g.asm_name,
            "assembly_level": g.level,
            "annotation": res.annotation,
            "source_url": res.url or g.link,
        }
        stats_path = os.path.join(out_dir, f"{stem}.assembly_stats.txt")
        if os.path.exists(stats_path):
            row.update({k: v for k, v in parse_assembly_stats(stats_path).items() if v})
        rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

MANIFEST_COLUMNS = ["CID", "taxonomy", "assembly_accession", "assembly_name",
                    "assembly_dir", "status", "file", "bytes", "url"]


def merge_manifest(manifest_path: str, out_dir: str, results: list[Result]) -> list[dict]:
    """Merge this run into the existing manifest and drop entries gone from disk."""
    rows: dict[tuple[str, str], dict] = {}
    if os.path.exists(manifest_path):
        with open(manifest_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("CID") and row.get("file"):
                    rows[(row["CID"], row["file"])] = row

    for res in results:
        for local_name, _ in res.files:
            path = os.path.join(out_dir, local_name)
            if not os.path.exists(path):
                continue
            rows[(res.genome.cid, local_name)] = {
                "CID": res.genome.cid,
                "taxonomy": res.genome.taxonomy,
                "assembly_accession": res.genome.accession,
                "assembly_name": res.genome.asm_name,
                "assembly_dir": res.genome.asm_dir,
                "status": res.status,
                "file": local_name,
                "bytes": str(os.path.getsize(path)),
                "url": res.url,
            }

    return [rows[key] for key in sorted(rows)
            if os.path.exists(os.path.join(out_dir, rows[key]["file"]))]


def parse_args(argv: list[str]) -> argparse.Namespace:
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(
        prog="download_genomes.py (Heliconius)",
        description="Download Heliconius genomes with annotations, metadata, stats and reports.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--tsv", default=os.path.join(here, "genomes.tsv"),
                        help="genome table (default: genomes.tsv next to this script)")
    parser.add_argument("--out", default=os.path.join(here, "genomes"),
                        help="output directory (default: genomes/ next to this script)")
    parser.add_argument("--jobs", type=int, default=4, help="parallel downloads (default: 4)")
    parser.add_argument("--only", default="", help="comma-separated CIDs to include")
    parser.add_argument("--exclude", default="", help="comma-separated CIDs to skip")
    parser.add_argument("--no-annotation", action="store_true",
                        help="skip gene annotations (GFF/GBFF)")
    parser.add_argument("--force", action="store_true",
                        help="re-download even if a verified file is present")
    parser.add_argument("--no-verify", action="store_true",
                        help="skip md5 verification of existing files")
    parser.add_argument("--dry-run", action="store_true",
                        help="resolve accessions and print the plan only")
    args = parser.parse_args(argv)

    if args.jobs < 1:
        parser.error("--jobs must be >= 1")
    for name in ("only", "exclude"):
        setattr(args, name, {t.strip() for t in getattr(args, name).split(",") if t.strip()})
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])

    if shutil.which("curl") is None:
        sys.exit("error: curl is required but was not found on PATH")

    genomes = read_genomes(args.tsv)
    known = {g.cid for g in genomes}
    for flag in ("only", "exclude"):
        unknown = sorted(getattr(args, flag) - known)
        if unknown:
            sys.exit(f"error: --{flag} has unknown CID(s): {', '.join(unknown)}")
    if args.only:
        genomes = [g for g in genomes if g.cid in args.only]
    if args.exclude:
        genomes = [g for g in genomes if g.cid not in args.exclude]
    if not genomes:
        sys.exit("error: no genomes selected")

    if not args.dry_run:
        os.makedirs(args.out, exist_ok=True)
    log(f"downloading {len(genomes)} genome(s) into {args.out} "
        f"with {args.jobs} parallel job(s)" + (" [dry run]" if args.dry_run else ""))

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda g: process(g, args), genomes))

    manifest_rows: list[dict] = []
    if not args.dry_run:
        manifest = os.path.join(args.out, "manifest.tsv")
        manifest_rows = merge_manifest(manifest, args.out, results)
        with open(manifest, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(manifest_rows)
        log(f"wrote {manifest} ({len(manifest_rows)} file(s) tracked)")

        meta_rows = build_metadata(args.out, results)
        if meta_rows:
            meta_path = os.path.join(args.out, "metadata.tsv")
            with open(meta_path, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=METADATA_COLUMNS, delimiter="\t",
                                   extrasaction="ignore")
                w.writeheader()
                w.writerows(meta_rows)
            log(f"wrote {meta_path} ({len(meta_rows)} genome(s))")

    failed = [r for r in results if r.status == "failed"]

    log("")
    log(f"{'CID':<30}{'accession':<18}{'status':<12}{'ann':<6}files")
    for res in results:
        log(f"{res.genome.cid:<30}{res.genome.accession:<18}{res.status:<12}"
            f"{res.annotation:<6}{len(res.files)}")
    if not args.dry_run:
        total = sum(int(r["bytes"]) for r in manifest_rows)
        log(f"\ntotal data on disk: {total / 1e9:.2f} GB across {len(manifest_rows)} file(s)")

    if failed:
        log(f"\n{len(failed)} genome(s) failed:")
        for res in failed:
            log(f"  {res.genome.cid} ({res.genome.accession}): {res.error}")
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit("\ninterrupted; re-run to resume partial downloads")
