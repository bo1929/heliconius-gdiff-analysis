#!/usr/bin/env python3
"""Build the uniform ANI-comparison tables and merge them into one long table.

    outputs/ + logs/  ->  tables/pwdist-<method>-<config>.tsv
                                              samples/samples-gdiff-<config>.tsv
                                              tables/resources.tsv
                                              merged-distances.tsv

``distance`` is 0-1 (0 = identical) and genome IDs are ``<CID>.<assembly>``.  dashing2
reports a similarity matrix and is converted here; the merge re-checks every other table
really is a distance by correlating it against ANIm.

    python3 ani-tables.py
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from itertools import groupby

SRC_OUT, SRC_FASTANI = "outputs", "outputs/fastani-results"
OUT_DIST = "dist"
OUT_TABLES, OUT_SAMPLES, MERGED = "tables", "samples", "merged-distances.tsv"
TABLE_HEADER = ["genome_a", "genome_b", "distance"]
COVERAGE = ["query_coverage", "subject_coverage"]
SUFFIXES = (".fna", ".fna.gz", ".fa", ".fasta")

#: gdiff presets to tabulate, and the logs each one wrote.  First entry is the chosen preset.
GDIFF_LOGS = [
    ("k25-l1000-n10000-h3", "logs/gdiff-logs/gdiff-sketch-k25-l1000.log",
     "logs/gdiff-logs/gdiff-dist-k25-l1000-n10000-h3.log"),
    ("k23-l666-n10000-h3", "logs/gdiff-logs/gdiff-sketch-k23-l666.log",
     "logs/gdiff-logs/gdiff-dist-k23-l666-n10000-h3.log"),
]
GDIFF_CONFIGS = [config for config, _, _ in GDIFF_LOGS]

SKANI_CONFIGS, DASHING_CONFIGS = ["default", "slow"], ["default", "v4"]
DASHING_INDEX = {"default": "jaccard", "v4": "containment"}

RESOURCE_JOBS = [
    ("mash", "k21-s10000", "sketch", "logs/mash-logs/mash-sketch-k21-s10000.log"),
    ("mash", "k21-s10000", "dist", "logs/mash-logs/mash-dist-k21-s10000.log"),
    ("skani", "slow", "triangle", "logs/skani-logs/skani-triangle-slow.log"),
    ("skani", "default", "triangle", "logs/skani-logs/skani-triangle-default.log"),
    ("dashing2", "default", "dist", "logs/dashing-logs/dashing-default.log"),
    ("dashing2", "v4", "dist", "logs/dashing-logs/dashing-v4.log"),
    *[("gdiff", c, stage, log)
      for c, sketch, dist in GDIFF_LOGS
      for stage, log in (("sketch", sketch), ("dist", dist))],
]

#: ANIm quotes first-contig headers, not assembly accessions.
ACCESSION_TO_ID = dict(zip(
    ("JBHLJH010000001.1 CM060809.1 JAKRZM010000001.1 JBCHJW010000001.1 JBCGHR010000001.1 "
     "JBTWUK020000026.1 JBTWUL010000030.1 JAATLQ010000001.1 JBHLJN010000001.1 JAATLR010000001.1 "
     "JAATLT010000001.1 CM167827.1 JBTWUJ010000013.1 JBHLJJ010000001.1 OU860780.1 "
     "JBHLJO010000001.1 JBHLJG010000001.1 JBDJHW010000001.1 FATV01000001.1").split(),
    ("Heliconius_besckei.GCA_044589335.1 Heliconius_charithonia.GCA_030704555.1 "
     "Heliconius_cydno.GCA_036379175.1 Heliconius_demeter.GCA_040967235.1 "
     "Heliconius_doris.GCA_041107715.1 Heliconius_erato_1.GCA_056362945.2 "
     "Heliconius_erato_2.GCA_056362835.1 Heliconius_erato_x_himera-F1.GCA_013421005.1 "
     "Heliconius_hecalesia.GCA_044589225.1 Heliconius_himera_1.GCA_013421025.1 "
     "Heliconius_himera_2.GCA_013421065.1 Heliconius_melpomene.GCA_056956725.1 "
     "Heliconius_numata.GCA_056362935.1 Heliconius_pardalinus.GCA_044589375.1 "
     "Heliconius_sara.GCA_917862395.2 Heliconius_telesiphe.GCA_044589245.1 "
     "Heliconius_timareta_1.GCA_044589445.1 Heliconius_timareta_2.GCA_040966655.1 "
     "Heliconius_timareta_3.GCA_900068125.1").split()))
COL_QUERY_DESC, COL_SUBJECT_DESC, COL_IDENTITY = 3, 5, 6
COL_QUERY_COV, COL_SUBJECT_COV = 7, 8


def short_id(name: str) -> str:
    """``genomes/Foo.GCA_1.1.fna`` -> ``Foo.GCA_1.1`` (idempotent)."""
    base = name.strip().rsplit("/", 1)[-1]
    for suffix in SUFFIXES:
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def rows_of(path: str, skip_header: bool = False):
    """Tab-split rows, dropping blanks and ``#`` comments."""
    skipped = not skip_header
    for line in open(path, encoding="utf-8"):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if not skipped:
            skipped = True
            continue
        yield line.rstrip("\n").split("\t")


def write_table(path: str, header, rows) -> int:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(header) + "\n")
        n = 0
        for row in rows:
            handle.write("\t".join(row) + "\n")
            n += 1
    return n


# --------------------------------------------------------------------------- tables


def build_anim() -> None:
    """Ground truth: distance = 1 - %identity.  Both directions kept, since ANIm is directed."""
    directed = {}
    for f in rows_of(os.path.join(SRC_OUT, "results.tsv"), skip_header=True):
        if len(f) <= COL_SUBJECT_COV:
            continue
        pair = (ACCESSION_TO_ID.get(f[COL_QUERY_DESC].split()[0]),
                ACCESSION_TO_ID.get(f[COL_SUBJECT_DESC].split()[0]))
        if None in pair:
            raise SystemExit(f"unmapped accession in {f[COL_QUERY_DESC][:40]!r}")
        directed[pair] = (float(f[COL_IDENTITY]), float(f[COL_QUERY_COV]), float(f[COL_SUBJECT_COV]))
    rows = []
    for a, b in sorted({tuple(sorted(p)) for p in directed}):
        if (a, b) not in directed or (b, a) not in directed:
            raise SystemExit(f"missing one direction of {a} / {b}")
        identity, cov_a, cov_b = directed[(a, b)]
        rows.append((a, b, f"{1 - identity:.10g}", f"{1 - directed[(b, a)][0]:.10g}",
                     f"{cov_a:.10g}", f"{cov_b:.10g}"))
    rows.sort(key=lambda r: (float(r[2]), r[0], r[1]))
    n = write_table(f"{OUT_TABLES}/pwdist-ANIm-default.tsv",
                    TABLE_HEADER + ["distance_ba", "query_coverage", "subject_coverage"], rows)
    print(f"{OUT_TABLES}/pwdist-ANIm-default.tsv     {n:>5} pairs (ground truth)")


def build_gdiff() -> None:
    for config in GDIFF_CONFIGS:
        rows = [(short_id(f[0]), short_id(f[1]), f[2])
                for f in rows_of(f"{OUT_DIST}/pwdist-gdiff-{config}.tsv", skip_header=True)]
        n = write_table(f"{OUT_TABLES}/pwdist-gdiff-{config}.tsv", TABLE_HEADER, rows)
        print(f"{OUT_TABLES}/pwdist-gdiff-{config:<20} {n:>5} pairs")


def build_mash(config: str = "k21-s10000") -> None:
    """mash emits both directions plus self-comparisons; the merge collapses them."""
    rows = [(short_id(f[0]), short_id(f[1]), f[2])
            for f in rows_of(f"{SRC_OUT}/pwdist-mash-{config}.tsv")]
    n = write_table(f"{OUT_TABLES}/pwdist-mash-{config}.tsv", TABLE_HEADER, rows)
    print(f"{OUT_TABLES}/pwdist-mash-{config}.tsv {n:>5} rows (self + both directions)")


def build_skani(config: str) -> None:
    """Lower-triangular ANI matrix; distance = 1 - ANI / 100."""
    lines = list(rows_of(f"{SRC_OUT}/pwdist-skani-{config}.tsv"))
    names, rows = [], []
    for f in lines[1:]:  # line 0 is the genome count
        name, values = short_id(f[0]), f[1:]
        if len(values) != len(names):
            raise SystemExit(f"skani/{config}: expected {len(names)} ANI values for {name!r}")
        rows += [(other, name, f"{1 - float(ani) / 100.0:.10g}") for other, ani in zip(names, values)]
        names.append(name)
    rows.sort()
    n = write_table(f"{OUT_TABLES}/pwdist-skani-{config}.tsv", TABLE_HEADER, rows)
    print(f"{OUT_TABLES}/pwdist-skani-{config:<20} {n:>5} pairs")


def build_dashing(config: str) -> None:
    """Upper-triangular similarity with ``-`` on self and the lower triangle.

    Containment *is* the shared k-mer fraction; Jaccard is not -- for equal-size sets
    ``J = s / (2 - s)``, so ``s = 2J / (1 + J)``.  Then ``distance = -ln(s) / k``.
    """
    source = f"{SRC_OUT}/dashing-{config}.tsv"
    index, names, rows, k = DASHING_INDEX[config], [], [], None
    for line in open(source, encoding="utf-8"):
        if line.startswith("#Dashing2Options"):
            k = int(m.group(1)) if (m := re.search(r"\bk:(\d+)", line)) else k
        elif line.startswith("#Sources"):
            names = [short_id(f) for f in line.rstrip("\n").split("\t")[1:]]
        elif not line.startswith("#") and line.strip():
            f = line.rstrip("\n").split("\t")
            i = names.index(short_id(f[0]))
            for j in range(i + 1, len(names)):
                s = float(f[j + 1])
                shared = s if index == "containment" else 2 * s / (1 + s)
                d = 1.0 if shared <= 0 else min(1.0, -math.log(shared) / k)
                rows.append((names[i], names[j], f"{d:.10g}"))
    if not names or k is None:
        raise SystemExit(f"{source}: missing '#Sources' or k-mer size in the dashing header")
    rows.sort()
    n = write_table(f"{OUT_TABLES}/pwdist-dashing-{config}.tsv", TABLE_HEADER, rows)
    print(f"{OUT_TABLES}/pwdist-dashing-{config:<18} {n:>5} pairs")


def build_fastani() -> None:
    """One file per pair; empty files (below the ANI floor) are simply absent."""
    distances, empty = {}, []
    for path in sorted(glob.glob(f"{SRC_FASTANI}/*.txt")):
        f = next(rows_of(path), None)
        if f is None:
            empty.append(os.path.basename(path))
            continue
        distances.setdefault(tuple(sorted((short_id(f[0]), short_id(f[1])))), []).append(
            1 - float(f[2]) / 100.0)
    rows = [(a, b, f"{sum(v) / len(v):.10g}") for (a, b), v in sorted(distances.items())]
    n = write_table(f"{OUT_TABLES}/pwdist-fastani-default.tsv", TABLE_HEADER, rows)
    print(f"{OUT_TABLES}/pwdist-fastani-default.tsv  {n:>5} pairs")
    if empty:
        print(f"  note: {len(empty)} empty fastani file(s) skipped, e.g. {empty[0]}", file=sys.stderr)


def build_samples() -> None:
    """Strip FASTA suffixes from the gdiff window samples, in place and idempotently."""
    for config in GDIFF_CONFIGS:
        path = f"{OUT_SAMPLES}/samples-gdiff-{config}.tsv"
        out, changed, count = [], False, 0
        for line in open(path, encoding="utf-8"):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != "config":  # header
                a, b = short_id(f[1]), short_id(f[2])
                changed |= (a, b) != (f[1], f[2])
                f[1], f[2] = a, b
                count += 1
            out.append(f)
        if changed:
            write_table(path, out[0], out[1:])
        print(f"{OUT_SAMPLES}/samples-gdiff-{config:<20} {count:>7} window samples")


# ------------------------------------------------------------------------ resources

GNU_WALL = re.compile(r"Elapsed \(wall clock\) time.*?([0-9]+:[0-9]+:[0-9.]+|[0-9]+:[0-9.]+)\s*$",
                      re.MULTILINE)
GNU_RSS = re.compile(r"Maximum resident set size \(kbytes\):\s*([0-9]+)")
TOOL_WALL = re.compile(r"elapsed:\s*([0-9.]+)\s*sec")


def log_metrics(text: str) -> tuple[float | None, int | None]:
    """Wall seconds and peak RSS.  ``time -v`` wall is preferred when the run carried it."""
    wall = None
    if m := GNU_WALL.search(text):
        wall = 0.0
        for part in m.group(1).split(":"):
            wall = wall * 60 + float(part)
    elif m := TOOL_WALL.search(text):
        wall = float(m.group(1))
    rss = GNU_RSS.search(text)
    return wall, int(rss.group(1)) if rss else None


def build_resources() -> None:
    """``tables/resources.tsv``: wall summed over stages, RSS the maximum over them."""
    rows, missing = [], []
    for (tool, config), jobs in groupby(RESOURCE_JOBS, key=lambda j: (j[0], j[1])):
        stages, walls, rsses = list(jobs), [], []
        for _, _, stage, relative in stages:
            name = f"{tool}/{config}/{stage}"
            if os.path.exists(relative):
                wall, rss = log_metrics(open(relative, encoding="utf-8", errors="replace").read())
            else:
                wall, rss, name = None, None, f"{name} ({relative} missing)"
            rows.append([tool, config, stage, "" if wall is None else f"{wall:.6g}",
                         "" if rss is None else str(rss), relative])
            if wall is None:
                missing.append(name)
            else:
                walls.append(wall)
            if rss is not None:
                rsses.append(rss)
        if len(stages) > 1:
            rows.append([tool, config, "total", "" if not walls else f"{sum(walls):.6g}",
                         "" if not rsses else str(max(rsses)), "sum of stages"])
    # fastani ran once per pair: wall summed (the real cost), RSS the high-water mark.
    logs = sorted(glob.glob("logs/fastani-logs/*.log"))
    if logs:
        metrics = [log_metrics(open(p, encoding="utf-8", errors="replace").read()) for p in logs]
        walls = [w for w, _ in metrics if w is not None]
        rsses = [r for _, r in metrics if r is not None]
        rows.append(["fastani", "default", "dist", "" if not walls else f"{sum(walls):.6g}",
                     "" if not rsses else str(max(rsses)), f"{len(logs)} per-pair runs, summed"])
    n = write_table(f"{OUT_TABLES}/resources.tsv",
                    ["tool", "config", "stage", "wall_seconds", "max_rss_kb", "source"], rows)
    print(f"{OUT_TABLES}/resources.tsv{'':<27} {n:>5} rows (wall time + peak RSS)")
    if missing:
        print(f"  note: no wall time in {len(missing)} log(s): {', '.join(missing)}", file=sys.stderr)


# ---------------------------------------------------------------------------- merge

TABLE_RE = re.compile(r"^pwdist-(?P<method>[A-Za-z0-9]+)-(?P<config>.+)\.tsv$")
METHOD_LABELS = {"anim": "ANIm", "fastani": "FastANI", "dashing": "dashing2"}


def read_table(path: str):
    """``pair -> mean estimate`` and, where present, ``pair -> coverage``."""
    reader = csv.DictReader(open(path, encoding="utf-8"), delimiter="\t")
    has_coverage = set(COVERAGE) <= set(reader.fieldnames or [])
    records, coverage = {}, {}
    for row in reader:
        if row["genome_a"] == row["genome_b"]:
            continue
        pair = "|".join(sorted((row["genome_a"], row["genome_b"])))
        records.setdefault(pair, []).append(float(row["distance"]))
        if has_coverage:
            coverage[pair] = {f: row[f] for f in COVERAGE}
    return {p: sum(v) / len(v) for p, v in records.items()}, coverage


def pearson(xs, ys) -> float | None:
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sxx**0.5 * syy**0.5)


def merge(table_dir: str, truth_token: str, out: str) -> int:
    """Stack every ``tables/pwdist-*.tsv`` into one long table carrying the ANIm truth."""
    loaded, truth_path, reference, truth_coverage = [], None, {}, {}
    for name in sorted(os.listdir(table_dir)):
        if match := TABLE_RE.match(name):
            method = match.group("method").lower()
            path = os.path.join(table_dir, name)
            estimates, coverage = read_table(path)
            loaded.append((method, METHOD_LABELS.get(method, method), match.group("config"), estimates))
            if method == truth_token:
                truth_path, reference, truth_coverage = path, estimates, coverage
    if not loaded:
        raise SystemExit(f"no pwdist-<method>-<config>.tsv files in {table_dir!r}")
    if truth_path is None:
        raise SystemExit(f"no table for the ground truth {truth_token!r} in {table_dir!r}")

    rows, inverted = [], []
    for method, label, config, estimates in loaded:
        for pair, value in estimates.items():
            a, b = pair.split("|")
            covered = truth_coverage.get(pair, {})
            rows.append([a, b, "" if pair not in reference else f"{reference[pair]:.10g}",
                         covered.get("query_coverage", ""), covered.get("subject_coverage", ""),
                         label, config, f"{value:.10g}"])
        if method != truth_token:
            shared = sorted(set(estimates) & set(reference))
            r = pearson([reference[p] for p in shared], [estimates[p] for p in shared])
            if r is not None and r < 0:  # a similarity would invert every downstream comparison
                inverted.append(f"{label}/{config}")
        print(f"  {f'{label}/{config}':<28} {len(estimates):>6} pairs", file=sys.stderr)
    if inverted:
        raise SystemExit("error: these hold a similarity (higher = closer), not a distance: "
                         + ", ".join(inverted))
    rows.sort()
    return write_table(out, ["genome_a", "genome_b", METHOD_LABELS.get(truth_token, truth_token)]
                       + COVERAGE + ["method", "config", "estimate"], rows)


# ----------------------------------------------------------------------------- main


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--table-dir", default=OUT_TABLES)
    parser.add_argument("--truth", default="anim", help="ground-truth method token [anim]")
    parser.add_argument("-o", "--out", default=MERGED)
    args = parser.parse_args()

    build_anim()
    build_fastani()
    build_gdiff()
    build_mash()
    for config in SKANI_CONFIGS:
        build_skani(config)
    for config in DASHING_CONFIGS:
        build_dashing(config)
    build_samples()
    build_resources()

    n = merge(args.table_dir, args.truth.strip().lower(), args.out)
    print(f"wrote {args.out} ({n} rows)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
