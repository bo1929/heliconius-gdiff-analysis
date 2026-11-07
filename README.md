# Heliconius analysis using gdiff
```
ani-analysis/          LEVEL 1 — ANI comparison + dip test
  build.py                 builds tables/ + merged-distances.tsv from outputs/
  plot-butterfly-ani.R     reads samples/ + dist/
  outputs/                 per-tool outputs: ANIm raw, mash, skani, dashing2, fastani
  tables/                  pwdist-<method>-<cfg>.tsv      (built by build.py)
  samples/                 samples-gdiff-<cfg>.tsv        (per-window distances)
  dist/                    pwdist-gdiff-<cfg>.tsv         (per-pair distances)
  logs/                    gdiff-{sketch,dist}-<cfg>.log + other tools' logs
  figures/                 *.pdf
  merged-distances.tsv     every table stacked with the ANIm reference
genomes.list, genomes.tsv   manifests — the sequences themselves live on the remote
```
