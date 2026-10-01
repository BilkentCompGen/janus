# bench

Everything that runs snape on the real tumor/normal pairs, measures it, and scores it. The reads are in `hifi/<pair>/{n,t}.fastq.gz`
(`hifi/fetch.sh`), the truth sets in `hifi/truth/`.

## Scripts

| script | what it does | writes |
|---|---|---|
| `run.sh <experiment> [pair ...]` | runs the whole pipeline (`pipeline.sh`, 16 threads, cores 0-15), one run after the other, profiled with external tools. Experiments: `full` (all reads of each pair) and `coverage` (each pair at 5x, 10x, ... 60x, as deep as both samples go; 1x = 3.1 Gb; pack keeps the first N x 3.1 Gb of each sample, `pack -b`, and the run continues from those svbs; depth by depth over the pairs). Intermediates are measured, then moved to `~/trash` (full) or deleted (coverage). | `runs/<experiment>/...` |
| `compete.sh [pair ...]` | the competitors on the same reads: per pair, align (minimap2 + samtools sort, measured), every tool in `callers/` on the full BAMs, the core tools (Severus, SAVANA, nanomonsv, Sniffles2, SVision-pro) on 5x/10x/20x/40x subsets of exactly snape's coverage reads; cores 0-15, 16 threads; docker tools measured through a cgroup slice per run; BAMs and intermediates deleted after the pair | `runs/align/<pair>/{n,t}/`, `runs/competitors/<tool>/<pair>[_<N>x]/` |
| `measure_tool.py <run folder>` | one competitor or alignment run -> run.json (wall, CPU = host + containers, peak memory, calls, file sizes), files.tsv (called by `compete.sh`) | `run.json`, `files.tsv` |
| `callers/` | one wrapper per competitor (`<tool>.sh <tumor.bam> <normal.bam> <outdir> <threads>`), `COMPETITORS.md` (the literature), `README.md` (images, commands, status) | |
| `measure.py <run folder>` | turns one run into structured measurements (called by `run.sh`) | `run.json`, `steps.tsv`, `files.tsv` |
| `verdict.py [seed]` | held-out verdict: boosted trees (a learned component) fitted on the other pairs' calls, applied to the held-out pair | `verdict/<pair>.vcf` |
| `eval.sh <pair> <vcf> <name>` | PASS records scored with minda (CASTLE consensus rule, or NIST truthset mode for HG008) | `eval/<pair>/<name>/` |
| `score.py <truth.vcf> <calls.vcf>` | the event matcher (500 bp at both ends, insertions by position) used by verdict.py and the figures | stdout |

Run an experiment detached: `tmux new-session -d -s queue bash bench/run.sh coverage`. Its progress is in
`runs/<experiment>/queue.log`.

## Folders

```
runs/full/<pair>/          one run of snape on all reads of a pair, from fastq.gz to calls
runs/coverage/<pair>_<N>x/ the same at N x coverage (N = 5, 10, ... 60); same files; monitor/params.json has the
                      experiment, coverage and base limit, run.json "params" too
  calls.vcf           the calls (all the verdict and evaluation need)
  run.json            every measurement in one file: pair, commit, machine, reads (bytes, reads, bases), totals,
                      per step (wall, peak RSS, user/sys/CPU s, drive read/written GB, drive utilization, max drive /
                      controller / CPU temperature, snape's own report lines), and the size of every file the run made
  steps.tsv           the per-step part of run.json as a table, last row = total
  files.tsv           every file the run made, its size, and whether it is kept, went to trash or was deleted
  logs/               each step's own output
  monitor/            raw inputs of run.json: pidstat (2 s), iostat, vmstat, temperatures (5 s), start.txt and
                      total.txt (/usr/bin/time of the whole run), note.txt when the run was made differently from
                      run.sh (also in run.json as "note")
runs/<experiment>/queue.log  the queue's progress
runs/align/<pair>/{n,t}/      read alignment for the competitors: run.json (time, CPU, memory, BAM bytes)
runs/competitors/<tool>/<pair>[_<N>x]/   one competitor run: calls.vcf, log.txt, run.json, files.tsv
verdict/              the verdict of record: <pair>.vcf (PASS / LOW / DUPL), log.txt (cut and PASS per pair),
                      time.txt / times.txt (its cost), commit.txt (the code commit it ran on)
eval/<pair>/<name>/   minda outputs: snape = the verdict's PASS records; severus, sniffles = other callers on HG008
callers/              other callers run by us from read alignments (Severus, Sniffles2): run.sh, times.txt, repro.txt
```

`steps.tsv` units: seconds, KiB, GB (10^9 bytes), degrees C. CPU s = user + sys, measured by `/usr/bin/time`.
