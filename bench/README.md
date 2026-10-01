# bench

Everything that runs snape on the real tumor/normal pairs, measures it, and scores it. The reads are in `hifi/<pair>/{n,t}.fastq.gz`
(`hifi/fetch.sh`), the truth sets in `hifi/truth/`.

## Scripts

| script | what it does | writes |
|---|---|---|
| `run.sh [pair ...]` | runs the whole pipeline (`pipeline.sh`, 16 threads, cores 0-15) on each pair, one pair after the other. It profiles each run with external tools and moves the intermediates to `~/trash` after measuring them. It waits while less than 400 GB is free. | `runs/<pair>/` |
| `measure.py <runs/pair>` | turns one run into structured measurements (called by `run.sh`) | `run.json`, `steps.tsv`, `files.tsv` |
| `verdict.py [seed]` | held-out verdict: boosted trees (a learned component) fitted on the other pairs' calls, applied to the held-out pair | `verdict/<pair>.vcf` |
| `eval.sh <pair> <vcf> <name>` | PASS records scored with minda (CASTLE consensus rule, or NIST truthset mode for HG008) | `eval/<pair>/<name>/` |
| `score.py <truth.vcf> <calls.vcf>` | the event matcher (500 bp at both ends, insertions by position) used by verdict.py and the figures | stdout |

Run the queue detached: `tmux new-session -d -s queue bash bench/run.sh`. Its progress is in `runs/queue.log`.

## Folders

```
runs/<pair>/          one run of snape, from fastq.gz to calls
  calls.vcf           the calls (all the verdict and evaluation need)
  run.json            every measurement in one file: pair, commit, machine, reads (bytes, reads, bases), totals,
                      per step (wall, peak RSS, user/sys/CPU s, drive read/written GB, drive utilization, max drive /
                      controller / CPU temperature, snape's own report lines), and the size of every file the run made
  steps.tsv           the per-step part of run.json as a table, last row = total
  files.tsv           every file the run made, its size, and whether it is kept or went to trash (svbs, counts, graph, ...)
  logs/               each step's own output
  monitor/            raw inputs of run.json: pidstat (2 s), iostat, vmstat, temperatures (5 s), start.txt and
                      total.txt (/usr/bin/time of the whole run), note.txt when the run was made differently from
                      run.sh (hcc1954: steps paused 30 s apart; also in run.json as "note")
runs/queue.log        the queue's progress
verdict/              the verdict of record: <pair>.vcf (PASS / LOW / DUPL), log.txt (cut and PASS per pair),
                      time.txt / times.txt (its cost), commit.txt (the code commit it ran on)
eval/<pair>/<name>/   minda outputs: snape = the verdict's PASS records; severus, sniffles = other callers on HG008
callers/              other callers run by us from read alignments (Severus, Sniffles2): run.sh, times.txt, repro.txt
```

`steps.tsv` units: seconds, KiB, GB (10^9 bytes), degrees C. CPU s = user + sys, measured by `/usr/bin/time`.
