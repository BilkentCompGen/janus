# bench

Everything that runs snape on the real tumor/normal pairs, measures it, and scores it. The reads are in `hifi/<pair>/{n,t}.fastq.gz`
(`hifi/fetch.sh`), the truth sets in `hifi/truth/`.

## Scripts

| script | what it does | writes |
|---|---|---|
| `all.sh` | every benchmark of the final snape (tag snape-final), one after the other: the full runs, the verdict, minda's scoring, the coverage series, the thread scaling, the competitors; each stage resumes where it stopped | `runs/all.log`, `runs/minda.txt` |
| `run.sh <experiment> [pair ...]` | runs the whole pipeline (`pipeline.sh`, 16 threads, cores 0-15), one run after the other, profiled with external tools. Experiments: `full` (all reads of each pair), `coverage` (each pair at 1x, 3x, 5x, 10x, 15x, ... 60x, as deep as both samples go; 1x = 3.1 Gb; pack keeps the first N x 3.1 Gb of each sample; up to 10x the support threshold is 2 reads instead of 3), `threads` (HCC1954 on all reads with 1, 2, 4, 8, 16 and 32 threads, pinned to as many cores). Intermediates are measured, then deleted. | `runs/<experiment>/...` |
| `compete.sh [pair ...]` | the competitors on the same reads: per pair, align (minimap2 + samtools sort, measured), every tool in `callers/` on the full BAMs, the core tools (Severus, SAVANA, nanomonsv, Sniffles2, SVision-pro) on 5x/10x/20x/40x subsets of exactly snape's coverage reads; cores 0-15, 16 threads; docker tools measured through a cgroup slice per run; BAMs and intermediates deleted after the pair | `runs/align/<pair>/{n,t}/`, `runs/competitors/<tool>/<pair>[_<N>x]/` |
| `measure_tool.py <run folder>` | one competitor or alignment run -> run.json (wall, CPU = host + containers, peak memory, calls, file sizes), files.tsv (called by `compete.sh`) | `run.json`, `files.tsv` |
| `callers/` | one wrapper per competitor (`<tool>.sh <tumor.bam> <normal.bam> <outdir> <threads>`), `COMPETITORS.md` (the literature), `README.md` (images, commands, status) | |
| `measure.py <run folder>` | turns one run into structured measurements (called by `run.sh`) | `run.json`, `steps.tsv`, `files.tsv` |
| `verdict.py [seed]` | the verdict, always held out: for each of the six labeled pairs (CASTLE x4, HG008, COLO829) boosted trees (a learned component, 5 seeds averaged) are fitted on the other five pairs' records, the cut is the best F1 on those records, and the held-out pair is scored; the pairs without a truth set (HG009, NCI-H209, NCI-H2126) by trees fitted on all six. HCC1395 is left out (see below). | `verdict/<pair>.vcf`, `verdict/folds.json`, `verdict/log.txt` |
| `eval.sh <pair> <vcf> <name>` | PASS records scored with minda (CASTLE consensus rule, or NIST truthset mode for HG008) | `eval/<pair>/<name>/` |
| `score.py <truth.vcf> <calls.vcf>` | the event matcher (500 bp at both ends, insertions by position) used by verdict.py and the figures | stdout |
| `compare_groups.py` | the decision record of 2026-10-02: snape with and without select's per-distance-group pieces on the CASTLE pairs + HG008, same verdict and scoring (F1 .8651 vs .8654, minda .8785 vs .8783); its run folders were cleared with the old results | `runs/groups/compare/` |
| `truth.py` | writes the truth sets that are not chr-named VCFs (HCC1395 from the SEQC2 table, COLO829 with chr names) in one form | `hifi/truth/<pair>/truth.vcf` |
| `dumpd.sh <pair> <dir>`, `groups_study.py <d.tsv>` | diagnostics of select's distance groups: the pipeline up to select with every measured anchor distance dumped, and how the distances of a site spread | `<dir>/d.tsv`, stdout |

Run everything detached: `tmux new-session -d -s all bash bench/all.sh` (progress in `runs/all.log`), or one experiment with `bash bench/run.sh <experiment>`; its progress is in
`runs/<experiment>/queue.log`.

## Folders

```
runs/full/<pair>/          one run of snape on all reads of a pair, from fastq.gz to calls
runs/coverage/<pair>_<N>x/ the same at N x coverage (N = 1, 3, 5, 10, 15, ... 60); same files; monitor/params.json has
                      the experiment, coverage, base limit, support reads and threads, run.json "params" too
runs/threads/<pair>_<N>t/  the same on all reads with N threads
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
runs/all.log          all.sh's progress; runs/minda.txt minda's output for the verdict's PASS records
runs/align/<pair>/{n,t}/      read alignment for the competitors: run.json (time, CPU, memory, BAM bytes)
runs/competitors/<tool>/<pair>[_<N>x]/   one competitor run: calls.vcf, log.txt, run.json, files.tsv
verdict/              the verdict: <pair>.vcf (PASS / LOW / DUPL, INFO VP = the held-out probability), log.txt (cut and
                      PASS per pair), folds.json (per pair: training pairs, cut, training and inference time, every
                      record's probability and label), stdout.txt
eval/<pair>/snape/    minda's scoring of the verdict's PASS records (CASTLE pairs: consensus rule; HG008: NIST truthset mode)
callers/              the competitor wrappers (callers/README.md)
```

`steps.tsv` units: seconds, KiB, GB (10^9 bytes), degrees C. CPU s = user + sys, measured by `/usr/bin/time`.

## Pairs and truth sets

Ten pairs are run: HCC1954, NCI-H2009, HCC1937, NCI-H1437 (truth: the CASTLE consensus without snape), HG008 (NIST draft
somatic SV benchmark), COLO829 (its published truth set), HCC1395, HG009, NCI-H209, NCI-H2126. The six with a truth set
used by the verdict are the four CASTLE pairs, HG008 and COLO829; HG009, NCI-H209 and NCI-H2126 have none. HCC1395 is run
and measured but left out of the verdict and the accuracy: its truth set (SEQC2, Talsania et al. 2022) disagrees with
every long-read caller (Feng et al. 2025).
