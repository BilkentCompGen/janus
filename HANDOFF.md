# Handoff — 2026-08-28 evening

## What changed in the code

`src/graph.cpp` only. Uncommitted (nothing was committed; CLAUDE.md says not to without being asked).
Baseline copies of the file at each step are under `tmp/prof/`.

1. **BUG FIX — reverse-strand cores stored the wrong bytes.** `extract2` returns a reverse-strand
   core's `[st, en)` in the reverse complement's own frame; the code sliced the FORWARD read at those
   offsets, so ~half of every graph's stored core sequence was a real but unrelated span from the
   mirror position. Proven with `tmp/prof/t_extract.cpp` (0/430 reverse cores hashed back to their
   own signature, against 430/430 for the mirrored interval) and on real data (50.5% of stored cores
   reproduced their signature before the fix, 100.0% after). Lengths were never affected; what was
   wrong is `bubbles`' minimizers, hence `drop.sigs`.
2. **LUT-driven packing.** ~21% off the `graph` stage. Details and the negative result about the
   removed memcpy are in `docs/NOTES.md`.
3. **Profiling instrumentation**, additive: a second per-read CPU counter, a sampled timer on the
   packing block (with the clock's own cost calibrated and subtracted), exact per-structure byte
   accounting per sample, a 200 ms RSS trace written to `<dir>/rss_trace.tsv`, and a signature count
   histogram. Prints under `[prof]` and `[bytes]` beside the existing `[eff]` table.

Verification protocol used for every change, and the one to keep using: run the unmodified binary
twice on the same input, compare `graph.bin` section by section, and require the changed binary to
sit in the same equivalence class — `sig` / `cov` / `E` / `esrc` / `edst` byte-identical, `slen` and
`seq` differing no more than two baseline runs differ from each other (~1.15% of `slen` bytes).

## What it cost in accuracy — READ THIS BEFORE TRUSTING ANY OLD NUMBER

| pair | before (RESULTS.md) | after the fix | ΔF1 |
|---|---|---|---|
| h2009 | 0.854 / 0.876 / 0.865 | 0.850 / 0.872 / 0.861 | −0.004 |
| h1437 | 0.727 / 0.783 / 0.754 | 0.665 / 0.783 / 0.719 | −0.035 |

Every accuracy number recorded before 2026-08-28 was measured with the defect. `bubbles`'
thresholds were calibrated against broken minimizers, and `max_unmatched` has been ruled out as the
knob that recovers this — the sweep is in `docs/NOTES.md`. **This is a decision, not a bug: left for
yuio.**

## What is on disk

- `tmp/measure/h2009/snape_out/` and `tmp/measure/h1437/snape_out/` — complete pipelines with the
  fixed binary, every stage timed under `/usr/bin/time -v` into `results/measure/<pair>/`.
- `results/measure/<pair>/minda/` — scored against CASTLE `minda_4_2`.
- `tmp/night/clean/{H2009,H1437,HCC1954,HCC1937}_minda_4_2.vcf` — truth sets with the archive's
  malformed `cmd:` line removed. Record counts match `docs/RESULTS.md`: 925 / 221 / 1106 / 703.
- `tmp/night/castle_archive/` (13 GB) — all four pairs' caller VCFs, intact.
- Reads: `hifi/h2009/{n,t}.fastq` (820 GB) and `hifi/h1437/{n,t}.svb` (98 GB). hcc1954 and hcc1937
  have no reads on disk; re-measuring those two needs a re-download.

## Running when this was written

- `tmp/night/fetch_h1437.sh` — resuming the h1437 `.gz` download (~165 GB left), log
  `tmp/night/fetch_h1437.log`. It keeps the `.gz` because the timing arm must be measured from plain
  FASTQ, not `.svb`.
- `tmp/night/exp2_fixed.sh` — experiment 2 conditions A and C on both pairs with the fixed calls,
  into `results/exp2_fixed/`.

## Next

1. Decide what to do about the `bubbles` calibration. Nothing else should be measured against the
   paper until that is settled.
2. h2009's alignment-first arm is still unmeasured: `mm2_n`, severus, sniffles, cuteSV, dysgu. Those
   need a quiet machine, so they must wait for the download to finish.
3. `tools/` (longcallD, sawfish, pbsv, SVDSS) was lost from the working tree in the `.git` cleanup
   and is untracked-empty; re-fetch from upstream releases if those callers are wanted.
4. `docs/RESULTS.md`'s cost table is stale in at least the `bubbles` row (04:02 there, 7–9 s
   measured on two pairs today). The new per-stage table is in `docs/NOTES.md`.
