# CLAUDE.md

This file guides Claude Code when working in this repository.

## Project

`svp` is a somatic structural-variant caller. The project is being rebuilt from scratch on the **`vg-rewrite`** branch around a **variation graph** model. `main` currently holds only the bare scaffolding; everything below describes where prior work lives.

Until the new implementation lands, do **not** try to reconstruct the old pipeline on this branch. Read the relevant prior branch instead.

## Branch map

| Branch | What's there |
|---|---|
| `main` | clean slate — this branch |
| `vg-rewrite` | active development for the variation-graph rewrite |
| `backup` | previous monolithic LCP-graph implementation (last known-good `svp` binary, full `src/`, tests, scripts, docs) |
| `experimental/modular-pipeline` | `svp2/` — 6-stage modular pipeline (s1_lcp → s6_finalize). Authoritative entry: `svp2/STATUS.md` on that branch |
| `kmer` | k-mer-based experiment |
| `no-compression` | variant of the LCP pipeline without depth-deepening |

To inspect prior code, `git checkout <branch>` — do not try to port files into `main` without an explicit instruction.

## What stays on `main`

- `docs/GRAPH_MODELS.md` — conceptual reference on de Bruijn graphs and variation graphs (formal definitions, variant families, GFA, ecosystem). Background reading for the rewrite.
- `scripts/truvari.sh` — Truvari benchmark wrapper (kept; new pipeline will still emit VCF and benchmark against the same simulation truth).
- `simulation/` (gitignored) — insilicoSV truth data; configs are regenerated as needed on whatever branch needs them.
- `lib/` (gitignored) — bundled dependencies. The rewrite is expected to keep `lib/lcptools` (LCP cores) and `lib/insilicoSV` (truth-data generator).

## Dependencies (carried over)

- LCP tools (`lib/lcptools`) — locally consistent parsing core extraction
- insilicoSV (`lib/insilicoSV`) — somatic / germline simulation generator
- truvari — SV-level VCF benchmarking (`pip install truvari`, or via the project venv)

Anything else (parasail, htslib, zstr, the old C++ source tree) is **not** on `main` — pull it from `backup` if it's needed.

## When the rewrite begins

Build system, source layout, and run commands will be defined on `vg-rewrite` and documented here once stable. Until then, this file intentionally contains no build/run instructions for the new pipeline.
