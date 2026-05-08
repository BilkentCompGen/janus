# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`svp` is a somatic structural-variant caller. Active development is on **`vg-rewrite`** (this branch). It is a pairwise normal-vs-tumor caller anchored on **LCP cores** — no graph construction, no rGFA, single binary. Output is a flat bubble TSV; mapping back to per-contig coordinates / VCF is a downstream concern.

## Branch map

| Branch | What's there |
|---|---|
| `main` | clean slate, scaffolding only |
| `vg-rewrite` | this branch — active LCP-anchor rewrite |
| `backup` | previous monolithic LCP-graph impl. Source of `zstr.hpp`, `MurmurHash3.{h,cpp}` |
| `experimental/modular-pipeline` | `svp2/` 6-stage modular pipeline; entry: `svp2/STATUS.md` |
| `kmer`, `no-compression` | dead-end experiments |

Don't port files across branches without explicit instruction.

## Build & run

```
make install-deps       # one-time: builds lib/lcptools (liblcptools.a)
make                    # builds ./svp
make clean              # moves build/, svp, outputs/ → ~/trash/svp-make-clean/
make run-<id>           # ids: del dup ins inv mixed snp
                        # auto-generates simulation/truth/<id>/ if missing,
                        # writes outputs/svp/<id>/bubbles.tsv, /usr/bin/time -v
make run-list           # shows ids and overrides
```

Overrides on `make run-<id>`:
- `LEVEL=N` lcptools deepening level (default 5)
- `FLANK=M` ±M bp window for core signature (default 10)
- `NORMAL=path/to/normal.fa` defaults to `simulation/ref.fasta`. Use a germline truth FASTA when calling against simulated germline+somatic.

Direct invocation: `./svp <normal.fa[.gz]> <tumor.fa[.gz]> <out.tsv> [--level N] [--flank M]`.

## Pipeline

`main.cpp` orchestrates this exact sequence — one stage per module, no cross-module state:

1. **`read_fasta_flat`** (`io.cpp`) — zlib `gzopen`/`gzgets`. Concatenates **all records** of one FASTA into a single `std::string`, separated by `2*flank+1` `N`-bumpers (so neither LCP cores nor signature flanks bridge records). Emits `rec_start[]`, `rec_end[]` for downstream per-record processing. Record names are dropped.
2. **`extract_cores`** (`cores.cpp`) — calls `init_lps_offset` + `lps_deepen` from `lib/lcptools`, **per record** (because lcptools takes `int len` and full-genome FASTAs > INT_MAX). Cores carry the lcptools `label` for use as a hash seed.
3. **`compute_sigs`** (`cores.cpp`) — `signature(core) = MurmurHash64(seq[start-m .. end+m], seed=core.label)`. Single 64-bit key collapses both label-equality and flank-equality checks into one `unordered_map` lookup.
4. **`unique_anchors`** (`match.cpp`) — keeps cores whose signature appears exactly once on **both** sides. Two `unordered_map<uint64_t,int32_t>` first-sighting passes; `-1` marks "seen ≥2 times". Output anchors are pre-sorted by `n_idx` because we walk normal first.
5. **`chain`** (`match.cpp`) — patience-sort LIS on tumor index, O(K log K), with `parent[]` reconstruction. Drops out-of-order anchors (translocations, inversions). Verified on the user's `ABCDE`/`ADBCE` example.
6. **`emit_bubbles`** (`match.cpp`) — between consecutive chained anchors (with virtual head/tail), greedy byte-equality `ext_r` then `ext_l` to expand matching flanks. Emits `Bubble{n_start, n_end, t_start, t_end, n_seq, t_seq, type}`. Type is classified at construction (`SNP`/`INS`/`DEL`/`CPLX`) from `(nl, tl)` lengths.
7. **`write_bubbles`** (`io.cpp`) — TSV with header `#type\tn_start\tn_end\tt_start\tt_end\tn_seq\tt_seq`. Empty seqs print as `.`.

Shared types are POD-only in `src/svp.h` (`Core`, `Anchor`, `Bubble`). No constructors, no methods, no inheritance.

## Coordinates

- 0-based half-open everywhere (matches lcptools).
- Bubble positions are **flat offsets into the concatenated buffer**, not into any FASTA record. Translating back is a separate downstream tool's job — that's why `read_fasta_flat` exposes `rec_start[]`/`rec_end[]`.
- N-bumper width is `2*flank+1` exactly so LCP cores and signature flanks cannot span records.

## Code style — strict

Read `docs/CODE_STYLE.md` before editing. Hard rules: C++17, no OOP, no user templates, no try/catch (fail-fast with `[FATAL]` + `std::abort()`), avoid if-else trees, ≤30 LOC per function, ≤700 LOC total across `src/*.cpp`. Comments only when the WHY is non-obvious. Single-line function signatures (`.clang-format` has `ColumnLimit: 0` — don't reformat).

`.o` and `.d` artifacts go to `build/`, never `src/`.

## Simulation tooling

`simulation/Makefile` is independent (delegated by the top-level Makefile). Layout is **flat** — no per-chromosome subdirs:

```
simulation/
  ref.fasta              base genome (gitignored)
  ref.fasta.fai
  beds/{chr19,genome}.bed   yaml's overlap_regions decides scope
  configs/<id>.yaml      insilicoSV config (id == filename stem)
  truth/<id>/            generated: out.fa, out.fa.fai, out.vcf, [novel.fa]
```

Common targets: `make -C simulation <id>` (or `all`), `make -C simulation list`, `make -C simulation index`, `make -C simulation clean-<id>`. Whether a sim is chr19-only or whole-genome depends solely on which BED its yaml points at — there are no separate folders per scope.

Truvari benchmark: `make -C simulation bench/<id> CALL=path/to/detected.vcf` (uses `scripts/truvari.sh`). Note: svp's flat-buffer bubble TSV is **not** a VCF — Truvari benchmarking requires a downstream coordinate-mapper that does not yet exist on this branch.

## Dependencies

- `lib/lcptools` — locally consistent parsing core extraction. Build with `make install-deps`. The static lib is stale-prone; if you see segfaults on struct fields after a `git pull` of `lib/`, rebuild it.
- `lib/insilicoSV` — truth simulator (driven by `simulation/Makefile`).
- `zlib` — for `gzopen`/`gzgets` FASTA reading.
- `truvari` — for VCF-level benchmarking (when a VCF emitter exists).

`lib/`, `simulation/`, `outputs/` are all gitignored.

## Current accuracy focus

Per the user, only **SNP / INS / DEL** classification matters for the moment. INV / DUP / CPLX may be miscalled — don't sink time into those without explicit direction. Single-bp position-equivalence ambiguity (homopolymer left/right alignment) within ±5 bp is acceptable.
