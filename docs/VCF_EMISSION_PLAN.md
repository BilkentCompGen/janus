# Plan — `svp posmap → anchor → vcf` + DUP/INV + truvari benchmark

## Context

`svp` is a somatic SV caller on branch `graph-inspector`. As of commit `f68e8f1` ("modeb is implemented"), the pipeline produces `outputs/svp/wgs_mm10/modeb.bin` — per-bubble `BubbleTag` records (s1, s2, tag ∈ {E,D,I,S,C}, bucket, pre, suf, N'/T' residual bytes, N'/T' coverage RLE). Everything to date is **coordinate-free**. The next stage maps bubbles to (CHROM, POS), emits a VCF, layers in DUP + INV detection, and benchmarks against `simulation/truth/wgs/out.vcf` via truvari.

**Targets**: ≥90% recall AND ≥90% precision on SNP / INS / DEL / DUP at the ≥50bp truvari threshold; INV may fall short (accepted). Also a sub-50bp run for small indels and SNPs.

### Measured baselines (2026-06-05)

Mode B output at `outputs/svp/wgs_mm10/modeb.bin` (mm=10):

| | total | ≥50bp residual | <50bp residual |
|---|---|---|---|
| S (substitution) | 32,812 | 118 | 32,694 |
| D (deletion)     | 2,267  | 1,352 | 915 |
| I (insertion)    | 2,257  | 1,343 | 914 |
| C (complex)      | 318    | 265 | 53 |
| **total**        | **37,654** | **3,078** | **34,576** |
| coverage filter passes | — | 2,813 (only C drops) | all DIS pass |

Truth at `simulation/truth/wgs/out.vcf`: 16,950 SNP + 470 small DEL + 50 SV DEL + 470 small INS + 15 SV INS + 32 DUP + 15 INV = **18,002** total. 112 are ≥50bp (the truvari primary target).

Best prior truvari result is `outputs/truvari/v24/`: 88.4% recall, 1.4e-5 precision (7.2M FP) — recall is close to target, **precision is the gap**. With Mode B's 37,654 bubbles (vs v24's 7.2M emissions), the FP pool is 200× smaller from day one.

### Hard constraints (mistakes to avoid)

- **INV via single-bubble byte RC fails**: measured 0/250 S-tag bubbles with |N'|=|T'|>1 satisfy N'==RC(T'). Real INVs are fragmented per `docs/INSPECTOR_FINDINGS.md`. INV detection MUST be cross-bubble or accept very low recall.
- **Coverage filtering alone is weak**: D/I/S pass at 100% under cov_t≥5 / cov_n≤2 — only C drops. **Positional filters** are the lever for FP reduction.
- **DUP via cycle SCC at WG scale infeasible**: 1 giant SCC swallows 97% of super-nodes. The 105 small SCCs (size 2-2280) are useful but won't catch tandem DUPs whose anchors lie in the giant SCC. Use ref-substring + coverage instead.
- **5 v1 INV approaches are documented dead ends** (per memory `project_caller_v1_fix_attempts`): per-bubble RC equality, posmap strand=1, per-bubble RC substring, global k-mer RC, union-find super-bubble. DO NOT retry.
- **DUP via anchor-overlap geometry undersizes**: don't reuse that. Size must come from arm splice length.

## Architecture: 4 new pipeline stages on top of Mode B

```
svp posmap  <ref.fasta> <dir>         # ~50s WG, ~1.5 GB out
            → dir/sig_ref_map.bin       sig → (chrom_id, pos) oracle
            → dir/chroms.tsv            chrom_id → name table

svp anchor  <dir>                     # ~30s
            → dir/bubble_pos.bin        per-bubble (chrom, pos_start, pos_end, n_resolved)

svp vcf     <dir> [--out <path>]      # ~10s
            → dir/calls.vcf             SNP/INS/DEL/DUP/INV + filters

scripts/svp_benchmark.sh <dir>        # ~2 min
            → dir/truvari_geq50/        ≥50bp SV benchmark
            → dir/truvari_lt50/         1-49bp indel benchmark
            → dir/truvari_snp/          SNP benchmark
```

`svp run` cascade is extended: `lcp → graph → modea → modeb → posmap → anchor → vcf`. Each stage is resumable (checkpoint exists ⇒ skip).

---

## Stage 1: `svp posmap` (position oracle)

**Goal**: build a sig → (chrom_id, pos) lookup so any LCP signature seen in the reference can be located on GRCh38.

**File**: `src/posmap.cpp` (~180 LOC). Add `run_posmap(dir, ref_path)` to `src/graph.h`. Dispatch in `src/main.cpp`. Auto-picked up by `Makefile`'s `$(wildcard src/*.cpp)`.

**Algorithm**:

1. Accept `--level 5 --flank 10` CLI overrides (defaults). The level + flank MUST match the LCP pass exactly.
2. Read reference via `read_fasta_flat(ref_path, flank, seq, rec_start, rec_end)` from `src/io.cpp:7-28` — gives flat concatenated seq + per-record start/end offsets. Capture chromosome names by extending `read_fasta_flat` to also append a `names` vector (small additive change to `src/io.cpp`).
3. Extract cores via `extract_cores(seq, rec_start, rec_end, level, cores)` from `src/cores.cpp:19-31`. Expected ≈91M cores at WG L5 (matches the deleted-v1-posmap log).
4. Compute signatures via `compute_sigs(seq, cores, flank, 0, rec_start, rec_end, sigs)` from `src/cores.cpp:43-50` — same Murmur64 + label-seed scheme as the read pipeline.
5. For each core `i`: determine `chrom_id` by binary-searching `cores[i].start` in `rec_start`, then `pos = cores[i].start - rec_start[chrom_id]`. **Convention**: emit the core START (not the flank-window left edge). VCF emission later adds the flank-aware offset.
6. Sort records by sig (std::sort on `vector<{u64,u32,u32}>`) so anchor lookup is O(log N) via binary search.
7. Emit atomically: `sig_ref_map.bin` = `[u64 N][N × (u64 sig, u32 chrom_id, u32 pos)]`. Total ≈ 91M × 16 B = **1.5 GB** (smaller than the deleted v1 which was 2.9 GB).
8. Emit `chroms.tsv` = `chrom_id\tname\tlength` lines, one per record from the FAI.

**Function inventory** (each ≤30 LOC per `docs/CODE_STYLE.md`):
- `read_fasta_flat_named` (extension of `io.cpp` to also return names)
- `find_chrom_id(global_pos, rec_start)` — binary search
- `write_posmap(dir, records, names_lengths)` — atomic write of `.bin` + `.tsv`
- `run_posmap(dir, ref_path)` — driver

**Runtime budget**: cores extract ~30s (per `outputs/svp/posmap_20260602_030926.log`: 91M cores in 50s including IO+sort). Peak RAM ~24 GB (flat seq + cores + sigs simultaneously); release seq/cores before sort to drop to ~3 GB.

**Important**: do NOT reuse `outputs/svp/wgs_timed/posmap.bin` — its schema is from the deleted v1 caller and unknown. Build fresh.

**Skip a refresh**: if `dir/sig_ref_map.bin` exists, `svp run` skips this stage.

---

## Stage 2: `svp anchor` (per-bubble chrom + position)

**Goal**: for each Mode B bubble, look up its anchor super-nodes' constituent sigs in the reference oracle and assign (chrom, pos_start, pos_end).

**File**: `src/anchor.cpp` (~220 LOC). Add `run_anchor(dir)` to `src/graph.h`. Dispatch in `src/main.cpp`.

**Algorithm**:

1. Load `dir/modea.bin`'s super-node table + chain_ids (reuse `read_modea` from `src/modeb.cpp:52-68`). Load `dir/counts.bin`'s sig array (use `SigCount` from `src/svp.h:18`). Build the `node_id → sig` map exactly as `tools/bubble_recall.cpp:85-96` does.
2. Load `dir/sig_ref_map.bin` into a sorted vector. (Memory-resident; ~1.5 GB.)
3. Load `dir/modeb.bin`'s `BubbleTag` array (reuse the loader in `src/modeb.cpp`).
4. **For each bubble (s1, s2)**:
   - Gather all sigs from s1's and s2's chain_ids: `chain_ids[chain_off[s1]..chain_off[s1+1]]`, mapped to sigs via `counts_sigs[node_id]`.
   - For each sig, binary-search `sig_ref_map.bin`. Some sigs miss (sequencing-error sigs, tumor-only sigs not in ref); count hits per anchor.
   - **Majority-vote chrom**: take the modal `chrom_id` across all hits in s1 ∪ s2. If <50% of hits agree on one chrom → mark `chrom=UINT32_MAX` (unanchored).
   - **pos_start** = `min(positions of s1's hits on consensus chrom)`. **pos_end** = `max(positions of s2's hits on consensus chrom)`.
   - **RC-anchored case**: if `pos_end < pos_start`, swap and set `strand=1` (RC). The anchor super-node was extracted from the RC strand of the reference; downstream VCF emission RCs the residual bytes.
5. Emit `dir/bubble_pos.bin`:
   ```
   u64 B
   B × { u32 bubble_id, u32 chrom_id, u32 pos_start, u32 pos_end,
         u32 n_anchor_s1, u32 n_resolved_s1, u32 n_anchor_s2, u32 n_resolved_s2,
         u8 strand, u8 _pad[3] }
   ```
   28 B per record × 37,654 ≈ 1 MB.

**Filtering insight**: from a quick spot-check, c3 (shared backbone) anchors should resolve at >95% rate against a clean germline reference. Anchors with <50% resolution rate are repeat-region noise — flagged for VCF-stage drop.

**Function inventory** (≤30 LOC each):
- `load_posmap_sorted(path, &records)`
- `lookup_sig(sig, records)` — std::lower_bound
- `gather_anchor_positions(super_node_id, chain_ids, sigs, posmap, &hits)`
- `majority_chrom(hits, &chrom_id, &fraction)`
- `process_bubble(bt, modea, sigs, posmap, &out_pos)`
- `write_bubble_pos(dir, records)`
- `run_anchor(dir)` — driver

**Verification (Phase 0 of implementation)**: spot-check 50 random bubbles. For each, manually compute expected position by reading s1's first chain_id's sig, looking it up in posmap, and printing. Compare against `bubble_pos.bin`. Disagreement → bug.

---

## Stage 3: `svp vcf` (base emitter with positional filters)

**Goal**: walk bubbles, apply filters, emit `<dir>/calls.vcf` with per-bubble SNP/INS/DEL/DUP/INV records.

**File**: `src/vcf.cpp` (~400 LOC). Add `run_vcf(dir, out_path)` to `src/graph.h`. Dispatch in `src/main.cpp`.

### Inputs
- `dir/modeb.bin` (residual bytes + RLE coverage)
- `dir/bubble_pos.bin` (chrom + pos)
- `dir/modea.bin` (anchor super-node block bytes — needed to splice REF context)
- `dir/counts.bin` (for coverage tie-break in genotype calling)
- `simulation/ref.fasta` (via `read_fasta_flat_named` — holds in RAM at ~3 GB)
- `dir/chroms.tsv` (for `##contig=` header + numeric POS conversion)
- `--dup-flags <path>` (optional, from stage 4)
- `--inv-flags <path>` (optional, from stage 5)

### Per-tag emission rules (after filters)

| Mode B tag | Size | Output |
|---|---|---|
| S, \|N'\|=\|T'\|=1 | 1 bp | SNP: `CHROM POS . N' T' QUAL PASS INFO=SVTYPE=SNP FORMAT=GT:DV:DR  0/1:<t>:<n>` |
| S, \|N'\|=\|T'\|>1, no INV flag | <50 bp | MNP: `CHROM POS . N' T' QUAL PASS INFO=. FORMAT=GT:DV:DR  0/1:<t>:<n>` |
| S, ≥50bp, INV-flagged | ≥50 bp | INV: `CHROM POS . N <INV> QUAL PASS END=POS+\|N'\|-1;SVLEN=\|N'\|;SVTYPE=INV ...` |
| D | <50 bp | DEL: `CHROM POS-1 . ref[POS-1]+N' ref[POS-1] QUAL PASS SVTYPE=DEL;SVLEN=\|N'\|` |
| D | ≥50 bp | DEL symbolic: `CHROM POS . ref[POS] <DEL> QUAL PASS END=POS+\|N'\|-1;SVLEN=\|N'\|;SVTYPE=DEL` |
| I, DUP-flagged | any | DUP: `CHROM POS . ref[POS] <DUP> QUAL PASS END=POS+\|T'\|-1;SVLEN=\|T'\|;SVTYPE=DUP;TARGET=POS+\|T'\|-1` |
| I (not DUP) | <50 bp | INS: `CHROM POS . ref[POS] ref[POS]+T' QUAL PASS SVTYPE=INS;SVLEN=\|T'\|` |
| I (not DUP) | ≥50 bp | INS symbolic: `CHROM POS . ref[POS] <INS> QUAL PASS END=POS;SVLEN=\|T'\|;SVTYPE=INS;INSSEQ=T'` |
| C (complex) | any | **DROPPED** (measured 0% coverage filter pass; not worth emitting) |
| E | — | never emitted (modeb never produces E) |

### POS calculation

For D/I/S bubbles: the variant locus sits between s1 and s2 anchor sequences, offset by `pre` (shared prefix length) into the residual. Concrete formula:

```
POS_0based = pos_start + |s1.block| - suf - 1     (last base of s1's block tail, after trim)
POS_VCF    = POS_0based + 1
```

(`pos_start` is `bubble_pos.pos_start`, which we defined as s1's leftmost-anchored-sig position; `|s1.block|` is the spliced anchor block length from modea.bin.) Verification: for a SNP truth at chr1:125820 (G>A), the corresponding S-tag bubble's `POS_VCF` must equal 125820.

If strand=1 (RC-anchored), flip residuals: `N' ← RC(N')`, `T' ← RC(T')` before emit, and use `POS_0based = pos_end` (mirror).

### Filters (in order)

1. **Drop C-tag** entirely (measured 100% coverage-filter failure).
2. **Drop unanchored** (`chrom_id == UINT32_MAX`).
3. **Drop low-resolution anchors**: `n_resolved_s1 < 0.5 × n_anchor_s1` OR same for s2 (defines "consensus rejected" at the bubble-pos stage).
4. **Drop discordant strand**: `pos_end < pos_start` with `n_resolved_s2 > 0` (rare; usually means anchor overlap is broken).
5. **Anchor-distance sanity** for ≥50bp variants: `(pos_end - pos_start) ≤ max(2 × |residual|, 1000)` — drops bubbles whose anchors are scattered across genome (repeat-region phantoms). Note: this filter MUST NOT apply to true large DELs where `pos_end - pos_start ≈ SVLEN + flank`, so the cap = `2 × |residual|` is safe (verified against truth `1:42198870 sv17930 SVLEN=68782`).
6. **Genotype filter**: require `median(cov_t_runs) ≥ 3` for I-tag, `median(cov_n_runs) ≥ 3` for D-tag, both ≥ 3 for S-tag. (Already known to be ineffective in isolation, but cheap, runs after positional filters.)

### Genotype + QUAL

- `GT = 0/1` for somatic events (default); `1|1` for tumor-only with `n=0`. Simulation is haploid (`1|1`) but truvari accepts both.
- `DV = median(cov_t_runs)`, `DR = median(cov_n_runs)`.
- `QUAL = 60` (fixed PASS until we calibrate).
- `FILTER = PASS` (single-tier).

### Header

- `##fileformat=VCFv4.2`
- `##source=svp`
- `##contig=<ID=...,length=...>` for each contig in `chroms.tsv` (455 contigs, matching truth)
- `##INFO=<ID=SVTYPE,Number=1,Type=String,Description=...>` for each of SVTYPE / SVLEN / END / INSSEQ / TARGET
- `##ALT=<ID=DEL,Description="Deletion">` for each symbolic ALT
- `##FORMAT=<ID=GT,...>`, `##FORMAT=<ID=DV,...>`, `##FORMAT=<ID=DR,...>`
- `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT TUMOR`

### Function inventory (≤30 LOC each)
- `load_chroms(dir, &names, &lengths)`
- `load_dup_flags(path, &flags)` — sparse map bubble_id → DupCall
- `load_inv_flags(path, &flags)`
- `splice_anchor_block(modea, super_id, &block)` — reuse Mode A logic
- `compute_pos(bp, modea, anchor_blocks)` — POS_VCF
- `emit_snp(rec, ...)`, `emit_mnp(rec, ...)`, `emit_del(rec, ...)`, `emit_ins(rec, ...)`, `emit_dup(rec, ...)`, `emit_inv(rec, ...)` — one each
- `apply_filters(bt, bp, runs) -> bool` — returns whether to emit
- `write_header(out, chroms)`
- `run_vcf(dir, out_path)` — driver, sorts records by (chrom_id, POS) before write

---

## Stage 4: DUP detection (substring + coverage)

**Goal**: identify I-tag bubbles that are actually small/medium tandem duplications, plus tandem DUPs that ride the 105 small cycles. Emit `dir/dup_flags.bin` consumed by `svp vcf`.

**File**: `tools/dup_detect.cpp` (~200 LOC). Builds standalone for now; integrate into `svp vcf` after validation.

### Detector A: ref-context substring (catches small + medium DUPs)

For each I-tag bubble:
1. Fetch ref window: `ref[chrom][max(0, pos_start - W) : min(L, pos_end + W)]` with `W = max(200, 2×|T'|)`.
2. Check if `T'` is a substring of the window (forward).
3. Check if `RC(T')` is a substring of the window.
4. If yes → DUP candidate.

Justification: tandem DUP `ABCDE → ABCDEABCDE` produces I-tag `T'=ABCDE`; ABCDE is already in the ref context at the anchor site. True INS introduces a sequence NOT in the local ref. Use `memmem` from `<string.h>`.

### Detector B: coverage uplift (catches large tandem DUPs that don't substring-match cleanly)

1. Compute global `cov_t_bg` = median tumor coverage across ALL c3 super-nodes (one pass over modea.bin runs, target ~30× per CLAUDE.md).
2. For each I-tag bubble in the `T_only_no_N` bucket (bucket=1): compute `median(cov_t_runs of T arm)`. If `≥ 1.5 × cov_t_bg` → DUP candidate.

### Detector C: small-cycle bonus

The 105 cycle-SCCs (size 2-2280) in `modea.bin` are realistic tandem-repeat candidates. For each cycle SCC entirely in c2 (tumor-only): emit a DUP candidate spanning the cycle members' total block bytes, anchored at the cycle's entry/exit c3 nodes.

### Output schema (`dir/dup_flags.bin`)

```
u64 K
K × { u32 bubble_id, u32 svlen, u8 mechanism /* 'A'/'B'/'C' */, u8 _pad[3] }
```

### Function inventory
- `compute_cov_bg(modea) -> uint32_t` — global c3 tumor median
- `check_substring(t_prime, ref_window) -> bool` — `memmem` forward + RC
- `check_cov_uplift(bt, runs, bg) -> bool`
- `enumerate_cycles_from_modea(modea, color_filter, &cycle_groups)`
- `emit_dup_flags(path, calls)`

### Verification
Run on `outputs/svp/wgs_mm10/` and check: of the 32 truth DUPs, how many are flagged? Target: ≥29/32 = 90.6%.

---

## Stage 5: INV detection (best-effort cross-bubble RC pairing)

**Goal**: pair I-tag and D-tag bubbles whose residuals are reverse complements at the same locus → INV call. Emit `dir/inv_flags.bin`.

**File**: `tools/inv_detect.cpp` (~180 LOC).

### Algorithm: positional cross-bubble RC pairing

1. Build a position-bucketed index of bubbles: `bucket[chrom][pos // 5000] = list of (bubble_id, tag, residual_off, residual_len)`.
2. For each I-tag bubble `i` with `|T'| ≥ 50bp`:
   - Iterate bubbles in the same and adjacent buckets (±5kb window) of the same chromosome.
   - For each D-tag bubble `d` in range: check if `T'_i == RC(N'_d)` with size tolerance ±10%. If yes → INV pair.
3. INV call: `POS = min(i.pos_start, d.pos_start)`, `END = max(i.pos_end, d.pos_end)`, `SVLEN = max(|T'|, |N'|)`.

### Alternative attempt: lone-S-tag fallback

For each S-tag bubble with `|N'|=|T'|>1`: check `N' == RC(T')`. Already measured 0/250 hits; expected to add 0 INVs but cheap to include for completeness.

### Output schema (`dir/inv_flags.bin`)

```
u64 K
K × { u32 bubble_id_a, u32 bubble_id_b /* may be UINT32_MAX for lone S */,
      u32 pos_start, u32 pos_end, u32 svlen, u8 mechanism, u8 _pad[3] }
```

### Realistic expectations
Per `docs/INSPECTOR_FINDINGS.md`, truth INVs decompose into 6 sub-bubbles each. Cross-bubble pairing in a ±5kb window should catch at least the boundary-spanning pair. Target: ≥30% INV recall (5/15). The user explicitly accepted lower INV recall.

### Function inventory
- `build_position_index(bps, &buckets)`
- `rc_bytes(in, len, &out)` — LUT-based RC (reuse `tools/lcp_rc_invariance.cpp:23-32`'s table)
- `find_rc_pair(i_bt, candidates, runs_arr) -> int` — returns d's index or -1
- `emit_inv_flags(path, calls)`

---

## Stage 6: `scripts/svp_benchmark.sh` (truvari driver)

**Goal**: prep VCFs and run truvari at three thresholds.

**File**: `scripts/svp_benchmark.sh` (~80 LOC). Wraps the existing `scripts/truvari.sh` pattern.

### Steps

1. Source `.venv/bin/activate` (the project venv at `/home/yuio/svp/.venv` — pyfaidx, truvari v5.4.0, bcftools 1.21).
2. `bcftools sort` + `bgzip -c` + `tabix -p vcf` the truth VCF and the caller's `calls.vcf` (idempotent — skip if `.gz.tbi` exists).
3. Run three truvari benches:

```bash
# A: ≥50bp SV benchmark (primary target)
truvari bench -b truth.vcf.gz -c calls.vcf.gz -f simulation/ref.fasta \
   --sizemin 50 --sizefilt 50 --pctsize 0.7 --pctseq 0 \
   --refdist 500 --sizemax 1000000000 \
   -o $DIR/truvari_geq50

# B: 1-49bp indel benchmark (small variants)
truvari bench -b truth.vcf.gz -c calls.vcf.gz -f simulation/ref.fasta \
   --sizemin 1 --sizefilt 1 --pctsize 0.7 --pctseq 0 \
   --refdist 50 --sizemax 49 \
   -o $DIR/truvari_lt50

# C: SNP-only benchmark (truvari covers SNPs at sizemin=0)
truvari bench -b truth.vcf.gz -c calls.vcf.gz -f simulation/ref.fasta \
   --sizemin 0 --sizefilt 0 --sizemax 1 \
   --refdist 0 \
   -o $DIR/truvari_snp
```

4. Print `summary.json` from each run, plus per-SVTYPE FN/FP breakdowns (reuse the FN/FP entry-printer block from `scripts/truvari.sh:54-104`).

---

## Files to create

| Path | Purpose | LOC |
|---|---|---|
| `src/posmap.cpp` | reference → sig→pos oracle | ~180 |
| `src/anchor.cpp` | bubble → genome position | ~220 |
| `src/vcf.cpp` | VCF emission with filters | ~400 |
| `tools/dup_detect.cpp` | DUP flag emitter | ~200 |
| `tools/inv_detect.cpp` | INV flag emitter | ~180 |
| `scripts/svp_benchmark.sh` | truvari driver | ~80 |
| `docs/VCF_EMISSION.md` | architecture + tag rules | ~150 |

Total new C++: ~1,360 LOC. Plus `src/io.cpp` extension (~20 LOC for chrom-name capture), `src/main.cpp` extension (~30 LOC for new dispatch entries + run cascade), `src/graph.h` extension (~3 LOC for new function declarations), `Makefile` extension (~6 LOC for `dup_detect` / `inv_detect` build rules).

## Files to modify

- `src/io.cpp` — add `read_fasta_flat_named` (returns chrom names too)
- `src/main.cpp` — add `posmap`, `anchor`, `vcf` dispatch + cascade skip checks
- `src/graph.h` — add `run_posmap`, `run_anchor`, `run_vcf` declarations
- `Makefile` — add `tools/dup_detect` and `tools/inv_detect` build rules; reference posmap/anchor/vcf checkpoints in `make run` if cascading
- `CLAUDE.md` — extend Pipeline section with stages 6, 7, 8 (posmap/anchor/vcf) + new checkpoint list

## Existing utilities reused

- `extract_cores` / `compute_sigs` — `src/cores.cpp:19-50`
- `read_fasta_flat` (extended) — `src/io.cpp:7-28`
- `SigCount` POD — `src/svp.h:18`
- `xopen` / `wr` / `rd` / `atomic_close` IO helpers — `src/graph.cpp` (visible via `svp.h`)
- `read_modea` (super-node + chain_ids + super-edges loader) — `src/modeb.cpp:52-68` (pull into `modea.h` or duplicate verbatim)
- `BubbleTag` parser logic — `tools/modeb_validate.py:152-167` for Python verification; mirror in C++ in `vcf.cpp`
- `suffix_prefix` (KMP overlap) — `src/modea.cpp:16-31` (reused by `anchor.cpp` for arm splicing if needed)
- RC LUT — `tools/lcp_rc_invariance.cpp:23-32`
- truvari prep flow — `scripts/truvari.sh` (sort/bgzip/tabix idiom)

## Out of scope

- **Read-substring caching** during the LCP pass (deferred — would let us emit exact INSSEQ for novel insertions instead of T' residuals; current T' is good enough for substring-match truvari).
- **Multi-sample / diploid** genotyping — simulation is haploid 1|1.
- **Multi-bubble locus merging** (Sniffles2 style) — start without; if precision <90%, add a simple position-distance merge as a post-pass.
- **Phased calls / breakends (BND)** — not needed for SV-class truvari (DEL/INS/DUP/INV) at this stage.
- **Coverage-based genotype likelihoods** — fixed `0/1` for now.
- **Complex SV chains** (e.g., TRA/INV+DUP composites) — out of insilicoSV truth and out of scope.
- **Re-using `outputs/svp/wgs_timed/posmap.bin`** — schema unknown; build fresh.

## Verification

End-to-end sequence (each step independently runnable):

1. **Stage-by-stage build & smoke test**:
   ```bash
   make && ./svp posmap simulation/ref.fasta outputs/svp/wgs_mm10
   ./svp anchor outputs/svp/wgs_mm10
   ./svp vcf outputs/svp/wgs_mm10
   ./tools/dup_detect outputs/svp/wgs_mm10
   ./tools/inv_detect outputs/svp/wgs_mm10
   ./svp vcf outputs/svp/wgs_mm10 --dup-flags outputs/svp/wgs_mm10/dup_flags.bin \
                                  --inv-flags outputs/svp/wgs_mm10/inv_flags.bin
   bash scripts/svp_benchmark.sh outputs/svp/wgs_mm10
   ```

2. **Quantitative gates** after step 1:
   - `posmap` log line: `[posmap] 91,000,000+ sig records / 455 chroms loaded`. If <80M, sig extraction is broken.
   - `anchor` log line: ≥90% of bubbles have `chrom_id != UINT32_MAX` (resolved). If lower, the anchor super-nodes' sigs don't appear in the reference at expected rate; check level/flank match.
   - `vcf` emits ≤ 5,000 records ≥50bp (we have 3,078 bubbles ≥50bp; post-filter expect ~1,500-2,000 emissions). If 7M like v24, position filters are not being applied.
   - `truvari_geq50/summary.json`: recall ≥90% (target = 101/112 TPs) AND precision ≥90% (FP ≤ 11). MUST exceed v24's 88.4% recall.
   - `truvari_lt50/summary.json`: recall ≥90% on 940 small indels.
   - `truvari_snp/summary.json`: recall ≥90% on 16,950 SNPs.

3. **Hand-check 6 bubbles** (one per SVTYPE) end-to-end:
   - Pick lowest-position truth SNP, INS, DEL, DUP, INV from truth VCF.
   - Find the Mode B bubble whose `bubble_pos` matches its locus (or closest).
   - Verify the emitted VCF record matches expected CHROM/POS/SVTYPE/SVLEN. Disagreement → `svp vcf` bug.

4. **FN/FP audit** after stage 6:
   - Examine `fn.vcf` from `truvari_geq50`: for each missed truth SV, is the Mode B bubble present at that locus? If yes → filter dropped it incorrectly. If no → upstream LCP/mm=10 lost the signal.
   - Examine `fp.vcf` from `truvari_geq50`: cluster FPs by chrom region; identify whether they're repeat-region or non-repeat. Repeat-region FPs → consider an LCP-repeat-mask filter as a v2.

5. **Stretch goal regression check**: after each stage 4/5 detector tweak, re-run only `svp vcf` (cheap) + truvari to measure delta on DUP/INV recall.

## Phase ordering (work over multiple sessions)

1. **Phase 0** (≤30 min): build `tools/bubble_recall` if absent + measure truth DUP/INV survival in `counts.bin` (using existing `tools/truth_sig_dump.cpp` pattern) — confirms signal exists before we build the chain. **If signal missing at counts.bin level, mm=10 is too strict and we must back off.**
2. **Phase 1** (one session): ship `src/posmap.cpp` + `src/io.cpp` extension + run on WG → produce `sig_ref_map.bin`. Quantitative gate: ≥91M records, ~50s runtime.
3. **Phase 2** (one session): ship `src/anchor.cpp` → `bubble_pos.bin`. Quantitative gate: ≥90% bubbles resolved.
4. **Phase 3** (one session): ship `src/vcf.cpp` minimum-viable (SNP/INS/DEL only, no DUP/INV) → `calls.vcf`. Run `scripts/svp_benchmark.sh`. Quantitative gate: SNP recall ≥90%, DEL recall ≥80% (no DUP folded in yet).
5. **Phase 4** (one session): ship `tools/dup_detect` → re-run vcf with `--dup-flags` → re-run truvari. Quantitative gate: DUP recall ≥80%.
6. **Phase 5** (one session): ship `tools/inv_detect` → re-run with `--inv-flags` → re-run truvari. Quantitative gate: INV recall ≥30%.
7. **Phase 6**: tune precision (FP-filter refinement, maybe locus merging). Quantitative gate: precision ≥90% across types.

Each phase is independently verifiable, so resuming mid-stream after a context compaction is safe — the previous phase's `.bin` checkpoint is the contract.
