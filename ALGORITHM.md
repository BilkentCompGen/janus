# Algorithm & Complexity

Notation used throughout:
- **N** = genome length (e.g. ~58M for human chr19)
- **d** = LCP depth parameter (`--depth`)
- **C** = number of LCP cores after level-1 parsing (~N/1.5)
- **C_d** = cores after d levels of deepening (~N/1.5^d)
- **V** = graph nodes (unique cores across all genomes), **E** = graph edges
- **B** = number of bubbles, **S** = number of bubble starts
- **L** = bubble path sequence length (varies per bubble)

---

## Phase 1: LCP Parsing

Each genome is parsed into **LCP cores** — locally consistent segments identified by local extrema patterns in the alphabet ordering.

**`parse1()`** scans the sequence left-to-right in a single pass. At each position it checks for three core types:
- **RINT** (run-interrupted): two equal characters flanking a run of identical characters, ended by a different character.
- **LMIN** (local minimum): pattern `a > b < c` — a valley in alphabet order.
- **LMAX** (local maximum): pattern `a < b > c` with additional boundary constraints.
- **SSEQ** (spanning sequence): gap-filling cores between identified extrema.

This produces ~N/1.5 cores. **O(N) time, O(N) space.**

**`lps_deepen()`** iteratively compresses cores to higher levels:
1. **DCT compression** (`lcp_dct`): Deterministic Coin Tossing — compares adjacent cores' bit representations, compressing the right core by encoding only the position/value of the first divergent bit. Single pass over the core array. **O(C) per level.**
2. **Re-parse** (`parse3`): Applies the same RINT/LMIN/LMAX/SSEQ pattern detection, but now on the core array instead of characters. Produces ~C/1.5 new cores. **O(C) per level.**

Each level reduces the core count by factor ~1.5. After d levels: **C_d ~ N / 1.5^d** cores.

**Total parsing cost:** O(N + N/1.5 + N/1.5² + ... + N/1.5^d) = **O(N)** (geometric series converges). Done 3 times (ref, normal, tumor).

---

## Phase 2: Graph Construction

**`LCPGraph::addString()`** — for each core in the parsed sequence:
1. Hash the core's underlying string (64-bit hash via `hash_string`).
2. Insert/merge into `nodes` HashMap, setting color bit (0x1=normal, 0x2=tumor).
3. Add forward/reverse edges between consecutive cores.

Shared cores (identical sequences in normal and tumor) map to the same node with both color bits set. Divergent regions produce single-color nodes.

**O(C_d) per genome, O(N/1.5^d) amortized.** Edge sets use `unordered_set` for O(1) amortized insert/lookup.

---

## Phase 3: Reference Indexing

**`RefIndex::addString()`** — indexes reference cores into a position-lookup structure:
- `cores[]`: ordered array of (hash, start, end) for each reference core.
- `index[]`: hash → list of positions (multimap). A core with exactly one position entry is **unique** and can serve as an anchor.

**`anchor(h)`**: O(1) hash lookup. Returns position only if the core hash is unique in the reference. **`anchorSpan(h)`**: same, returns (start, end).

**O(C_d) time and space.**

---

## Phase 4: Bubble Detection

**`detectBubbles()`** — identifies divergence boundaries in the graph by scanning all nodes:

A node is a **bubble start** if:
- It is shared (colors & 0x3 == 0x3), AND
- At least one successor is NOT shared (has only one color).

A node is a **bubble end** if:
- It is shared, AND
- At least one predecessor is NOT shared.

This captures the transition points: shared → divergent (start) and divergent → shared (end).

**O(V + E)** — each node examined once, each edge checked once.

---

## Phase 5: Path Finding & Bubble Extraction

**`extractBubbles()`** — for each bubble start, finds both the normal-only and tumor-only paths to the same bubble end.

**`findColoredPath(start, color_bit)`** — BFS from `start`, following only edges to nodes with the target color and NOT the other color:
1. Enqueue successors of `start` that are single-color.
2. At each step, **merge** the visited node's sequence into an accumulator using **KMP overlap detection** (`Node::merge()`).
3. Stop when reaching a bubble end node. Extract the divergent substring (between start node suffix and end node prefix overlaps).
4. **Underflow guard**: If `before < start_len + end_overlap` (the start and end node overlaps fully cover the intermediate path), the extracted sequence is set to empty rather than producing garbage from unsigned underflow. This occurs for small tandem duplications where the single normal-path core is entirely overlapped by the start and end nodes.

**`Node::merge()`** — computes the longest suffix-prefix overlap between the accumulator and the next node's sequence using the KMP failure function. This correctly handles the overlapping nature of LCP cores (adjacent cores share boundary characters).

- KMP failure function: O(|pattern|) to build.
- Matching against accumulator suffix: O(min(|acc|, |pat|)).

**Per bubble:** O(path_length) for the BFS + O(L) for sequence accumulation, where L is the total path sequence length.

A bubble is **valid** only if normal and tumor paths converge to the **same** end node. Invalid (open-ended) paths are discarded.

**Total: O(sum of all path lengths across all bubbles).**

---

## Phase 6: Reference Anchoring

**`anchorBubble(bubble, graph)`** — maps a bubble's position back to reference genome coordinates:

1. **Direct anchor**: Try `anchor(start_node)`. If the start node's hash is unique in the reference, return its end position. **O(1).**

2. **Backward walk**: BFS through predecessor edges from the start node. For each predecessor, check if it has a unique reference anchor. If found, scan forward in the reference core array to locate the start node's position. **O(predecessors_visited + C_d) worst case.**

3. **Forward fallback**: Try `anchor(end_node)` and compute position by subtracting the normal sequence length. **O(1).**

**Typical case: O(1)** (most bubble start nodes anchor directly). **Worst case per bubble: O(C_d)** for the linear scan.

---

## Phase 7: SV Calling

### Quick classification — `detectSV(sequence, query)`

Scans inward from both ends of the longer and shorter sequence simultaneously to find the divergent region:
- If boundaries cross: **no_call** (sequences are identical within the bubble).
- If single position: **snv**.
- If shorter sequence boundaries collapse: **indel** (clean insertion or deletion).
- Otherwise: **multi** (complex variant requiring alignment).

**O(L) where L = max(|normal_seq|, |tumor_seq|).**

### Simple calls (snv, indel)

Directly emit a single SVCall from the detected region. **O(1).**

### Complex calls — `callMulti(normal_seq, tumor_seq)`

Uses **parasail Needleman-Wunsch** (SIMD-accelerated global alignment) with:
- Scoring: match=5, mismatch=-4, gap_open=50, gap_extend=1.
- High gap-open penalty biases toward fewer, larger indels rather than fragmented small ones.

**NW alignment: O(L_n * L_t) time and space** (parasail uses SIMD striped vectorization for practical speedup, but asymptotic bound is quadratic).

CIGAR parsing then groups operations into **variant blocks** separated by anchor matches (>=50bp):
- **Indel blocks** (contain D/I ops): consolidated into a single call per block. Prevents fragmentation of nearby DEL+INS into hundreds of tiny alternating calls.
- **SNP blocks** (only mismatches + short matches): each mismatch emitted as an individual 1bp SNP call, preserving SNP resolution.

**CIGAR processing: O(|CIGAR|) = O(L_n + L_t).**

---

## Phase 8: Duplication Detection

Tandem duplications are detected through two complementary mechanisms: **pattern-based flipped DUP detection** for medium/large DUPs, and **Smith-Waterman INS→DUP reclassification** for small DUPs.

### 8a: Flipped DUP Detection (`detectFlippedDUPs`)

A tandem duplication `[A][B] → [A][B][B]` creates a distinctive graph signature: the extra copy produces tumor-only nodes at the DUP boundaries, while the duplicated content itself remains shared (it exists in both genomes). This manifests as bubble start/end nodes with a specific color pattern.

**Observation from data**: Medium/large DUPs (>~300bp) produce bubble starts near the DUP **end** and bubble ends near the DUP **start** — positions are **flipped** relative to normal bubbles. Both the start's successors and the end's predecessors show the `T=1 S≥1 N=0` pattern (one tumor-only neighbor, one or more shared neighbors, zero normal-only neighbors). This pattern does not appear at DEL, INS, or SNP sites.

**Algorithm** (replaces the previous cycle detection BFS):

1. **Filter**: For each bubble start and end, call `anchorSpan()`. Discard any without a unique reference anchor (non-unique nodes create meaningless bubbles).

2. **Classify**: Count successor colors for each anchored bubble start (`countSuccessorColors`). A DUP candidate start has `tumor_only ≥ 1`, `shared ≥ 1`, `normal_only = 0`. Same for bubble ends with predecessor colors (`countPredecessorColors`).

3. **Pair**: Sort candidate starts by ref position descending, candidate ends ascending. For each start, find the nearest unmatched end where `end_ref < start_ref` (flipped order). Skip pairs already handled by normal bubble extraction.

4. **Emit**: For each pair: `pos = end_ref_start`, `svlen = start_ref_start - end_ref_start`. Emit as DUP call.

**Complexity**: O(S_dup · E_dup) for the greedy pairing where S_dup and E_dup are the number of DUP candidate starts and ends (typically small subsets of total bubble boundaries).

### 8b: INS→DUP Reclassification (Smith-Waterman local alignment)

Small DUPs (<~300bp) are detected as normal bubbles with `N=1 T=1` pattern. The bubble extraction produces both normal and tumor paths, and `callSV` classifies them as INS. However, the inserted sequence is a **rotated fragment** of the duplicated reference region (due to LCP core boundary alignment), not the exact duplicated sequence. This means simple character-by-character flank comparison fails.

**Algorithm**: For each INS call with `length < 2000bp` and non-empty `tumor_seq`:

1. Extract a reference window around the insertion point: `ref[pos - pad .. pos + pad]` where `pad = max(300, 2 × tumor_seq.size())`.
2. Run **parasail Smith-Waterman** (local alignment) of `tumor_seq` against the reference window.
   - Scoring: match=2, mismatch=-4, gap_open=5, gap_extend=1.
3. Compute normalized score: `sw_score / (tumor_seq.size() × match_score)`.
4. If normalized score ≥ 0.70: reclassify as DUP.

The Smith-Waterman local alignment finds the best matching substring regardless of rotation or positional shift, which is why it succeeds where character-by-character comparison fails.

**Complexity**: O(|tumor_seq| × |window|) per INS call. Since both are bounded by constants (tumor_seq < 2000bp, window < ~4300bp), this is O(1) per call.

---

## Phase 9: Output

- **VCF**: O(number_of_calls * avg_allele_length).
- **GFA**: O(V + E) for segments/links. Walk generation finds entry nodes per color (nodes without same-color predecessors), then follows single-color edges. O(V + E) total across all colors.

---

## Summary Table

| Phase | Operation | Time Complexity | Space |
|-------|-----------|----------------|-------|
| LCP parsing | `init_lps` + `lps_deepen` (×3 genomes) | O(N) | O(N) |
| Graph construction | `addString` (×2 genomes) | O(C_d) | O(V + E) |
| Ref indexing | `RefIndex::addString` | O(C_d) | O(C_d) |
| Bubble detection | `detectBubbles` | O(V + E) | O(B) |
| Path finding | `extractBubbles` | O(sum of path lengths) | O(L_max) |
| Anchoring | `anchorBubble` (×B) | O(B) typical, O(B·C_d) worst | O(V) for BFS |
| SV calling (simple) | `detectSV` + emit | O(L) per bubble | O(1) |
| SV calling (complex) | `callMulti` (NW alignment) | **O(L_n · L_t)** per bubble | O(L_n · L_t) |
| Flipped DUP detection | `detectFlippedDUPs` | O(S_dup · E_dup) | O(S_dup + E_dup) |
| INS→DUP check | SW local alignment | O(1) per INS call | O(window) |
| Output | VCF + GFA | O(V + E + calls) | O(V + E) |

**Dominant costs:**
1. **LCP parsing**: O(N) — linear in genome size, run once per genome.
2. **NW alignment**: O(L_n · L_t) — quadratic in bubble sequence lengths. This is the bottleneck for large complex variants, mitigated by parasail's SIMD vectorization.

**Overall: O(N + B · L²)** where L is the average bubble sequence length. DUP detection adds negligible cost.

**Space: O(N + V + E)** dominated by genome storage and the graph.

---

## Benchmark Results (depth=5, chr19, 2026-03-26)

### All variants (deep_analysis.py — strict type matching, 200bp tolerance)

| Type | Sim | Detected | TP | Wrong type | Ghost FP | FN | Precision | Recall | F1 |
|------|-----|----------|-----|-----------|----------|-----|-----------|--------|--------|
| DEL | 195 | 187 | 187 | 0 | 0 | 8 | 1.0000 | 0.9590 | 0.9791 |
| INS | 195 | 192 | 181 | 10 | 1 | 4 | 0.9427 | 0.9784 | 0.9602 |
| SNP | 200 | 199 | 199 | 0 | 0 | 1 | 1.0000 | 0.9950 | 0.9975 |
| DUP | 195 | 167 | 149 | 18 | 0 | 28 | 0.8922 | 0.8418 | 0.8663 |
| Mixed | 3281 | 3103 | 3030 | 71 | 2 | 180 | 0.9765 | 0.9439 | 0.9599 |

### SVs only (>=50bp, deep_analysis.py)

| Type | GT SVs | Det SVs | TP | FP | FN | Precision | Recall | F1 |
|------|--------|---------|-----|-----|-----|-----------|--------|--------|
| DEL | 143 | 135 | 135 | 0 | 8 | 1.0000 | 0.9441 | 0.9712 |
| INS | 147 | 147 | 146 | 1 | 1 | 0.9932 | 0.9932 | 0.9932 |
| DUP | 142 | 121 | 109 | 12 | 23 | 0.9008 | 0.8258 | 0.8617 |
| Mixed | 394 | 391 | 369 | 22 | 23 | 0.9437 | 0.9413 | 0.9425 |

### Truvari Results (SV-level, >=50bp)

| Type | TP | FP | FN | Precision | Recall | F1 |
|------|-----|-----|-----|-----------|--------|--------|
| DEL | 135 | 1 | 8 | 0.9926 | 0.9441 | 0.9677 |
| INS | 145 | 2 | 2 | 0.9864 | 0.9864 | 0.9864 |
| DUP | 28 | 94 | 114 | 0.2295 | 0.1972 | 0.2121 |
| Mixed | 288 | 105 | 106 | 0.7328 | 0.7310 | 0.7319 |

### Key observations

- **DEL, INS, SNP**: Excellent accuracy. DEL and INS SV-level F1 > 0.97. SNP F1 = 0.998. Zero ghost false positives for DEL and SNP.
- **DUP deep_analysis F1 = 0.87**: The flipped DUP detection catches 92 medium/large DUPs. The SW-based INS→DUP reclassification catches 57 small DUPs. 18 DUPs remain misclassified as INS (SW score below 0.70 threshold). 28 DUPs are completely missed (25 of those are <50bp and produce no bubble at all).
- **DUP truvari F1 = 0.21**: Truvari uses different matching logic that does not handle the DUP↔INS type ambiguity well. The deep_analysis results are more representative of actual detection accuracy.
- **Mixed**: 3030/3281 correct type calls. The 71 wrong-type calls include 50 real INS falsely reclassified as DUP (SW local alignment matches reference by chance for small insertions), 18 real DUPs still called as INS, and 3 spurious long-range DUP pairings (SVLEN > 3M) that matched near SNP sites.
- **INS→DUP false reclassifications**: 10 real insertions in the INS-only run are reclassified as DUP (all <25bp). At the SV level (>=50bp), this noise disappears — INS SV F1 = 0.993.

### DUP detection breakdown by mechanism

| Mechanism | DUP-only run | Mixed run |
|-----------|-------------|-----------|
| Flipped DUP detection (medium/large) | 92 | 123 |
| INS→DUP reclassification (small) | 57 | 74 |
| Still called as INS (should be DUP) | 18 | 18 |
| Not detected at all | 28 | 135 (132 are <50bp) |

### Known issues

1. **Spurious long-range DUP pairings**: The greedy nearest-match in `detectFlippedDUPs` can pair distant `T=1 S=1` bubble endpoints, producing DUP calls with SVLEN in the millions. 3 such calls in the mixed run. A max-distance cap would fix this.
2. **Small INS→DUP false positives**: The 0.70 SW threshold catches small insertions that happen to match nearby reference by chance. Affects ~10 calls in the INS-only run, all <25bp.
3. **Small DUPs (<50bp) missed**: 25 of 28 missed DUPs in the DUP-only run are <50bp. These produce no bubble at all — the LCP core tiling at depth 5 does not create divergent paths for variants smaller than a single core.
