# RC edges for partial-overlap strand pairing

**Status**: future work. Not implemented.

**Motivation**: Stage 0 validation (`tools/rc_lcp_equiv`, run on 20 reads /
~287 kb / 3982 cores from `simulation/reads/normal/hifi.fastq.gz`) found that
two simpler approaches to pooling forward and RC observations into a single
sig count both fail:

| approach | pooling rate | notes |
|---|---:|---|
| exact RC-mirror cores: forward (st,en,lb) ↔ RC (N−en, N−st, lb) | **0.0000 %** (0 / 3982) | LCP chiral asymmetry is total — no forward core ever sits at the RC-mirror position with the same label |
| canonical-window hash: `murmur64(min(W, RC(W)), 0)` | **0.1005 %** (4 / 3980) | matches the prior WG measurement (`feedback_rc_pairing_byte_test_fails`: 0.07 %) |

The forward LCP cores in read R and the LCP cores in RC(R) do not occupy
byte-RC-mirror positions, so neither their bytes nor their hashes coincide.
Strand pooling cannot be achieved at the per-core hashing level alone.

What does empirically hold (confirmed by Stage 0): the flank window of a
forward core and the flank window of an RC-extracted core covering the same
biological locus **do overlap** as byte strings when we map back into the
original read coordinates — they share a substring even though neither is
the other's exact byte-RC. That overlap is the signal RC edges would exploit.

## Design

After the existing `lcp` stage builds `counts.bin` (with per-sig flank
payload), run a post-process that for each sig X computes RC(payload[X])
and finds sigs Y whose payload partially overlaps RC(payload[X]):

```
for each sig X with payload W_X:
    R_X := RC(W_X)
    find sigs Y where payload[Y] occurs as a prefix/suffix/substring of R_X
         (or vice versa: R_X occurs inside payload[Y])
    emit an "RC edge" record {X → Y1, Y2, ..., Yk}
```

Output: `dir/rc_edges.bin` — packed `(src_sig: u64, dst_sig: u64)` pairs,
undirected (store each pair once with src < dst).

Downstream `run_graph` extends to treat each connected component of RC-edges
as a single canonical node, summing the component members' `n` and `t`
counts into the canonical node. This is the strand-pooling step that
canonical-hash hashing cannot achieve.

## Open questions

1. **Overlap threshold**. Minimum overlap length? An exact byte match of ≥k
   contiguous bases, with k chosen empirically? Hamming-1 to absorb
   sequencing errors? Affix-only (prefix or suffix of length ≥k) is probably
   the right starting point since LCP chirality tends to shift cores by a
   few bases, not embed them inside each other.
2. **Edge directionality**. RC is symmetric (RC(RC(W)) == W), so RC edges
   are undirected. Store once, by canonical key (min(src, dst)).
3. **Performance**. Whole genome has ~89 M nodes × 173 B avg payload (per
   `reference_seqs_bin_schema`). Naive all-pairs RC overlap is O(N²) = 8e15
   comparisons — infeasible. The standard fix is a generalized suffix array
   on the concatenated payload, query each `R_X` in O(|R_X| + occurrences).
   `libsais` or `divsufsort` plus an LCP array gets us to O(N log N) total.
   Expected wall: ~minutes on the WG counts.bin.
4. **Multi-component noise**. If many sigs end up in one giant RC-edge
   connected component (a "hairball"), summing collapses biology along with
   noise. Cap component size (e.g., reject components > some threshold,
   leave them un-pooled). Tune on the GT recall/precision table.
5. **Bytes that survive after pooling**. Each canonical node represents a
   strand-pair pool; downstream consumers (modea, inspect) need a canonical
   byte payload. Pick the lex-min over component members so the choice is
   deterministic.

## How this connects to the rest of the caller

- Stage `lcp` is unchanged. Same `counts.bin` / `edges.bin` layout.
- Stage `graph` gains a pre-step: load `rc_edges.bin`, union-find the sigs
  into RC components, rewrite the node table to use canonical-component
  representatives.
- Stage `modea` / `bubbles` are unchanged (they consume the rewritten
  graph).
- New stage between `lcp` and `graph`:
  `svp rc-edges <dir> [--min-overlap K] [--max-component M]`. Idempotent and
  resumable (`dir/rc_edges.bin` is the checkpoint).

## Expected outcome

If RC pooling works as hoped: per-sig coverage roughly doubles for sigs that
are part of a real strand pair (the ones currently giving GT mean count
≈ 12 → expected ≈ 24-25). Singleton phantom-error sigs are unlikely to
overlap anything via RC — they stay singleton and get pruned at higher
`min_mult`. The CLAUDE.md baseline `min_mult=3` could move to 8-9 with
matched or better recall.

If RC pooling collapses too much (component-size cap repeatedly triggers,
or recall drops despite higher counts): tune the overlap threshold up,
narrow to suffix-only matches, or scope RC edges to within-bubble pairs
only.
