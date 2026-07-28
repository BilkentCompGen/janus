# Truvari — what it actually checks, and how well svp currently does against it

_Written 2026-07-28, immediately before the direction change (blacklist + refvote removal,
DEL/INS/DUP/INV/**CPX** calling, core-level coordinates). This is the record of the
**reference-anchored** era: the method that scores best today and the numbers it scores._

Everything below is read out of the installed source (`.venv`, **Truvari v5.4.0**) or out of the
`summary.json` files the runs wrote. Nothing here is quoted from memory.

---

## 1. How Truvari decides a call is a TP

### 1.1 Chunking, then all-vs-all, then a picker
`bench.py` streams both VCFs in position order and cuts them into **chunks** — a chunk ends when the
next record starts more than `refdist` (500 bp) past the last one. Inside a chunk every base record is
compared against every comp record (`Bench.compare_chunk`), producing a match matrix. The default
picker is **`single`** (`pick_single_matches`): sort all pairs by score, greedily take the best pair,
retire both sides, repeat. So **one truth SV can be credited to at most one call** — a caller that
splits an event into three correct pieces gets 1 TP + 2 FP, not 3 TP.

### 1.2 The gates (`VariantRecord.var_match`, `variant_record.py:778`)
A pair is a match only if **all** of these pass:

| # | gate | default | what we pass |
|---|---|---|---|
| 1 | **same SVTYPE** (`same_type`) | on | on — `--typeignore` turns it off |
| 2 | **position**: `overlaps(bstart-refdist, bend+refdist, cstart, cend)` | `refdist=500` | default |
| 3 | **size similarity**: `min(a,b)/max(a,b) ≥ pctsize` | `pctsize=0.70` | default |
| 4 | **reciprocal overlap** ≥ `pctovl` | `pctovl=0.0` → **off** | default (off) |
| 5 | **sequence similarity** ≥ `pctseq`, only if BOTH sides are resolved | `pctseq=0.70` | **`--pctseq 0` → off** |

Genotype comparison runs but does **not** gate the match; it only fills `gt_match` (our VCFs carry no
GT, hence `gt_concordance: 0.0` in every summary — that number is meaningless for us, not a failure).

Score, used only for ranking inside the picker:
`score = (seqsim + sizesim + ovlpct) / 3 × 100`, ties broken by start- then end-distance.

### 1.3 Size filters — why 258,767 calls become 683 comparisons
- `sizemin=50` — base records smaller than this are **not counted at all** (not FN).
- `sizefilt=30` — comp records below this are dropped; comp calls between 30 and 50 that match nothing
  are **not** counted as FP either (`bench.py:390`).
- `sizemax` — default 50,000. **We pass `--sizemax 1000000000`**, otherwise every SV above 50 kb would
  silently vanish from the denominator.

This is why `hard`'s 258,767-record VCF is scored as `comp cnt = 683`: everything else is SNV/short-indel
sized. **The precision Truvari reports is SV-scale precision only** — our ~250 k small calls are neither
credited nor punished.

### 1.4 Consequences that matter for us
- **`--pctseq 0` is not a nicety, it is mandatory**: our ALT fields are symbolic (`<DEL>`), never
  resolved sequence, so a 0.70 sequence gate would reject everything.
- With sequence off, a match is decided by **TYPE + within 500 bp + 70 % size agreement**. That is
  exactly what `corelen`'s bp-exact sizing buys us and what core-count × average-step sizing loses.
- **`-f/--reference` is not needed when `pctseq=0`** (verified: the same bench runs, same numbers,
  with no `-f`). Useful, because `simulation/ref.fasta` no longer carries the ALT contigs
  (`22_GL383582v2_alt` …) that the July 12 VCFs still reference, and passing it now raises a pysam
  `KeyError` on those records.
- **A DUP is not an INS** unless `--dup-to-ins` is given, and a complex event written as several
  records competes with itself under `pick=single`.

---

## 2. The best method today

```
svp graph  →  svp refine  →  tools/corelen  →  tools/call_v2on3  →  tools/anchor_v2on3  →  truvari
```

The `_v2on3` pair is the **bp-exact** path: `corelen` re-reads the reads and recovers, per node, the
core's length in bp and, per edge, the native overlap between adjacent cores — the two numbers the v2
payload used to carry. `call_v2on3` then types and sizes every variant from the spliced bp length of a
path, instead of `call_v3`'s core-count × global average step (77.56 bp).

**That single difference is worth ~35× in precision** (`hard`, PASS: 0.010 → 0.693). Everything else in
the two paths is the same graph and the same calls (identical 199,801-call census on `complex`).

Scoring invocation (`scripts/truvari.sh`, plus a second `--passonly` run):

```bash
truvari bench -b truth.vcf.gz -c calls.vcf.gz -f ref.fa \
  --pctseq 0 --sizemax 1000000000 [--passonly] [--typeignore] -o <outdir>
```

---

## 3. Measured baseline

### 3.1 `hard` and the three translocation sims — reference record
All from July 12; **graphs built BEFORE refidx** (dual-strand, no vote — `hard/graph.bin` is 3.12 GB).
`--passonly`:

| sim | caller | TP | FP | FN | precision | recall | F1 | truth |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| hard | v3 | 449 | 42,516 | 151 | 0.010 | 0.748 | 0.021 | 600 |
| **hard** | **v2on3** | **473** | **210** | 127 | **0.693** | **0.788** | **0.737** | 600 |
| tra_intra | v3 | 471 | 42,711 | 249 | 0.011 | 0.654 | 0.021 | 720 |
| tra_intra | v2on3 | 492 | 261 | 228 | 0.653 | 0.683 | 0.668 | 720 |
| tra_inter | v3 | 483 | 42,463 | 237 | 0.011 | 0.671 | 0.022 | 720 |
| tra_inter | v2on3 | 501 | 226 | 219 | 0.689 | 0.696 | 0.692 | 720 |
| tra_inv | v3 | 473 | 42,553 | 237 | 0.011 | 0.666 | 0.022 | 710 |
| tra_inv | v2on3 | 496 | 249 | 214 | 0.666 | 0.699 | 0.682 | 710 |

`hard`, v2on3, PASS, per type: **DEL 269/288 (0.93) · INV 56/78 (0.72) · INS 58/78 (0.74) · DUP 90/156 (0.58)**.
Adding `--typeignore` changes **nothing** (473/210/0.737, identical) — every call Truvari matched on
`hard` already carried the right SVTYPE.

### 3.2 `complex` — measured 2026-07-28, first time ever
Same pipeline, run today, **graph built WITH refidx** (vote + repeat merge; 1.58 GB, K = 48.3 M → refined
43.3 M, `graph` stage 12 min 33 s at 81.6 % of 16 cores). Truth = 26,407 records, of which 1,407 are
SV-sized.

| caller | filter | TP | FP | FN | precision | recall | F1 |
|---|---|---:|---:|---:|---:|---:|---:|
| v3 | all | 286 | 7,616 | 1,121 | 0.036 | 0.203 | 0.061 |
| v3 | PASS | 286 | 6,810 | 1,121 | 0.040 | 0.203 | 0.067 |
| v2on3 | all | 302 | 2,133 | 1,105 | 0.124 | 0.215 | 0.157 |
| **v2on3** | **PASS** | **302** | **268** | 1,105 | **0.530** | **0.215** | **0.306** |
| v2on3 | PASS + typeignore | 469 | 101 | 938 | **0.823** | 0.333 | 0.474 |

Per type (v2on3, PASS): **DEL 207/260 (0.80) · INS 95/130 (0.73) · INV 0/130 · DUP 0/130 · all 19 complex
grammars 0/887.**

Two separate facts hide in that table:

1. **The complex classes are found but mistyped.** Type-ignore moves 167 calls from FP to TP and lifts
   precision to 0.82. The material is in the graph and anchored to the right place; the caller has no
   vocabulary above DEL/INS/BND to name it. → this is exactly the argument for a single **CPX** class.
2. **INV and DUP collapse to zero**, and type-ignore does not rescue them. On `hard` they were 0.72 and
   0.58. The two runs differ in **both** the sim and the graph build (pre-vote dual-strand vs voted), so
   this is a confounded measurement, not yet a verdict on the strand vote. The clean experiment —
   re-running `hard` through today's binary — has **not** been done.

---

## 4. What this baseline depends on, and what the direction change breaks

Every number above is **reference-anchored**: `anchor_v3`/`anchor_v2on3` re-extract the reference, match
graph anchor sigs to reference positions, and emit CHROM/POS/END/SVLEN so Truvari can compare intervals.
Truvari cannot score anything else — its whole gate set (refdist, size similarity, contig identity) is
defined on reference coordinates.

So if calls move to **normal-genome / LCP-core coordinates**, Truvari stops being usable as-is. The
options, in the order they cost:

1. Keep an anchor step whose only job is core-ID → reference-position, and keep scoring with Truvari
   (what happens today; the anchor is already a separate binary, not part of the caller).
2. Project the core-level call back to bp once, at report time, and score with Truvari — same thing,
   but the caller itself never sees a coordinate.
3. Score natively against the truth's **novel adjacencies** (`simulation/truth/*/adjacencies.bed`, already
   emitted by the sims) or against a core-interval BED derived from the normal graph, with our own
   bench. Then `refdist`/`pctsize` have to be redefined in core units, and the numbers become
   incomparable with everything in §3.

Whatever is chosen, §3 is the line to beat: **F1 0.737 on `hard`, 0.306 (0.474 type-blind) on `complex`.**
