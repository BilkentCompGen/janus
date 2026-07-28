# Ground-truth cores — a truth VCF re-expressed in LCP signatures

`tools/gt_cores <ref.fasta> <truth.vcf> <out.jsonl> [K=10] [threads=16]`

The graph has no coordinates. A truth VCF has nothing but coordinates. This file is the join: every truth
variant, re-expressed as the **core signatures** a caller would have to produce, so a call can be checked
against ground truth by sig instead of by position.

Built 2026-07-28 for `hard` (125,830 records) and `complex` (26,407) →
`outputs/svp/gt/{hard,complex}.jsonl` (435 MB / 106 MB). Runtime: 4.5 s and 1.3 s on 16 threads.

## Format — one JSON object per VCF record

```
id, svid, chrom, pos, end, svtype, svlen, grammar, symbol   passed through from the VCF
span: [s0, e0]                        the variant's reference span, 0-based half-open
                                      (an insertion has an EMPTY span at its insertion point)

pre  / pre_rc     [[st, en, sig], ...]   the K cores of the NORMAL sequence ending before the variant
post / post_rc    [[st, en, sig], ...]   the K cores of the NORMAL sequence starting after it
                                          st/en are reference coordinates on `chrom`
var  / var_rc     [[st, en, sig], ...]   cores of the variant's OWN sequence, coords relative to it
alt  / alt_rc     [[st, en, sig, novel]] the ALT allele rebuilt with its flanks and re-parsed;
                                          coords relative to the rebuilt window
var_truncated: 1  present when `var` had >200 cores (first/last 100 kept)
alt: null + alt_skip: "..."             the grammar is multi-block / dispersed and cannot be spliced
                                          from a single window (all 19 complex grammars)
```

**`_rc` is not decoration.** svp inserts both strands now, so every locus has a mirror twin with
different sigs and no core in common. A call may land in either family; both are given.

**`novel=1` means the sig does not occur anywhere in the reference window** — forward or reverse
family. Those are precisely the cores that should be tumour-only in the graph.

## Why `alt` exists and `var` is not enough

A deletion has no sequence of its own. What the tumour graph carries is the core that **spans its two
joined flanks**, and that core is in neither `pre`, nor `post`, nor `var`. It appears only in `alt`.
`var` answers "what sequence does this variant consist of"; `alt` answers "what does the tumour read
look like here" — the second is what the caller sees.

## What the file already shows (measured over both sets)

| type | records (hard / complex) | mean `var` cores | mean **novel** `alt` cores |
|---|---|---|---|
| SNV | 118,650 / 20,000 | 0.0 | 2.3 |
| DEL | 3,578 / 2,760 | 6.6 / 6.4 | 2.3 |
| INS | 3,368 / 2,630 | 0.8 / 1.4 | 3.4 / 4.2 |
| DUP | 156 / 130 | 149.5 / 83.0 | 2.3 |
| INV | 78 / 130 | 85.4 / 91.4 | **4.7** |
| 19 complex grammars | — / 887 | 120–140 | n/a (`alt_skip`) |

- **A SNV has no core of its own and still breaks 2.3 cores.** `var` is empty below ~152 bp because an
  L5 core is longer than the variant. That is the resolution limit, stated as a number.
- **INV gives 4.7 novel cores — twice the 2.3 of a DEL, in two separate clusters**, with the inverted
  interior scoring `novel=0` because it matches the reference's *reverse* family. The graph-level claim
  "a large inversion's interior is shared; what is callable is its two junctions" is visible directly in
  the data, without running the caller.
- **DUP's junction costs the same 2.3 cores as a DEL's**, while its `var` carries the whole duplicated
  segment — the signal is the junction and the doubled coverage, not new sequence.

## Known gap

The 887 complex-grammar records (`dupINVdup`, `delINVdel`, `rTRA`, `dDUP`, …) get `alt: null`. Their
blocks are dispersed across the genome, so the ALT allele cannot be spliced from one window. Under the
new **CPX** vocabulary these should be assembled **per `SVID`** (all of an event's records at once) into
one alt allele — that is the next thing this tool needs, and it needs `TARGET`/`INSORD` handling.

## Method notes

- Parse is byte-identical to the graph's: `extract()` at LEVEL 5 / MARGIN 10 over uppercased sequence.
  No repeat merge, no blacklist, no strand vote — all three are gone from svp.
- Cores come from a ±4,000 bp window, not a whole-genome parse. The LCP parse is locally consistent, so
  an interior core is the same either way. `selfcheck` re-cuts 200 sampled variants at double the window
  and asserts the emitted sigs are unchanged: **0/199 differ**, both sets.
- **The reference span is verified against the simulator's own tumour genome**, not asserted:
  `scripts/verify_gt_cores.py` rebuilds each ALT allele, takes the 120 bp centred on the junction, and
  looks for it in `out.fa`. **200/200 on `hard`, 200/200 on `complex`.** One base of error anywhere in
  the span would break every probe. Two failure modes belong to the verifier, not to gt_cores: its
  ref→tumour offset model ignores translocations and dispersed dups (so it falls back to scanning the
  whole contig — 71 of the 200 complex probes land that way), and a locus inside a larger inversion
  appears reverse-complemented (2 probes, found as RC).
- **This check caught a real off-by-one.** The insertion point of a symbolic `<INS>` is `POS-1`, BEFORE
  the POS base, not `POS`. The first version spliced at `POS`: 0/40 probes matched. At `POS-1`, 39/40
  (the 40th sits inside a 929 kb inversion and matches as RC). Every INS junction core in the first
  generation of these files was wrong by one base; both sets were regenerated.
- 12 of 125,830 (`hard`) and 2 of 26,407 (`complex`) records get fewer than K cores on one side — they
  sit within 10 cores of a contig end.
