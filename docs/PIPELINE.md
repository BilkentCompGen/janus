# Pipeline — the stages, and what each one does

`snape` calls somatic structural variants from **one normal and N tumour** HiFi read sets. The
somatic decision is made in a graph built from the reads themselves, with no reference and no
assembly; a reference enters only at the end, to give the surviving candidates coordinates a VCF can
carry.

    snape run --samples <sheet.tsv> <ref.fa> <dir> [--threads 16] [--level 5] [--margin 10]
              [--min-mult 3] [--supp-discover 5] [--supp-assign 3]

The sheet is a three-column TSV with a header row, one row per read set, exactly one `normal`:

    name    role    path
    N01     normal  hifi/normal.fastq
    T01     tumor   hifi/t01.fastq
    T02     tumor   hifi/t02.fastq

It is copied to `<dir>/samples.tsv` on the first run and checked row by row on a resume, so a
restart cannot silently swap a line. The names become the VCF's sample columns and the prefix on
every synthesised read name, which is how a contig still knows its line after alignment.

One command, eight stages, each skipped when its output is already in `<dir>`, so a kill and a
restart pick up where they left off. L is the number of tumour lines:

| # | stage | in | out |
|---|---|---|---|
| 1 | `graph` | all read sets | `graph.bin` (S = L+1 count planes) |
| 2 | `bubbles` | `graph.bin` + all read sets | `drop.sigs` |
| 3 | `junctions` | `graph.bin` + tumour reads | `ctg.<name>.fa`, one per line |
| 4 | minimap2 + samtools | each `ctg.<name>.fa` + `ref.fa` | `ctg.<name>.bam`, L independent jobs |
| 5 | `call` (no normal) | all contig BAMs | `candidates.vcf` |
| 6 | `pick` | all contig FASTAs + `candidates.vcf` + normal reads | `normal_at_candidates.fastq` |
| 7 | minimap2 + samtools | that FASTQ + `ref.fa` | `normal.bam` — ONE alignment, whatever L is |
| 8 | `call` (with normal) | all contig BAMs + `normal.bam` | `calls.vcf` + `calls.<name>.vcf` |

**A variant is extracted once, not once per line.** Stages 1–3 share one graph and one chunk set, so
a junction two lines carry is one chunk with one id, and which lines carry it falls out of the count
planes rather than out of intersecting coordinates afterwards. What is per line is the evidence:
each line cuts its own contigs from its own reads and aligns them on its own. The normal is aligned
once for the whole fleet, which is what keeps L tumours cheaper than L pairwise runs.

`call` counts support **per line**: discovery needs `--supp-discover` contigs from ONE line, and a
line is then genotyped present at the lower `--supp-assign` floor. With one tumour the assignment
floor cannot bite, because discovery already demanded more — so a one-tumour run produces what it
always produced.

Every stage takes `.svb` or uncompressed FASTQ directly, so nothing has to convert between them —
minimap2 is the only thing in the pipeline that cannot read `.svb`, and `run` streams `snape dump` into
it rather than writing a FASTQ out. The only outside programs are the aligner and samtools, which snape
is not and does not reimplement. Every default `run` passes is the one measured in `RESULTS.md`.

## 1. `graph` — reads → `graph.bin`

One pass over both samples. Each read is parsed into LCP cores (level 5, `lib/lcptools` `biolcp`)
**both ways**, its own bytes and their reverse complement, because the parse is not
reverse-complement symmetric (`FOUNDATIONS.md`). Cores within `margin` (10 bp) of a read end are
dropped — those are the ones a read boundary can truncate.

A core's identity is `core_sig` = **MurmurHash64A over the core's homopolymer-compressed bytes**
(`src/extract.cpp`). Compression before hashing is why a HiFi run-length slip does not change a
locus's identity; start/end stay in raw read coordinates, so positions and margins are unaffected.

Counts are exact — one `sig -> slot` map and one `(prev_sig, cur_sig) -> slot` map, 64-way sharded,
merged from per-worker maps, with the S counts of each slot in one contiguous arena rather than in
the map value. A signature becomes a node if it was seen `min_mult` (3) times in **any** sample; an
edge survives iff both of its endpoints do. Node id is the index into the sig-sorted survivor array.

```
char[4] "SNPG" | u32 version=2 | u32 S | u32 reserved   # S = samples, normal included
u64 K                                              # nodes
u64 sig[K]                                         # MurmurHash64A over the HPC core bytes
u8  cov[S][K]                                      # reads of sample s carrying node i, saturating
u64 E
u32 esrc[E], u32 edst[E]                           # directed adjacency, SoA
u8  ecov[S][E]                                     # per-sample traversal counts
```

Plane 0 is always the normal, planes 1..S-1 the tumour lines in sheet order. No colour array and no
sequence is stored. **Colour is derived**: `cov[0][i] == 0` is somatic, and *which* lines carry it
is which planes are non-zero — the support pattern, without a coordinate anywhere. That derived
somatic set is the entire somatic hypothesis; everything downstream only narrows it.

Counts saturate at 255 because no stage reads their magnitude: three of the four uses are zero
tests or a floor test, and the fourth (the chunk VAF) is behind a filter that ships disabled.
Measured on the real hcc1954 graph, clipping at 255 leaves the surviving core set and the chunk set
identical. The width is one typedef (`cov_t`, `src/snape.h`).

The magic and version exist because v1 had neither — it began with `K` — so a stale v1 file and a
v2 file were indistinguishable except by size arithmetic, and `run` resumes on nothing but the
presence of the file. A version mismatch is an abort, not silent garbage.

`graph.bin` is not byte-reproducible in its junk-core tail (map merge order), but the surviving core
set is deterministic.

## 2. `bubbles` — `graph.bin` + reads → `drop.sigs`

Most somatic material in the graph is not a junction: it is SNVs and sub-SV edits, and each one
would otherwise become a chunk, a contig, and a minute of alignment. A **bubble** is a source with
out-degree ≥ 2 whose two node-disjoint linear arms converge on one sink, with exactly one arm
somatic. An L5 core spans ~217 bp, so a core signature cannot see a 50 bp change. Two things separate
a local edit from a junction and both are needed, because they are orthogonal:

- **how much sequence changed** — the arms' length difference, in BASES. An arm of n cores spans
  `sum(core len) - overlap*(n-1)`, where the overlap of consecutive level-5 cores is measured on each
  run's own reads (122.5 bp on hcc1954, h2009, hcc1395 and colo829; 122.8 on the simulation);
- **whether what is there is the same sequence in the same order** — the arms opened up with
  minimizers (k=21, w=22, on homopolymer-compressed bytes) and compared by longest common
  subsequence, which collapses for a novel adjacency and does not for a local edit.

An arm is dropped into `drop.sigs` when `|span_S - span_G| <= max_dbp` AND `|m_S| - lcs <= 20`.
A balanced junction changes no length at all, so the span clause alone would delete it; a
substitution changes no order, so the LCS clause alone would keep it.

**`max_dbp = 50` is not a swept number: it is `call`'s own `min_sv`.** Two arms differing by less
than the smallest event we would report differ by nothing we would report. An arm containing a core
no read gave bytes for has no measurable span and is KEPT rather than judged.

This replaced a length test in MINIMIZER COUNTS (`max_dlen` 12). A minimizer is one sample per 16.6
raw bases — measured — so that was a ~200 bp test guarding a 50 bp floor, and it is where the small
insertions and deletions were going. Over four CASTLE panels the change takes insertion false
negatives from 138 to 101 and deletion from 167 to 151, with breakend and duplication false negatives
identical on every panel. It drops 45.4–57.4% of somatic cores where the old rule dropped 48.6–63.5%.
`SNAPE_BUBBLE_TSV=<path>` dumps one line per bubble with every statistic a drop rule could be built
from, which is how any candidate rule is scored offline in two seconds instead of a 4-minute re-run.

## 3. `junctions` — `graph.bin` + tumour reads → `ctg.<name>.fa`

A **chunk** is a connected component of the somatic subgraph, minus everything `drop.sigs` named:
one junction candidate, found without a reference. **The chunk set is built once, from the one
shared graph**, so a junction two lines carry is one component with one id — this is where "a
variant is extracted once, not once per line" is actually enforced. Each chunk also records its
carrier mask, the lines whose count plane is non-zero on any of its cores.

The cutting pass is then per line: each tumour line scans its OWN reads and writes its own
`ctg.<name>.fa`, and the per-chunk budget is per (chunk, line) — otherwise the first line to reach a
chunk would eat all five slots and a chunk three lines carry would yield contigs from one. For each
chunk, the first `per_chunk` (5) reads of that line carrying it are cut down to

    [ start of (first somatic core − pad_cores) , end of (last somatic core + pad_cores) ]

with `pad_cores` = 50, which is the novel sequence plus a few kilobases of germline anchor either
side. No assembly and no cross-read merge is needed to get there: consecutive LCP cores of one read
are already adjacent and overlapping in that read's own bytes, so the span between two core
boundaries **is** the sequence. Runs closer together than the padding are merged into one window —
a mobile-element insertion is somatic only at its two joins, and cutting them apart would mean no
contig ever contains the insertion.

A window whose achieved anchor is under `min_flank` (300 bp) on either side is not emitted, and the
count of those is reported: they are exactly the chunks a cross-read merge would be for. Records are
named `<line>|<read>_<n>` — the `<n>` so a read carrying two junctions cannot inflate either one's
support, and the `<line>|` prefix so a contig still knows which line it came from after alignment.
`.svb` reads have no names of their own (they are synthesised `r0, r1, ...`), so without the prefix
two lines would produce the same names and their support would silently merge.

## 4, 7. minimap2

`minimap2 -ax map-hifi -Y --secondary=no` into `samtools sort`. Stage 4 places the junction contigs
of each line as an independent job; stage 7 places only the normal reads stage 6 picked, **once for
the whole fleet** — that single normal alignment is what makes N tumours cheaper than N pairwise
runs. Both write to `<out>.tmp`, pass `samtools quickcheck`, and only then get their real name:
`samtools sort` can exit 0 on a truncated stream, and a resume must not mistake a half BAM for a
finished one.

## 5, 8. `call` — aligned contigs → breakpoints

The part of a breakpoint caller that is doing work on this input, and nothing else: a MAPQ and
mismatch gate, split alignments chained in **query** order, CIGAR indels at or above the SV floor,
clustering, a support floor. Phasing, coverage histograms, VAF and multi-sample subtraction are idle
here by construction — a contig is not a read pileup, and the somatic decision was already made in
the graph.

All lines' contig BAMs are read into ONE junction list, and the clustering rule is untouched: it
does not know lines exist. What is per line is the support tally — a cluster's distinct contigs are
counted into one bucket per line. Discovery needs `--supp-discover` (5) from a single line, so five
contigs spread over five lines is not a call; a line is then genotyped present at `--supp-assign`
(3). Output is one `calls.vcf` with a column per line plus a `calls.<name>.vcf` projection each. The
germline check runs once, on the one normal BAM, over the combined call set.

A junction is two reference sides plus the direction each side leaves the reference in: `d=+1` the
alignment ends there going right, `d=-1` it starts there going left. That pair is the whole type
vocabulary — (+,−) deletion, (−,+) duplication, equal signs inversion, two contigs a breakend — so
the caller never has to know what an SV is in order to type one.

Three cases do not reach the aligner as a split, and each is read out of the contig's own bases:

- **short tandem duplication** — minimap2 writes the second copy as an insertion. The copy it
  duplicates is right there in the contig, on one side or the other of the inserted bases, and those
  flanking bases are aligned, so they *are* the reference at that locus. Which side matches says
  which copy the aligner kept, and that fixes the duplicated interval.
- **short inversion** — the aligner walks straight through it and pays with a dense run of tiny
  indels. A window that shredded is a breakpoint pair whatever the aligner called the pieces.
- **a VNTR that looks like one** — the same shred. It is separated from a real inversion by the
  tandem-repeat fraction of the window's own bases (`shred_rep`, 0.40); random sequence sits near
  0.25.
- **an insertion the aligner split instead of writing as a CIGAR I** — most short insertions are
  another copy of sequence the genome already has, so minimap2 may place that copy at its source
  locus, leaving three segments whose OUTER two are adjacent in the reference. Chaining pairwise
  would emit two breakends at the source locus and the event would score missed at its own
  coordinates. When two segments of one contig resume the reference within `min_sv` of where the
  previous one left off, same chromosome and strand, with at least `min_sv` of query sequence
  between them, that sequence is emitted as an INSERTION at the join. The length is every query base
  between the two anchors, so a middle segment the MAPQ gate already dropped is handled by the same
  rule. Measured: fires on 128 of h2009's 537,144 segments, worth +3 TP and 0 FP.

Support counts **distinct contig names**, so a contig reaching the same junction twice supports it
once. Positions are 1-based from the moment a junction is built, so the VCF writer converts nothing.

**The mismatch gate** (`max_mm`, 0.01) drops an alignment carrying far more mismatch than the reads
can contain — that is not this locus. It counts substitutions only: NM includes inserted and deleted
bases, so using it raw would score a contig carrying a real 400 bp deletion at 8% and throw away
exactly the events being looked for.

**The anchor-repeat gate** (`anchor_rep`, 0.50) drops a call whose anchor is a tandem repeat, taken
on the anchor side only, where the bases are aligned at under 1% mismatch and therefore are the
reference at that locus.

**The germline question** is the one thing the graph cannot ask. The graph decides somatic by exact
core sequence, and in a repeat that is not the same question as somatic: the tumour and normal reads
simply do not spell the locus with the same 217 bp core, so `ncov == 0` is true of the *string*
while the *variant* is shared. `call` asks the normal reads at the reference locus instead — same
event, or a locus the normal cannot place either (a quarter or more of its reads above 2% mismatch),
and the call is dropped. Breakends are asked at both ends.

Stage 5 runs this with no normal, to produce the candidate loci. Stage 8 runs it again with the
normal BAM stage 7 built.

## 6. `pick` — which normal reads are worth aligning

The germline question is the largest single source of precision in the pipeline: on the shipped build
at the shipped floor, without it hcc1954 leaves 744 false calls instead of 112, and precision falls
from 0.894 to 0.560 (measured 2026-08-13; the 1409/241 this line used to quote was a pre-fix build at
a different floor). Answering it by aligning the whole normal would make the
aligner the dominant cost of a caller whose premise is that the graph replaces alignment — and it
aligns far more than the question needs, which is only ever asked at the candidate loci.

Each candidate names the contig it came from (`CTG=` in the VCF). `pick` takes those contigs'
**germline anchors** — the flanks, which are sequence the normal has too, as opposed to the somatic
middle, which it does not — parses them into cores, and emits every normal read sharing at least
`min_shared` (5) of them. On hcc1954 that is about 1.8% of the sample, and it answers the question as
well as all of it: 0.887/0.858 against 0.891/0.858 (both arms at a floor of four, where the
comparison was made, and at the 1.56% the pre-fix build picked), and 16 minutes of alignment against
about two hours.

## Source layout

    src/samples.cpp     the sample sheet: parse, validate, lock to <dir>/samples.tsv
    src/graph.cpp       the graph stage, the .svb/FASTQ reader every stage shares, checkpoint IO
    src/extract.cpp     biolcp adapter: homopolymer compression + MurmurHash64A core_sig
    src/bubbles.cpp     sub-SV bubbles -> drop.sigs
    src/junctions.cpp   somatic chunks -> junction contigs
    src/call.cpp        aligned contigs (+ the normal at the locus) -> VCF
    src/pick.cpp        the normal reads the germline question needs
    src/svb_pack.cpp    FASTQ -> 2-bit .svb
    src/svb_dump.cpp    .svb -> FASTQ
    src/main.cpp        the CLI and `run`
    src/snape.h           shared structs, the reads interface, checkpoint IO decls
    src/unordered_dense.h   vendored flat hash map (ankerl)

One binary, one Makefile, no CMake, no CI.

`-DLARGE_GENOME` is load-bearing: it must match the flag `lib/lcptools` was built with (64-bit core
start/end), or core structs read back corrupt.
