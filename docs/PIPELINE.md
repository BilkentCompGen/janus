# Pipeline — the stages, and what each one does

`svp` calls somatic structural variants from a tumour/normal pair of HiFi read sets. The somatic
decision is made in a graph built from the reads themselves, with no reference and no assembly; a
reference enters only at the end, to give the surviving candidates coordinates a VCF can carry.

    svp run <normal.svb|fq> <tumour.svb|fq> <ref.fa> <dir> [--threads 16] [--level 5] [--margin 10] [--min-mult 3]

One command, eight stages, each skipped when its output is already in `<dir>`, so a kill and a
restart pick up where they left off:

| # | stage | in | out |
|---|---|---|---|
| 1 | `graph` | both read sets | `graph.bin` |
| 2 | `bubbles` | `graph.bin` + both read sets | `drop.sigs` |
| 3 | `junctions` | `graph.bin` + tumour reads | `ctg.fa` |
| 4 | minimap2 + samtools | `ctg.fa` + `ref.fa` | `ctg.bam` |
| 5 | `call` (no normal) | `ctg.bam` | `candidates.vcf` |
| 6 | `pick` | `ctg.fa` + `candidates.vcf` + normal reads | `normal_at_candidates.fastq` |
| 7 | minimap2 + samtools | that FASTQ + `ref.fa` | `normal.bam` |
| 8 | `call` (with normal) | `ctg.bam` + `normal.bam` | `calls.vcf` |

Every stage takes `.svb` or uncompressed FASTQ directly, so nothing has to convert between them —
minimap2 is the only thing in the pipeline that cannot read `.svb`, and `run` streams `svp dump` into
it rather than writing a FASTQ out. The only outside programs are the aligner and samtools, which svp
is not and does not reimplement. Every default `run` passes is the one measured in `RESULTS.md`.

## 1. `graph` — reads → `graph.bin`

One pass over both samples. Each read is parsed into LCP cores (level 5, `lib/lcptools` `biolcp`)
**both ways**, its own bytes and their reverse complement, because the parse is not
reverse-complement symmetric (`FOUNDATIONS.md`). Cores within `margin` (10 bp) of a read end are
dropped — those are the ones a read boundary can truncate.

A core's identity is `core_sig` = **MurmurHash64A over the core's homopolymer-compressed bytes**
(`src/extract.cpp`). Compression before hashing is why a HiFi run-length slip does not change a
locus's identity; start/end stay in raw read coordinates, so positions and margins are unaffected.

Counts are exact — one `sig -> (ncov, tcov)` map and one `(prev_sig, cur_sig) -> (ncov, tcov)` map,
64-way sharded, merged from per-worker maps. A signature becomes a node if it was seen `min_mult`
(3) times in **either** colour; an edge survives iff both of its endpoints do. Node id is the index
into the sig-sorted survivor array.

```
u64 K                                              # nodes
u64 sig[K]                                         # MurmurHash64A over the HPC core bytes
u32 ncov[K], u32 tcov[K]                           # normal / tumour read counts on the core
u64 E
CEdge{u32 src, u32 dst, u32 ncov, u32 tcov}[E]     # directed adjacency, per-colour traversal counts
```

No colour array and no sequence is stored. **Colour is derived**: `ncov == 0 && tcov > 0` is
somatic, `tcov == 0 && ncov > 0` normal-only, both non-zero shared. That derived somatic set is the
entire somatic hypothesis — everything downstream only narrows it.

`graph.bin` is not byte-reproducible in its junk-core tail (map merge order), but the surviving core
set is deterministic.

## 2. `bubbles` — `graph.bin` + reads → `drop.sigs`

Most somatic material in the graph is not a junction: it is SNVs and sub-SV edits, and each one
would otherwise become a chunk, a contig, and a minute of alignment. A **bubble** is a source with
out-degree ≥ 2 whose two node-disjoint linear arms converge on one sink, with exactly one arm
somatic. An L5 core spans ~217 bp, so a core signature cannot see a 50 bp change — the arms are
opened up with minimizers (k=21, w=22, ~13 per core, taken on homopolymer-compressed bytes) and
compared **in order**:

- same sequence with a local change → the arms are the same length and the longest common
  subsequence covers nearly all of the somatic arm → an edit, and its cores go into `drop.sigs`;
- a novel adjacency → shared minimizers do not stay in order and the LCS collapses → kept.

Thresholds: `max_dlen` 12 minimizers of arm-length difference, `max_unmatched` 20 somatic minimizers
outside the LCS. On `hard` that drops 83% of somatic cores and 84% of the contigs for 2 junctions of
recall.

## 3. `junctions` — `graph.bin` + tumour reads → `ctg.fa`

A **chunk** is a connected component of the somatic subgraph, minus everything `drop.sigs` named:
one junction candidate, found without a reference. The tumour reads are scanned once and, for each
chunk, the first `per_chunk` (5) reads carrying it are cut down to

    [ start of (first somatic core − pad_cores) , end of (last somatic core + pad_cores) ]

with `pad_cores` = 50, which is the novel sequence plus a few kilobases of germline anchor either
side. No assembly and no cross-read merge is needed to get there: consecutive LCP cores of one read
are already adjacent and overlapping in that read's own bytes, so the span between two core
boundaries **is** the sequence. Runs closer together than the padding are merged into one window —
a mobile-element insertion is somatic only at its two joins, and cutting them apart would mean no
contig ever contains the insertion.

A window whose achieved anchor is under `min_flank` (300 bp) on either side is not emitted, and the
count of those is reported: they are exactly the chunks a cross-read merge would be for. Records are
named `<read>_<n>`, so a read carrying two junctions cannot inflate either one's support.

## 4, 7. minimap2

`minimap2 -ax map-hifi -Y --secondary=no` into `samtools sort`. Stage 4 places the junction contigs;
stage 7 places only the normal reads stage 6 picked.

## 5, 8. `call` — aligned contigs → breakpoints

The part of a breakpoint caller that is doing work on this input, and nothing else: a MAPQ and
mismatch gate, split alignments chained in **query** order, CIGAR indels at or above the SV floor,
clustering, a support floor. Phasing, coverage histograms, VAF and multi-sample subtraction are idle
here by construction — a contig is not a read pileup, and the somatic decision was already made in
the graph.

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

The germline question is the largest single source of precision in the pipeline: without it hcc1954
leaves 1409 false calls instead of 241. Answering it by aligning the whole normal would make the
aligner the dominant cost of a caller whose premise is that the graph replaces alignment — and it
aligns far more than the question needs, which is only ever asked at the candidate loci.

Each candidate names the contig it came from (`CTG=` in the VCF). `pick` takes those contigs'
**germline anchors** — the flanks, which are sequence the normal has too, as opposed to the somatic
middle, which it does not — parses them into cores, and emits every normal read sharing at least
`min_shared` (5) of them. On hcc1954 that is 1.56% of the sample, and it answers the question as
well as all of it: 0.887/0.858 against 0.891/0.858, fourteen minutes of alignment against two hours.

## Source layout

    src/graph.cpp       the graph stage, the .svb/FASTQ reader every stage shares, checkpoint IO
    src/extract.cpp     biolcp adapter: homopolymer compression + MurmurHash64A core_sig
    src/bubbles.cpp     sub-SV bubbles -> drop.sigs
    src/junctions.cpp   somatic chunks -> junction contigs
    src/call.cpp        aligned contigs (+ the normal at the locus) -> VCF
    src/pick.cpp        the normal reads the germline question needs
    src/svb_pack.cpp    FASTQ -> 2-bit .svb
    src/svb_dump.cpp    .svb -> FASTQ
    src/main.cpp        the CLI and `run`
    src/svp.h           shared structs, the reads interface, checkpoint IO decls
    src/unordered_dense.h   vendored flat hash map (ankerl)

One binary, one Makefile, no CMake, no CI.

`-DLARGE_GENOME` is load-bearing: it must match the flag `lib/lcptools` was built with (64-bit core
start/end), or core structs read back corrupt.
