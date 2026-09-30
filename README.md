# snape

Somatic structural variants from a tumor/normal pair of HiFi read sets. Every read is cut into LCP cores, the cores of
all samples become one graph, combinations of graph edges that tumor reads carry and normal reads never do become
sites, and a short piece of one tumor read around each site is aligned to the reference and turned into calls.

## Build

The vendored LCP library needs its `config.h` and `liblcptools.a` in `lib/lcptools` (neither is tracked). Its own
default target installs into /usr/local, so build the two files directly, once:

```
cd lib/lcptools && make LABEL=0 POS=64 config.h encoding_s.o core_s.o lps_s.o && ar rcs liblcptools.a *_s.o && cd ../..
make -j16                                # build/snape
```

g++ (C++17), zlib and pthreads; `-Wall -Wextra -Wpedantic` builds without warnings.

## Running

```
snape pack -i reads.fastq.gz -o n.svb                  # once per sample (plain or gzipped FASTQ)
snape refidx -i hg38_parmask.fa -o hg38_parmask.cidx -t 32   # once per reference
minimap2 -x sr -t 16 -d hg38_parmask.sr.mmi hg38_parmask.fa       # once per reference: short-read index for the gap stretches
./pipeline.sh <dir> <normal.svb> <tumor.svb> [threads]      # default 16 threads
```

`pipeline.sh` runs the stages below into `<dir>` (reads are linked in as `n.svb`, `t.svb`), logs each one to
`log.<step>` with its time and memory in `steps.txt`, and leaves `done.<step>` behind so a rerun continues where it
stopped. The binary and the reference directory are the `B=` and `R=` paths at its top (`build/snape` and `ref/` of
the main checkout); `ref/` holds `hg38_parmask.fa` with `.fai`, its short-read minimap2 index `hg38_parmask.sr.mmi` and the core
index `hg38_parmask.cidx`, and `minimap2` must be on the path. It ends with `ev.paf`; `snape call` then makes the
VCF (see bench/ for the call and evaluation harness):

```
snape call -i ev.paf -r ev.tsv -d hg38_parmask.fa -g graph.gph -o calls.vcf -v vntr.trf.bed
```

Options: `-i` input, `-o` output, `-r` reads (comma-separated, normal first), core index, pieces or ev.tsv, `-d` sites
or reference FASTA, `-t` threads (at most 64), `-g` graph, `-v` VNTR bed.

| step | command | reads | writes |
|---|---|---|---|
| extract | `extract -i n.svb -o n.svb.exc` (and t) | svb | exc |
| count | `count -i n.svb -o n.cnt` (and t) | svb, exc | cnt |
| graph | `graph -i n.cnt,t.cnt -o graph.gph` | cnt of every sample | gph |
| prune | `prune -i graph.gph -r n.svb,t.svb -o prune` | gph, svb + exc of every sample | prune.drp, prune.sites |
| select | `select -i graph.gph -d prune.sites -r n.svb,t.svb -o ev` | gph, sites, svb + idx | ev.fa, ev.tsv |
| place | `place -i ev.fa -r hg38_parmask.cidx -d hg38_parmask.fa -o ev` | pieces, cidx, reference | ev.win.fa, ev.on.fa, ev.off.fa |
| local | `minimap2 -cx map-hifi --secondary=no ev.win.fa ev.on.fa` | windows, placed pieces | local.paf |
| lift | `lift -i local.paf -d hg38_parmask.fa -o aln.paf` | local.paf, .fai | aln.paf |
| clip | `clip -i ev -r aln.paf -o ev.gaps.fa` | pieces, aln.paf | ev.gaps.fa |
| gaps | `minimap2 -cx sr --secondary=no hg38_parmask.sr.mmi ev.gaps.fa` | stretches, short-read index | ev.gaps.raw.paf |
| patch | `patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf` | pieces, both alignment sets | ev.paf |
| call | `call -i ev -r n.svb,t.svb -d hg38_parmask.fa -g graph.gph -o calls.vcf -v vntr.trf.bed` | ev.*, svb + idx + exc, gph, reference, TRF bed | calls.vcf |

## Stages

**pack** packs a FASTQ into svb. A read with any base other than A, C, G, T (either case) is stored raw.

**extract** homopolymer-compresses each read (one base per run), parses both strands of the compressed sequence with
LCP up to level `LEV` and writes every core's interval in read coordinates to the exc file. Forward cores are on the
read, reverse cores on its reverse complement.

**count** hashes every core (64-bit hash of its compressed bases), skipping cores within `MARGIN` bases of a read end,
and counts them in one sharded table; consecutive distinct cores along a strand are counted as pairs in a second table.
Cores seen at least `MM` times are the nodes; a pair of two nodes is an edge.

**graph** unites the samples' nodes and edges (the first sample is the normal). Each node keeps every sample's count.
An edge is novel when the normal has no read on it and a later sample has at least `MM`. Nodes are renumbered into
segments, the maximal non-branching paths (u to v joins them when u has one successor and v one predecessor), so a
segment's nodes are consecutive and its internal edges are implicit.

**prune** finds sequencing errors and somatic sites.
- Tips: from a node with two or more successors (and, in a second pass, predecessors), a chain of at most `ML` nodes
  that all have the same sample set, are not in every sample, have one predecessor each, and ends in a dead end.
- Bubbles: from a node with two or more successors, each branch is traced through single-successor nodes, up to
  `ML + 1` nodes, to a node with two or more predecessors. Two branches that meet at the same node form a bubble; an
  arm is a variant arm when all its nodes have one sample set, not every sample, that no node of the other arm has.
  One pass over all reads takes each arm's bases from the read that carries its exact path; when the two arms' bases
  are fewer than `X` edits apart, the variant arm's inner nodes are sequencing errors.
- Combinations: in the same read pass, every read chain names its combinations, a pair (edge in, edge out) through a
  branching node or across a whole segment of two or more nodes, and counts them per sample in one table, keeping for
  each sample the read on which the combination is most central (largest smaller flank; smallest read index on ties).
  Only combinations at a node or segment touching a novel edge are counted.
- A combination is a site when the normal has no read with it, a tumor has at least `MM`, it touches a novel edge
  whose ends are not dropped, and none of its edges' nodes are dropped.

**select** reads each site's best read again (random access through the idx), finds the combination on it and up to
`NA` anchors on each side: nodes on the read that occur once in its chain, have reads in the normal and in the site's
tumor, and at most `MULT` times that sample's median node count (median over nodes that every sample has). Anchors are
scored by `sqrt(min(1, normal/median) * min(1, tumor/median))` times a weight that falls from 1 at the combination to 0
at `DL` bases, and listed innermost first. The site runs from its innermost left anchor to its innermost right anchor
(the combination's own cores on a side without anchors); its interval runs from the outermost anchors with a flank of
`max(PAD, FL * site length)` on each side. Overlapping intervals on one read merge into one piece, which is cut from
the read.

**refidx** cuts the reference into N-free windows of `WIN` bases, parses each with up to `MAR` bases of context on each
side and keeps the cores whose forward start lies in the window.

**place** looks up every core of a piece in the core index. Hits with exactly one reference match are grouped by
chromosome, strand and diagonal (neighbouring diagonals at most `TOL` apart); a group of at least `MINC` hits places the
piece and gives a reference window from its first to last hit, widened by the piece length plus `MAR` on each side.
Overlapping windows on a chromosome merge.

**lift** moves every local alignment back to chromosome coordinates (`aln.paf`). There is no whole-genome alignment:
a piece's unaligned sides, including the partner side of a junction that lies in repeated sequence and therefore
found no window, are the clip stage's job.

**clip** and **patch** (src/gaps.cpp) align the pieces' unaligned stretches on their own. A stretch of `MING` to `MAXG`
bases between two alignments or at a clipped end is too short to survive as a supplementary alignment of a 15 kb piece,
but on its own, in minimap2's short-read mode against `hg38_parmask.sr.mmi`, it places the partner breakpoint of a
complex junction. clip writes the stretches as `ev.gaps.fa` (named `piece|start|end`), patch writes the pieces' own alignments (`aln.paf`) and then
theirs, back in piece coordinates, as `ev.paf`.

**call** (src/call.cpp) reads the sites, the pieces and their alignments, the graph and the packed reads, makes its own
pass over every read of both samples, and writes `calls.vcf` with AA3's 36 evidence features per record (`F=`) for the
verdict. The pass fetches each site's tumor read, measures the distance between the site's innermost anchors on it and
counts the read's nodes between them (tumor-only, branch); then every read is streamed once, and each read carrying a
left and a right anchor of a site gives its distance, so every site gets per-sample spanning reads, their length
differences from the tumor allele (`D0`, `D1`), reads within `DSAME` of it, and the same over consistent reads (their
other anchors at the tumor read's spacing). Records come from the alignments: gaps of `MS` bases or more inside an
alignment, and the jumps between consecutive alignments (DEL, INS, DUP, INV, BND), attached to the site with most tumor
reads within `ROW` bases in the piece, merged within `TOL` bases; insertions of `NEAR` bases or more copying the adjacent
reference become tandem duplications. Each record's features are AA3's: tumor reads and their share of the depth, pieces,
tumor length clusters, spanning normal reads and those at the tumor length, the normal's length-difference spread,
tumor-only and branch node density, mapping quality and divergence, reference repeat content in `WIN` bases, the TRF
tandem-repeat bed, type and length. `ERR` is AA3's logistic error probability with weights `W0`; FILTER is PASS below
the `PCUT` cut, DUPL for the less likely of two records within `DUPT` bases. `bench/verdict.py` (boosted trees, held out
per pair) makes the reported PASS from `F=`.

## Constants

| constant | value | file | meaning |
|---|---|---|---|
| LEV | 5 | extract.h | LCP level of the cores |
| W | LEV + 4 = 9 | extract.cpp | bits of a core's start delta and length in exc |
| MM | 3 | count.cpp | occurrences a core needs to become a node |
| MARGIN | 10 | count.cpp | cores closer than this to a read end are not counted |
| NS | 65536 | count.cpp | shards of the core table (top 16 hash bits) |
| ZDIV, ZADD, ZPAD | 24, 120000000, 64 | count.cpp | core table size: total cores / ZDIV + ZADD slots, plus ZPAD per shard; the pair table is 5/4 of it |
| MM | 3 | graph.cpp | tumor reads an edge the normal lacks needs to be novel |
| WK | 16 | graph.cpp | segment walks run side by side |
| ML | 4 | prune.cpp | longest tip; bubble arms are traced up to ML + 1 nodes |
| MM | 3 | prune.cpp | tumor reads a combination the normal lacks needs to be a site |
| X | 20 | reads.cpp | edits at which two bubble arms stop being the same allele |
| ZF | 16 | keys.cpp | combination table slots per possible combination (power of two) |
| PAD | 1000 | select.cpp | least flank on each side of a site |
| FL | 3 | select.cpp | the flank is at least FL times the site length |
| NA | 3 | select.cpp | anchors per side |
| MULT | 2 | select.cpp | an anchor's count is at most MULT times the sample median |
| DL | 15000 | select.cpp | distance at which an anchor's weight reaches zero |
| WIN | 40000 | place.cpp | reference bases per refidx window |
| MAR | 2000 | place.cpp | context around a refidx window; margin around a placement |
| TOL | 500 | place.cpp | largest diagonal step inside one placement |
| MINC | 3 | place.cpp | unique core hits a placement needs |
| EXPL | 100 | lift.cpp | uncovered bases that send a piece to genome-wide alignment |
| IBUF | 64 MiB | io.cpp | stream buffer (doubled while one record does not fit) |
| NONE | 2^32 - 1 | io.h | no node, edge or slot |
| PF | 16 | io.h | prefetch distance |

## File formats

Binary files start with a 32-byte header `{char magic[4]; uint32 ver; uint64 n; uint64 total; uint8 resv[8]}`, all
little-endian; the magic is `SVB` plus one letter. Hashes are the 64-bit core hashes of count; node and edge indices are
positions in the graph's node list and in its edge list sorted by (source, target).

| file | magic, ver | n | total | resv | body |
|---|---|---|---|---|---|
| svb | SVB1, 1 | reads | bases | | per read `uint32 tag = length \| raw << 31`, then `length` ASCII bytes if raw, else `(length + 3) / 4` bytes of 2-bit bases (A 0, C 1, G 2, T 3, first base in the high bits) |
| svb.idx | SVBI, 1 | reads | longest read | | `n + 1` uint64 record offsets (the first is 32, the last the file size) |
| exc | SVBX, 3 | reads | cores | LEV, W | per read, in svb order: 6 bytes big-endian (cores, forward cores, record bytes including these 6), then per core, MSB first, W bits of start minus the previous core's start (0 before each strand's first core) and W bits of length; W ones escape to 32-bit start delta and 32-bit length |
| cnt | SVBC, 1 | nodes K | edges E | rb, db, MM, MARGIN | K uint64 hashes (ascending), K uint16 counts, E uint16 edge read counts, then a bit stream: K out-degrees of db bits, E edge targets of rb bits (edges by source, then target); counts saturate at 65535 |
| gph | SVBG, 2 | nodes K | external edges X | samples S | uint64 U segments, K uint64 hashes in segment order, S x K uint16 counts (sample-major), U + 1 uint32 segment starts, X uint64 external edges `source << 32 \| target` (sorted), X uint8 novel flags, K uint8 novel flags of each node's internal edge to the next node of its segment |
| drp | SVBP, 1 | dropped nodes | bubbles examined | | n uint64 hashes of the dropped nodes, ascending |
| sites | SVBS, 1 | sites | samples S | | per site 2 + S uint64: read `sample << 40 \| read index`, key `(edge in + 1) << 32 \| edge out`, reads of each sample with the combination; sorted by (read, key) |
| cidx | SVBR, 1 | entries | sequences | | 16-byte entries `{uint64 hash; uint32 pos; uint32 seq << 1 \| reverse}`, sorted by (hash, seq and strand, pos) |

Text outputs:

- `ev.fa`: the pieces, named `s<sample>_<read index>_<k>` (the k-th piece of that read), bases on the read's forward
  strand.
- `ev.tsv`, one row per site with a piece: piece name; site start and end in the piece (nearest anchors, or the
  combination itself on a side without anchors); site index in prune.sites; reads of each sample with the combination
  (comma-separated, normal first); combination key; tumor sample; left and right anchors as `node:start` in the piece,
  innermost first (`-` when none).
- `ev.win.fa`: reference windows named `<chromosome>:<start>` (0-based). `ev.on.fa`, `ev.off.fa`: pieces with and
  without a placement.
- `lift.paf`: minimap2 PAF of the placed pieces with the target renamed to the chromosome, its length from the `.fai`
  and its coordinates shifted by the window start. `aln.paf`: the local alignments in chromosome coordinates.
