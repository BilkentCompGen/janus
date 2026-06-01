# Variant calling from bubble candidates: research synthesis

Research goal: design the next stage of svp — turning **bubble candidates** (connected
components of an LCP-core read graph, currently emitted as `bubbles.tsv` with coverage but
no positions / no allele sequences) into accurate variant calls (SNP/INS/DEL/INV/DUP/BND
and complex events), without resorting to read alignment or assembly.

The old svp on `backup` had a calling layer (`src/sv.cpp`); the user said it was "mistaking
lots of variants" and asked to redesign it. This doc audits what went wrong, audits how
recent long-read SV callers (Severus, Sniffles2, SAVANA, cuteSV) do the job, and proposes a
concrete bubble-native architecture.

## 1. Old SVP postmortem (`backup` branch, `src/sv.{h,cpp}` ≈ 162 LOC)

The old caller was a one-shot iteration over `LCPGraph::CompressedNode`s:

```cpp
// sv.h
enum class SVType { DELETION, INSERTION, DUPLICATION };
struct SVCall { SVType type; uint64_t node_hash; std::vector<uint64_t> path;
                std::vector<size_t> reference_indices; std::string chrom;
                std::vector<uint64_t> genomic_starts, genomic_ends;
                std::vector<int> tumor_samples; uint32_t normal_count, tumor_count; };

// sv.cpp
for (const auto& [hash, node] : graph.getCompressedNodes()) {
    if (isInNormal(node)) detectDeletions(...);  // emits SVType::DELETION
}
```

**Concrete bugs that map onto the user's "mistaking lots of variants":**

| # | bug | effect |
|---|---|---|
| 1 | `detectDeletions` is the ONLY emitter. The enum has `INSERTION` and `DUPLICATION` but no function ever produces them. | Tumor-only nodes never emit calls. Recall collapse on INS/DUP. |
| 2 | Per-node emission with no clustering. Every normal-only compressed node becomes a DEL. | Massive FP inflation; one true deletion ≈ many overlapping core nodes → many duplicate calls. |
| 3 | `LinearReference::align(path)` is **first-hit only**: takes `path[0]`, looks up its positions, returns those (or falls through to a scattered `pos-i` scan over all positions). No chain consistency, no repeat disambiguation. | Repetitive regions get arbitrary positions. False coordinates. |
| 4 | No INV, no BND, no complex SVs. | Inversions and translocations silently dropped. |
| 5 | No filtering / QC. | Low-multiplicity sequencing-error nodes that survive the prune still emit calls. |
| 6 | No allele sequences in output (REF/ALT are absent). | Cannot benchmark against truvari/svdb without lossy synthesis. |
| 7 | No genotyping. | No GT/AF; nothing actionable downstream. |
| 8 | Loops over `compressed_nodes` (a hash map) — non-deterministic call order. | Run-to-run instability when piped. |

These are all design-level, not implementation bugs. The new caller has to start from a
different shape.

## 2. Severus (Keskus et al. 2024) — mechanism extract

Files audited: `severus/breakpoint_finder.py` (2 874 lines), `severus/build_graph.py`
(674 lines), `severus/vcf_output.py` (523 lines), `severus/bam_processing.py`,
`severus/main.py`, `severus/resolve_vntr.py`. Copies under
`outputs/svp/research/severus_*.py`.

**Signature primitive (`breakpoint_finder.py:34–112`).** A `Breakpoint` is `(ref_id,
position, dir_1, connections[], spanning_reads[genome][hap], qual, CI)`. A `DoubleBreak`
pairs two breakpoints into a candidate SV, carrying `bp_1`, `bp_2`, directionality,
`supp` (read support count), `hap`, and per-filter flags.

**Breakend graph (`build_graph.py:22–73`).** Nodes are **genomic segments** (maximal
contiguous ref intervals spanned by reads supporting an SV), not individual breakpoints.
Three edge types: `adjacency` (a DoubleBreak; carries read support + genotype),
`ref_adj` (consecutive-segment sanity edges), `cluster_conn` (added after clustering, to
link nodes of the same complex event).

**Clustering = connected components + indel heuristic
(`build_graph.py:306–384 cluster_adjacencies`).** Step 1: `nx.connected_components(G)`.
Step 2: per-component, classify simple / complex / indel. Step 3: complex iff ≥2 distinct
SV types in the component AND inter-breakpoint spacing ≤ 1.5 × median SV length (line 376).
Otherwise reclassified as simple / tandem to suppress false complex calls on tandem
repeats. **No Louvain, no agglomeration** — pure CC + post-filter.

**Complex-SV typing (`vcf_output.py:248–433 db_2_vcf`).** Templated insertions are
detected by graph topology: an insertion breakpoint sandwiched between two deletion
breakpoints whose read-support sets intersect (line 723: `[cn for cn in conn_1 if cn in
conn_2]`) and whose insertion length matches the gap between the deletion endpoints
(line 732, `np.std(has_ins) < MAX_STD`). DUP→INS collapse (lines 1380–1420). Chromothripsis
proxy: ≥3 junctions on same chromosome with alternating orientation (line 347), connected
via single-copy segments.

**Somatic filtering (`main.py:155–158`, `breakpoint_finder.py:738–755`).**
`is_pass='FAIL_LOWCOV_NORMAL'` iff zero supporting reads in normal AND normal spanning
reads at the breakpoint ≤ `COV_THR=3`. Strand-balance check: spanning reads on opposite
strand > 3× supporting reads ⇒ artifact.

**Genotyping (`vcf_output.py:146–246 call_genotype`).** Bayesian three-genotype model
with fixed priors `(0.90, 0.09, 0.01)`. Likelihood `P(data | GT) =
Π binomial(read | allele)` with allele probabilities `(1−ε, 0.5, ε)` for hom-ref / het /
hom-alt. Emit `GT`, `GQ = 10 log10(L_second / L_best)`. Phased GT (`1|0` / `0|1`) when
read carries HP tag. AF = `DV / (DV + DR)`, where `DV` = unique supporting reads,
`DR = max(0, spanning_reads - supp)`.

**VCF format.** Simple SVs as single record (`SVTYPE=DEL|INS|INV|DUP`, `END`, `SVLEN`,
`STRANDS`). Complex events as paired BND records with `CLUSTERID=severus_N` for
downstream Truvari/svdb to recognise. Critical INFO fields: `PRECISE`/`IMPRECISE`, `MAPQ`,
`DETAILED_TYPE=Templated_ins|reciprocal_inversion|…`, `INSIDE_VNTR=TRUE`, `SUPP_READS`,
`REF_READS`. FORMAT: `GT:VAF:hVAF:DR:DV` (`hVAF` = per-haplotype VAF triple).

**VNTR handling (`resolve_vntr.py`).** Per-read realignment of breakpoint coordinates to
the nearest tandem-repeat boundary (lines 108–170), and per-locus merge of insertion +
flanking-deletion calls when both straddle the same VNTR with compatible lengths (line
281). Annotate the call with `INSIDE_VNTR=TRUE` so downstream filters know confidence is
reduced.

## 3. Sniffles2 (Smolka et al. 2024) — mechanism extract

Files audited: `cluster.py`, `sv.py`, `consensus.py`, `leadprov.py`, `genotyping.py`,
`postprocessing.py` (3 139 lines total). Copies under `outputs/svp/research/sniffles_*.py`.

**Signature primitive: `Lead` (`leadprov.py:35–55`).** Per-read SV observation. Fields
`read_id`, `contig`, `ref_start/end`, `qry_start/end`, `strand`, `mapq`, `nm`, `source`
(`INLINE` for CIGAR indels vs `SPLIT_PRIM`/`SPLIT_SUP` for supplementary alignments),
`svtype` (INS/DEL/INV/BND/DUP/SINGLE_LEFT/SINGLE_RIGHT), `svlen`, cached `seq` (insertion
bases, capped to `dev_seq_cache_maxlen`), `bnd_info`, `hap`, `phase_set`.

**Clustering: Fast Adaptive Scanning Clustering (`cluster.py`).** Two stages:

1. **Binning** (219–276): bin Leads by reference position into `cluster_binsize` (default
   1 000 bp); seed cluster from each bin with ≥ `dev_min_leads_cluster` reads. Compute
   per-seed `mean_svlen` and `stdev_start` from a sampled subset.
2. **Hierarchical merge** (277–308): sort clusters by start; merge two if
   `inner_dist ≤ min(stdev_start_A, stdev_start_B) * cluster_r` (default `r=1.0`). **TR-aware
   relaxation** (line 295): inside an annotated TR, additionally allow merging if
   `outer_dist ≤ (|svlen_A|+|svlen_B|) * cluster_repeat_h` (default `0.5`). **BND merge**:
   allow if `inner_dist ≤ cluster_merge_bnd`. Then **re-split** (125–161, 337–353) by SV
   length / mate position so distinct events accidentally fused get separated.

**SVCall (`sv.py:87–130`).** `contig`, `pos`, `end`, `svtype`, `svlen`, `precise` (true iff
`stdev_pos + stdev_len` below threshold), `support`, `qual`, `filter`, `qc`, **five-point
coverage** (`upstream / start / center / end / downstream`) for transition detection, dict
of per-sample `genotypes`, `bnd_info`. Multi-sample merge via `SVGroup` (227–481).

**Genotyping (`genotyping.py:62–184`).** Binomial likelihood model with three genotypes,
`p = ε` for hom-ref, `0.5` for het, `1−ε` for hom-alt (ε ≈ 0.01 default). `GQ =
10 log10(L_secondbest / L_best)`. `GQ_z = 10 log10(L_homref / L_best)` for filtering
(`genotype_min_z_score`). **Per-type corrections**: long INS support rescaled by
`long_ins_rescale_mult × svlen / long_ins_length`; DUP coverage adjusted by `+0.75 × supp`
(extra-copy correction); INV reference coverage = upstream+downstream only.

**Consensus (`consensus.py:123–277`).** k-mer-based (k=15) MSA over per-read INS sequences,
seeded by the best Lead (closest to consensus pos+len). Majority vote at each anchor
position. **Does not transfer to us** — we have no per-read sequence.

**Filtering (`postprocessing.py:200–600`).** `min_support = base + mult × coverage_regional`
(auto) or fixed. `STDEV_POS > qc_stdev_abs_max` (20 bp) or `> 2.0 × svlen` ⇒ filter.
`STDEV_LEN / svlen > 1.0` ⇒ filter (non-BND). Coverage-transition tests for DEL/DUP. INS:
flag if upstream OR downstream coverage `< qc_coverage` (default 4).

**Mosaic mode (`postprocessing.py:467–589`).** Allow VAF in `[0.05, 0.5]`. Lower
`mosaic_min_reads = 2` but tighten precision: `(stdev_len/svlen > 0.1 OR stdev_pos > 5 bp)
AND svlen ≤ max_svlen_mosaic ⇒ filter`. Disable GT Z-score filter. Edge-proximity filter:
SV near read 5'/3' end < `dev_min_close_edge_dist` (10 bp) ⇒ filter.

## 4. Other tools — brief

**SAVANA** (Cortes-Ciriano lab, *Nature Methods* 2025; `savana/breakpoints.py` +
`classify.py`). Somatic-only LR caller. Pipeline: split-read breakpoint reconstruction →
proximity cluster → consensus → **ML classifier** (random forest, trained on Truvari-
labelled training set) producing per-call somatic / germline / artefact probability. CNA
pipeline on top: SV breakpoints seed circular binary segmentation; B-allele frequency at
LOH regions estimates tumour purity → absolute CN fit. **Useful idea:** ML classifier
trained on a labelled callset is the cleanest answer for "is this call real?" when no
clean rule exists. Same model can be retrained on our truth-labelled bubble calls.

**cuteSV** (`cuteSV_resolveINDEL.py`, `cuteSV_resolveINV.py`, `cuteSV_resolveDUP.py`).
Per-SV-type sliding-window scanning cluster (`max_cluster_bias` ≈ 200 bp for INS/DEL,
500 bp for DUP). Greedy single-pass: extend current cluster if new signature's pos within
bias, else flush. Simpler than Sniffles2's adaptive merge; harder to tune across coverage
regimes. **Useful idea:** per-SV-type dedicated resolver functions keep the per-type
logic isolated and easy to reason about.

**NanomonSV / dysgu / SVIM** — all in the same family. None do anything Severus +
Sniffles2 don't already cover for our purposes.

## 5. What transfers to our bubble pipeline

Our input is **fundamentally different**: no alignments, no contigs, no per-read
position. We have

- per-signature multiplicity `(n, t)` (from `counts.bin`),
- within-read adjacency edges (from `edges.bin`),
- bubble = connected component of non-shared nodes after `min_mult` pruning,
- a signature → reference-position oracle by extracting LCP cores from the truth/germline
  assembly (rule 1 allows this for analysis).

| concept | transfers? | notes |
|---|---|---|
| Severus `Breakpoint` (pos + dir + connections) | ✅ | Signature mapped to ref position; "direction" = sig's strand-mapping in ref (forward vs RC); "connections" = the bubble edges incident on this signature. |
| Severus `DoubleBreak` (pair of breakpoints) | ✅ | A bubble's two "external" anchors (shared backbone nodes flanking the bubble) ARE a DoubleBreak. |
| Severus `node = genomic segment` | ✅ | Collapse colinear ref signatures into a segment between two anchors. |
| Severus CC clustering + indel heuristic | ✅ | Apply directly to the bubble's internal edges (which we already have). |
| Severus templated-insertion / chromothripsis detection by topology + orientation parity | ✅ | Bubble's signature-pair orientation in the ref oracle gives parity. Need to compute. |
| Severus Bayesian GT (DV/DR + binomial) | ✅ | DV = sum of t-count over bubble's tumor-only nodes; DR = spanning normal reads at the flanking anchors. We have both. |
| Severus somatic filter (tumor support vs normal-spanning) | ✅ | Already what `colour` does at the node level; lift to bubble level. |
| Severus VCF schema (`CLUSTERID`, `DETAILED_TYPE`, `hVAF`, paired BND for complex) | ✅ | Emit identical schema. Compatible with Truvari/svdb without adapter. |
| Severus VNTR snap + insertion-deletion merge | ✅ | Need a VNTR BED for our reference. Snap signature-derived breakpoint to nearest TR boundary; merge co-located INS+DEL bubbles inside a VNTR. |
| Sniffles2 binned + adaptive merge cluster | △ | The within-bubble clustering is already done by graph BFS. **Cross-bubble** clustering (merging two adjacent bubbles representing the same real variant) is the same problem and the algorithm transfers — but the distance metric must use bubble-anchor positions, not raw read positions. |
| Sniffles2 TR-aware relaxed merge | ✅ | Same formula `(svlen_A + svlen_B) * cluster_repeat_h` with bubble-derived svlen estimate. |
| Sniffles2 binomial genotyping (`p = ε, 0.5, 1-ε`) | ✅ | Use `(k, n) = (variant-supporting bubble reads, total spanning reads at flank)`. Identical math. |
| Sniffles2 long-INS support rescale, DUP coverage adjust | ✅ | Same per-type corrections, applied to bubble-aggregated counts. |
| Sniffles2 STDEV_POS / STDEV_LEN filters | △ | We don't have per-read positions inside the bubble. Use per-signature reference-position spread (across the bubble's nodes) as a proxy. |
| Sniffles2 coverage-transition filter (DEL/DUP) | ✅ | Compute coverage from `n_count + t_count` at the flanking-anchor signatures. |
| Sniffles2 mosaic mode (VAF window, edge filter) | ✅ | Identical thresholds on `t_count / (n_count + t_count)`. |
| **Sniffles2 / Severus consensus assembly of allele sequence** | ✅ but bespoke | We **do not have read sequences**. Two options: (a) symbolic — emit a structural-only call with `<DEL>`/`<INS>` placeholders, sequence in INFO is bubble's signature path; (b) reconstructive — store a representative read substring per bubble node during the `lcp` stage and stitch them. (b) is the right answer and was already on the TODO from `[[project_svp_read_pipeline]]`. |
| SAVANA ML classifier | △ | Worth keeping in reserve. Requires labelled training (which we have via the truth VCF). Borrow if rule-based filters underperform. |
| cuteSV per-SV-type resolver functions | ✅ | Architectural style, not algorithm. Adopt: one resolver per SV type, each consuming bubble candidates filtered to that type. |

## 6. Proposed bubble→variant architecture

Three new stages on top of the existing `lcp / graph / bubbles`:

```
svp posmap   <ref.fa[.gz]>      <dir>     [--level 5 --flank 10]
                 build the signature → ref-position(s) oracle. Records:
                   sig_u64 -> { (chrom_id, pos_start, pos_end, strand) ... }
                 Per memory's lcp_read_graph::build_ref_posmap (line ~102)
                 generalised to multi-record / WG. Output: dir/posmap.bin.

svp resolve  <dir>                          [--min-supp 3 --vntr <bed>]
                 1. Load posmap.bin + graph.bin + bubbles (from bubbles.tsv +
                    re-derived node membership; per docs/PERFORMANCE.md the
                    BFS is cheap so we redo it instead of storing sigs_per_bubble).
                 2. For each bubble:
                    a) Map every node's signature -> {ref positions}. The bubble's
                       "consensus locus" = densest cluster of these positions.
                       Borrow Severus segment construction.
                    b) Anchor signatures = the shared-backbone (colour=3) nodes
                       adjacent to the bubble via the saved edges. Their positions
                       are the flanking breakpoints (left, right). Borrow
                       Severus DoubleBreak.
                    c) Classify type from the (left, right, internal) triple:
                         INS  -> left ≈ right (no ref gap), internal has cores not in ref
                         DEL  -> left < right by >50 bp, internal cores ⊂ ref-only colour
                         DUP  -> internal cores match ref segment whose copies in tumor > normal
                         INV  -> internal signatures map to ref with reverse strand
                                 vs flanking strand
                         BND  -> left and right map to different contigs OR
                                 inconsistent strand orientation
                         CPLX -> >1 of the above co-located, or ambiguous topology
                       (Per-type resolver functions, cuteSV style.)
                    d) VNTR snap (Severus / Sniffles2 rule): if either breakpoint
                       within 100 bp of a TR boundary in <bed>, snap to the boundary;
                       if INS+DEL bubbles co-locate inside a VNTR with compatible
                       lengths, merge to a single templated-insertion call.
                    e) Cross-bubble cluster (Sniffles2 adaptive merge): two bubbles
                       with overlapping flanking ranges and compatible svlen merge.
                 3. Genotype each call (Severus / Sniffles2 binomial model):
                       DV_tumor   = Σ t_count over the bubble's tumor-only nodes
                       DR_tumor   = max(0, spanning_tumor_reads_at_flank - DV_tumor)
                       VAF_tumor  = DV_tumor / (DV_tumor + DR_tumor)
                       GT, GQ     = arg max binomial likelihood under {0/0, 0/1, 1/1}
                 4. Filter (Sniffles2 stack adapted):
                       support_min   = base + mult * regional_cov
                       stdev_pos_max = 20 bp on signature ref-position spread
                       coverage_dip  = enforce on DEL flanks; coverage_jump on DUP
                       inside_vntr   = flag, do not auto-fail
                       single_strand = filter unless support >= mosaic threshold
                 5. Emit variants[] -> dir/variants.bin (checkpoint).

svp vcf      <dir>                          [--ref-fai <fa.fai>]
                 Read variants.bin, emit Severus-shape VCF:
                   simple SVs        -> single record  SVTYPE / END / SVLEN / STRANDS
                   complex events    -> paired BND records linked by CLUSTERID
                   INFO              -> PRECISE | IMPRECISE, MAPQ, DETAILED_TYPE,
                                        INSIDE_VNTR, SUPP_READS, REF_READS
                   FORMAT            -> GT:VAF:hVAF:DR:DV
                 Output: dir/somatic.vcf (passing somatic) + dir/all.vcf (all).
```

**Why this shape:**

- Posmap = the missing piece. The old caller's `LinearReference::align` was a buggy
  in-line attempt at the same thing; pulled out into its own stage and its own checkpoint
  makes it reusable for analysis tools (it's already what `tools/bubble_recall` would need
  to do per-bubble) and avoids redoing the LCP-on-reference cost.
- The bubble → DoubleBreak mapping (anchor signatures = colour-3 neighbours of the
  bubble's connected component) is **the bridge between coordinate-free bubbles and
  Severus's coordinate-aware schema**. Once we have that bridge, all Severus/Sniffles2
  logic above the breakpoint primitive applies.
- The 6 classifier steps (a–e + GT) plus the 5 filter steps are concrete enough to code
  without further design. None of them require alignment or assembly — only the posmap.
- VCF output is Severus-schema so Truvari / svdb / bcftools-merge work out of the box, and
  benchmarking against the truth VCF (18 002 variants in `simulation/truth/wgs/out.vcf`)
  is the existing recall pipeline with `truvari bench`.

**What's still missing for *exact* allele sequence (REF/ALT bases):**

This is the one piece neither Severus nor Sniffles2 can give us without per-read sequence.
The fix is in the `lcp` stage: store a representative read substring per signature
(roughly `[start-flank, end+flank]` of the first observed read carrying that signature).
~24 B header + ~40 B per representative seq = ~20 GB extra at WG scale, but only on the
output side — workers cache one substring per **new** signature, so per-emit cost is one
hashmap lookup. This was already a known TODO from `[[project_svp_read_pipeline]]`
("trim/calling + recall still TODO") and is the prerequisite for emitting precise
INS/SNP sequences. Without it, INS goes out as `<INS>` symbolic with `SVLEN`.

## 7. Concrete next steps (in order)

1. **`svp posmap`** stage + `dir/posmap.bin` format. Generalise `lcp_read_graph::
   build_ref_posmap` from chr19 to multi-record / WG. ~150 LOC.
2. **`svp resolve`** stage with steps a–d above (skip GT and filter first pass) emitting
   a minimal `variants.bin` containing `(chrom, pos, end, svtype, svlen, dv, dr,
   flanking_strand)`. ~400 LOC.
3. **`svp vcf`** emitter producing Severus-schema VCF. Validate that `truvari bench
   --base truth/wgs/out.vcf --comp dir/somatic.vcf` runs without parse errors. ~250 LOC.
4. **Filter + GT** (steps e and GT in §6): wire the Sniffles2 binomial model + the
   coverage/STDEV filters. Measure recall/precision delta. ~200 LOC.
5. **Allele-sequence emission**: add representative-read-substring caching to `process_read`,
   thread it through `posmap` + `resolve`. INS/SNP precise sequences. ~150 LOC in graph.cpp
   + ~100 LOC in resolve. (Per `[[project_svp_read_pipeline]]` next-steps.)
6. **VNTR snap + merge** (Severus rule). Requires a VNTR BED for the reference; can use
   UCSC's simpleRepeat track. ~100 LOC.
7. **Cross-bubble cluster** (Sniffles2 adaptive merge, TR-aware). ~150 LOC.
8. **Complex-SV typing** (templated insertion / chromothripsis proxy). ~150 LOC. Optional
   in v1; required to claim "complex SV support".

Steps 1–4 give a callable VCF; 5 makes it precise; 6–8 raise recall in repeats and
complex regions. Estimated ~1 500 LOC total across `src/posmap.{cpp,h}`,
`src/resolve.{cpp,h}`, `src/vcf.{cpp,h}`, plus `process_read` modifications. Within the
700-LOC src budget if old `match.cpp`-era code paths stay deleted.

## 8. References

- **Severus**: Keskus, A. G., *et al.* "Detection of complex structural variants from
  long reads using Severus." *Nature Biotechnology* (2024). Source:
  https://github.com/KolmogorovLab/Severus
- **Sniffles2**: Smolka, M., *et al.* "Detection of mosaic and population-level structural
  variants with Sniffles2." *Nature Biotechnology* 42, 1571–1580 (2024). Source:
  https://github.com/fritzsedlazeck/Sniffles
- **SAVANA**: Hartmann, H., *et al.* "SAVANA: reliable analysis of somatic structural
  variants and copy number aberrations in clinical samples using long-read sequencing."
  *Nature Methods* (2025). Source: https://github.com/cortes-ciriano-lab/savana
- **cuteSV**: Jiang, T., *et al.* "Long-read-based human genomic structural variation
  detection with cuteSV." *Genome Biology* 21, 189 (2020). Source:
  https://github.com/tjiangHIT/cuteSV

Source copies kept under `outputs/svp/research/` for offline re-reading.
