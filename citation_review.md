# Citation review of the janus main text (2026-10-23)

Inventory of every citing sentence: paper/citations.md (regenerated after the fixes). Four Opus reviewers checked each (sentence, work) pair against the cited paper, one per section; metadata of every cited entry was checked against Crossref by DOI (all resolve; Minda and scikit-learn have no DOI, the HG008 benchmark is an FTP release).

## Verdicts per cited pair (before the fixes)

### introduction

- S13 / chaisson: **OK**. Chaisson et al. 2019 abstract defines indels as <50 bp and SVs as >=50 bp (27,622 SVs >=50 bp per genome). The 50 bp threshold is the field convention and this paper uses it.
- S14 / cosenza: **OK**. Review of somatic SV roles and mechanisms; abstract: DNA repair and replication processes plus mitotic errors produce a spectrum from simple to highly complex rearrangements.
- S14 / chromothripsis: **OK**. Stephens et al. 2011 describes chromothripsis, an example of chromosome-scale somatic rearrangement in cancer.
- S15 / cosenza: **OK**. Abstract states verbatim that pan-cancer studies show SVs 'are the predominant class of driver mutation in most cancer types' (55% of PCAWG drivers are SVs).
- S15 / hanahan2011: **OK**. Establishes genome instability and mutation as an enabling characteristic that fuels acquisition of the hallmarks (second clause of the sentence).
- S15 / hanahan2022: **OK**. Updated hallmarks review; keeps genome instability as an enabling characteristic. It supports the second clause only, as hanahan2011 does.
- S16 / mcgranahan: **OK**. Review of clonal heterogeneity and branched tumor evolution, in which inherited and ongoing alterations drive clonal diversification.
- S16 / ith: **OK**. Dentro et al. 2021 show that subclonal diversification, including subclonal SVs, is pervasive across 2,658 cancers.
- S17 / schram2017: **OK**. Reviews the clinical relevance of gene fusions (SV products) for diagnosis and targeted therapy. Supports clinical importance.
- S17 / severus: **OK**. Severus abstract: somatic SV detection in heterogeneous, complex tumor genomes is not handled by most current methods, and it reports clinically relevant cryptic rearrangements. Supports both clauses.
- S17 / nanomonsv: **OK**. Tumor-normal long-read somatic SV method; documents the difficulty of complex somatic SVs. Acceptable as support for 'technically demanding'.
- S18 / manta: **OK**. Manta uses paired-read and split-read evidence with local assembly in short-read data.
- S18 / svaba: **OK**. SvABA detects SVs by genome-wide local assembly of short reads.
- S18 / gripss: **OK**. Shale et al. 2022 covers GRIDSS2/GRIPSS (assembly-based breakend calling) integrated with PURPLE copy number, i.e. depth-based evidence. Acceptable. GRIDSS itself (Cameron et al. 2017, Genome Res, doi:10.1101/gr.222109.117) would be the more direct source for assembly-based breakend calling, but this pairing is not wrong.
- S19 / chaisson: **OK**. Chaisson et al. show that short-read callsets miss most SVs, mainly insertions in repetitive or GC-rich regions, which the multi-platform approach recovers. logsdon (already in refs.bib) would also support this claim.
- S20 / chaisson: **OK**. Abstract: 'a three to sevenfold increase in SV detection compared to most standard high-throughput sequencing studies, including those from the 1000 Genomes Project'. The number is correct. Nuance: the gain comes from an integrated multi-platform callset (long reads, Strand-seq, optical maps), compared with standard short-read studies rather than with short reads of the same genomes. The sentence's word 'Integrating' fits this; it should not be read as long reads alone.
- S21 / sedlazeck: **OK**. NGMLR/Sniffles paper shows that long reads spanning repeats and SVs improve SV detection on PacBio and ONT data.
- S21 / logsdon: **OK**. Review of long-read sequencing (PacBio, ONT) and how it resolves repeats and SVs.
- S21 / Ahsan_2023: **OK**. Survey of long-read SV detection algorithms; covers the advantages of long reads for SVs.
- S22 / wenger: **WRONG**. The number is misstated. Wenger et al. 2019 report HiFi reads that are 'highly accurate (99.8%)' with 13.5 kb average length, which is a mean accuracy for that dataset. It is not a lower bound ('>='). HiFi reads are defined as >=Q20 (>=99%), so '>=99.8%' is not what the paper reports or what the platform guarantees. The qualitative claim (higher accuracy than ONT, shorter reads) is supported. Fix: Change to '...PacBio HiFi reads achieve ~99.8% average base accuracy [wenger]...', or '>=99% (Q20) per read, ~99.8% on average [wenger]'.
- S22 / logsdon: **OK**. Review compares PacBio HiFi and ONT: HiFi is more accurate, ONT gives longer reads.
- S22 / Delahaye_2021: **OK**. Characterizes ONT error rates (~6-8%) and systematic biases, including homopolymers and short repeats. Supports the claim that ONT is less suited for exact matching.
- S23 / pcawg: **OK**. PCAWG SV study (Li et al. 2020) analyzes bulk tumor whole genomes. It is an example of the norm, not a survey proving 'most studies', but that is acceptable for this statement.
- S24 / ith: **OK**. Dentro et al. find subclonal mutations in ~95% of tumors; many variants are present only in subpopulations.
- S25 / severus: **OK**. Severus addresses tumor heterogeneity and purity and evaluates low-VAF somatic SVs with matched normals. Acceptable.
- S25 / nanomonsv: **WEAK**. Nanomonsv presents a tumor-control somatic SV caller evaluated on 3 cancer cell lines. Its abstract does not discuss dilution by normal-cell contamination or VAF thresholds. Fix: Replace with 'ith' (Dentro et al. 2021, which models purity and subclonal VAFs explicitly) or 'savana' (SAVANA estimates tumor purity to handle dilution). Alternatively, keep severus alone.
- S26 / navin: **OK**. Navin et al. 2011 introduced single-cell sequencing to resolve tumor clones (first clause).
- S26 / gawad: **OK**. Review of single-cell genome sequencing that explicitly discusses amplification bias, coverage non-uniformity and allelic dropout (second clause).
- S28 / mcgranahan: **OK**. Review discusses multi-region sampling to capture heterogeneity.
- S28 / gerlinger: **OK**. The canonical multi-region sequencing study.
- S29 / severus: **OK**. Severus abstract: 'most current long-read SV detection methods are not developed for the analysis of tumor genomes characterized by complex rearrangements and heterogeneity'. Supports the claim that germline tools fail on low-VAF tumor variants. The specific 'half or all of the reads' mechanism is general knowledge. sniffles2 also notes that germline mode expects 50%/100% allele fractions.
- S30 / sniffles2: **OK**. Sniffles2 has a dedicated mosaic mode for low-frequency SVs and is not a tumor-normal somatic caller.
- S31 / severus: **OK**. Tumor-normal somatic SV caller.
- S31 / savana: **OK**. Somatic SV caller with matched germline control. It also runs without a control, but that does not contradict the sentence.
- S31 / nanomonsv: **OK**. Tumor/control paired somatic SV caller.
- S32 / Qin_2025: **OK**. Qin, Heinz and Li (2025) state that existing long-read SV callers only look at alignment against a single reference genome. Supports the sentence; the cited callers' own papers also confirm it.
- S33 / Ballouz_2019: **OK**. Opinion paper on the pitfalls of the GRCh38 reference for variant calling (reference bias). It supports the 'reference bias causes alignment and variant-calling errors' part, though not the somatic part. The two co-citations cover the somatic part, so the set is fine.
- S33 / Qin_2025: **OK**. Abstract: existing callers 'are susceptible to systematic false discovery caused by germline differences between the individual genome and the reference genome', i.e. germline variants called as somatic or mosaic.
- S33 / refbias: **OK**. Zhang et al. 2025 show that a donor-specific assembly as reference gives 'reduced germline artifacts' in somatic SV calling.
- S34 / hifiasm: **OK**. Proper source of the assembler used to build personal (donor-specific) assemblies. It does not itself discuss reference bias, but as a tool citation for 'personal assemblies' it is fine.
- S34 / Liao_2023: **OK**. Draft human pangenome reference: reduces small-variant errors and increases SV detection versus GRCh38, i.e. mitigates reference bias.
- S34 / Garrison_2018: **OK**. vg shows that mapping to a variation graph reduces reference bias, and it requires full read alignment to the graph.
- S34 / Alser_2020: **OK**. Primer showing that read mapping/alignment is a major computational bottleneck in genome analysis. Supports 'heavy computation'.
- S34 / Qin_2025: **OK**. Uses pangenome alignment plus a de novo personal assembly to reduce false somatic calls, and its pipeline requires full read alignments.
- S34 / refbias: **OK**. Donor-specific diploid assembly as reference, which requires full read mapping to the DSA.
- S35 / Pevzner_2001: **WEAK**. Establishes de Bruijn (Eulerian path) graphs for fragment assembly of a single genome. It does not show that k-mer/dBG approaches 'compare sequences' or 'capture sample differences directly'. Fix: Keep it as the dBG origin, but add cortex (Iqbal et al. 2012), which is the source showing dBGs capturing differences between samples. novobreak also fits.
- S35 / debruijn: **WEAK**. Compeau et al. 2011 is a tutorial on dBGs for genome assembly. It does not address comparing samples. Fix: Same fix as Pevzner_2001: add cortex to this sentence for the 'capture sample differences' claim.
- S36 / novobreak: **OK**. novoBreak identifies tumor-specific k-mers and locally assembles the reads containing them.
- S36 / cortex: **OK**. Cortex introduced colored de Bruijn graphs representing multiple samples.
- S37 / roberts: **OK**. Proper source of minimizers (named concept). The drawbacks are elaborated in S38.
- S37 / rustmdbg: **WEAK**. rust-mdbg is neither the source of minimizers nor of syncmers, and it does not establish their drawbacks. It is an application (minimizer-space dBG assembly). Fix: Remove from this sentence, or move it to S35/S36 as an example of a minimizer-space graph. roberts and syncmers already cover the named concepts.
- S37 / syncmers: **OK**. Proper source of syncmers. It also discusses minimizer shortcomings (conservation under mutation).
- S38 / roberts: **OK**. Window minimizers sample one k-mer per window, so selected markers leave gaps (up to w).
- S38 / syncmers: **OK**. Edgar 2021 analyzes k-mer conservation under mutations: a substitution destroys the k overlapping k-mers. It also covers the sparseness/gaps of submer-based selection.
- S38 / lcpan: **OK**. LCPan contrasts LCP with minimizers/syncmers ('minimizers are only guaranteed to have a maximal gap of size w between them') and stresses that every character is covered by a core.
- S39 / sv94a: **OK**. Original Sahinalp-Vishkin locally consistent parsing, which partitions the whole string into blocks.
- S39 / sv94b: **OK**. Original LCP/symmetry-breaking paper (STOC 1994).
- S40 / lcpan: **OK**. LCPan: 'identical cores share the same labels' (labeling and partition consistency).
- S40 / gencore: **OK**. GenCore: partition consistency and labeling consistency, i.e. 'cores that represent the same underlying sequence are assigned the same label'.
- S41 / lcpan: **WRONG**. The number '2.3 cores on average' affected by an isolated error does not appear in LCPan. I checked arXiv v1-v4 and the Genome Biology full text (PMC13270696). The ~2.3 figure in LCPan is c ~ 2.34, the per-level growth factor of core length and spacing (the core count shrinks by ~0.43 per level). It is not the number of cores disrupted by a mutation. The only bound I found is in GenCore (Sec. 2.2): 'a base substitution can change at most three consecutive cores at level 1 ... at most three level-2 cores are changed ... the number of cores changed at each level remains constant.' The second clause (iterative parsing yields longer cores across levels) is supported by LCPan. Fix: Rewrite as 'A single base substitution changes at most three cores at each level [gencore], and iterative parsing produces cores ~2.3 times longer at each level that capture structural context [lcpan].' If 2.3 is the authors' own empirical measurement, report it as such rather than citing lcpan.
- S47 / severus: **WEAK**. The Severus paper benchmarks on a multi-technology tumor/normal cell-line panel from the same consortium, so it is relevant for the cell lines and the multi-caller benchmarking practice. However, the CASTLE dataset (Cancer Standards Long-read Evaluation) was introduced and named in the DeepSomatic paper. The '18 call sets' consensus is the authors' own construction. Fix: Add: Park J, Cook DE, Chang P-C, Kolesnikov A, ..., Keskus AG, ... (2026). Accurate somatic small variant discovery for multiple sequencing technologies with DeepSomatic. Nature Biotechnology. doi:10.1038/s41587-025-02839-x (preprint doi:10.1101/2024.08.16.608331). Keep severus. If the consensus is built with minda, cite the existing key 'minda' where the consensus is described.

Uncited statements flagged:

- S27: Bulk sequencing therefore remains the standard approach. Proposed: pcawg (existing key), optionally with ith. A statement about current field practice. It reads as a conclusion from S26, but a reader may expect support. Low priority: adding pcawg/ith or leaving it uncited as a transition are both defensible.

### methods

- S52 / lcpan: **OK**. LCPan (Ashyralyyev et al. 2026) describes iterative LCP that produces variable-length cores at multiple levels and introduces the Lcptools implementation janus uses. It supports the LCP part of the sentence. The homopolymer-compression clause is not covered by either cited work (see uncited entry for S52/S68).
- S52 / sv94a: **OK**. Sahinalp & Vishkin 1994 (LNCS 778) is one of the two works that LCPan itself cites as introducing LCP ('Locally Consistent Parsing (LCP) [21,22]', where [21] is this paper and [22] is the STOC'94 suffix-tree paper). It is an appropriate original-method citation.
- S55 / cortex: **WEAK**. Cortex (Iqbal et al. 2012) is the right source for the colored de Bruijn graph and does clean errors from the graph (tip and low-coverage supernode removal, in its Supplementary Methods). But it does not pop bubbles: bubbles are its variant-calling signal, and it classifies them statistically as errors, repeats or variants. So it does not establish 'standard tip clipping and bubble popping'. In S55 it is used for the graph-cleaning operations and is not cited where the colored dBG is introduced (S54). Fix: For tip clipping and bubble popping, cite Velvet: Zerbino DR, Birney E. 'Velvet: algorithms for de novo short read assembly using de Bruijn graphs.' Genome Research 18(5):821-829, 2008. doi:10.1101/gr.074492.107 (tip removal and the Tour Bus bubble-merging algorithm). It is not in refs.bib. Move cortex to S54 (colored de Bruijn-like graph), or keep it here alongside Velvet.
- S57 / Li_2018: **OK**. minimap2 is named, and Li 2018 is its canonical citation. pipeline.sh confirms minimap2 (-cx map-hifi) is used for the snippet-to-window alignment.
- S68 / lcpan: **OK**. LCPan introduces Lcptools ('We present the first iterative implementation of LCP with Lcptools', github.com/BilkentCompGen/lcptools) and analyses cores per level. It is the proper source for LCPtools and leveled LCP cores. Homopolymer compression in the same sentence is uncited (see uncited entry).
- S105 / Li_2018: **OK**. minimap2 is named. Li 2018 is the correct source.
- S110 / severus: **OK**. Severus Methods: junctions are pairs of breakpoints, and 'Each breakpoint is directional, either head or tail (corresponding to 3' and 5' directions)'. This matches 'two oriented breakpoint coordinates'. Bib details check out: Keskus et al., Nat Biotechnol 44(2):247-257, 2026 issue (online April 2025), doi:10.1038/s41587-025-02618-8.
- S126 / sklearn: **OK**. janus uses sklearn.ensemble.HistGradientBoostingClassifier (src/janus/classify.py), so scikit-learn is the proper software citation. But the sentence never names scikit-learn, and the 2011 paper is a software source, not the source of the gradient-boosting method.

Uncited statements flagged:

- S52 / S68: Homopolymer compression of the reads (S52 'After homopolymer compression'; S68 'janus homopolymer-compresses each read'). Both sentences carry citations, but only for LCP/LCPtools. Proposed: Li_2018 (existing key). Homopolymer compression is an external technique. The minimap2 paper, Section 2.1.5 'Indexing with homopolymer compressed k-mers', defines HPC(s) as 'contracting homopolymers in s to a single base' and credits the idea to SmartDenovo, which has no paper. LCPan does not cover it.
- S54: janus constructs a colored de Bruijn-like graph with per-sample read counts on each node. Proposed: cortex (existing key); optionally debruijn (Compeau, Pevzner & Tesler 2011) for de Bruijn graphs in general. The colored de Bruijn graph is an external concept introduced by Iqbal et al. 2012 (Cortex). This is where the concept is first named in the Methods, but the citation currently sits on the next sentence.
- S59: A gradient-boosted classifier filters candidate records. Proposed: sklearn (existing key); optionally Friedman JH 2001, Annals of Statistics 29(5):1189-1232, doi:10.1214/aos/1013203451. The external ML method is named here, but it is cited only later (S126).
- S76: Thresholds were selected by a parameter sweep 'over the CASTLE runs'. Proposed: severus (existing key). This is the first mention of the CASTLE panel (an external dataset) in the Methods. The Severus paper is its source: 'cancer standards long-read evaluation', six tumor-normal cell-line pairs. The Results cite severus for CASTLE, but the Methods do not. S85 ('across the CASTLE pairs') is covered by the same fix.
- S109: SV signatures are taken from internal indels of at least 50 bp. This relies on the standard 50 bp SV definition. Proposed: chaisson (existing key; already used in the Introduction for 'SVs span at least 50 bp'). Standard definition. A citation is optional because the Introduction already cites chaisson for it, but the task lists the 50 bp threshold as citation-worthy.
- S111: Signatures of identical type and orientation whose breakpoints match within 50 bp are merged into one record. Proposed: severus (existing key), only if the 50 bp tolerance was adopted from Severus. Low priority. Severus clusters junctions with BP_DISTANCE = 50 bp by default, and the previous sentence says 'Following Severus'. If the authors took this convention from Severus, it should be credited. If it is their own choice, no citation is needed.
- S137: janus builds a neighbor-joining tree of the tumor samples. Proposed: Saitou N, Nei M. 'The neighbor-joining method: a new method for reconstructing phylogenetic trees.' Molecular Biology and Evolution 4(4):406-425, 1987. doi:10.1093/oxfordjournals.molbev.a040454 (not in refs.bib). Neighbor joining is a named external algorithm. src/janus/joint.py implements NJ itself (Q-criterion), so the R package ape (key 'ape', used in the Supplementary evaluation) is not the right source here; the original method paper is.

### results

- S147 / minda: **OK**. Minda is the KolmogorovLab tool for merging/benchmarking SV call sets (minda ensemble / minda truthset), as used here. The GitHub README asks users to cite the Severus paper (Keskus et al., Nat Biotechnol 2025), which introduces Minda in its Methods ('Minda algorithm').
- S151 / severus: **OK**. Published Severus paper (Keskus et al., Nat Biotechnol, doi 10.1038/s41587-025-02618-8) introduces the 'cancer standards long-read evaluation (CASTLE) panel': HCC1954, H2009, HCC1937, H1437, Hs578T and HCC1395 with Illumina, ONT R10 and PacBio HiFi; the five pairs used here are a subset. Data at SRA PRJNA1086849.
- S152 / minda: **OK**. Minda is the tool that built the ensemble (Severus paper: 'we created an ensemble of confident calls ... using Minda').
- S152 / severus: **OK**. Severus paper defines the CASTLE confident set: 'supported by at least two (of three) technologies and at least 4 (of 11) callers', built with Minda from Severus, SAVANA, nanomonsv, Sniffles2, SvABA, GRIPSS and Manta across HiFi, ONT and Illumina. Matches the sentence exactly.
- S152 / nanomonsv: **OK**. Shiraishi et al. 2023 NAR is the nanomonsv paper; nanomonsv is one of the callers of the 11 call sets (per Severus).
- S152 / manta: **OK**. Chen et al. 2016 Bioinformatics is the Manta paper; Manta is one of the three Illumina call sets in the consensus.
- S152 / svaba: **OK**. Wala et al. 2018 Genome Res is the SvABA paper; SvABA is one of the Illumina call sets in the consensus.
- S152 / gripss: **WRONG**. The key points to Shale et al. 2022 Cell Genomics ('Unscrambling cancer genomes...'), which is the LINX paper (GRIDSS-PURPLE-LINX interpretation). Its full text (Europe PMC XML, PMC9903802) never mentions GRIPSS; it is not the source of the GRIDSS/GRIPSS call set. The Severus paper itself cites GRIDSS (Cameron et al. 2017) for 'GRIDSS/GRIPSS'. Same misattribution in introduction.tex line 9 (\cite{manta,svaba,gripss}). Fix: Replace with the GRIDSS2 paper: Cameron DL, Baber J, Shale C, et al. 'GRIDSS2: comprehensive characterisation of somatic structural variation using single breakend variants and structural variant phasing', Genome Biology 22:202 (2021), doi 10.1186/s13059-021-02423-x (verified on Crossref; describes the GRIDSS2 somatic caller and its somatic filter). GRIPSS itself has no paper; add a software citation to https://github.com/hartwigmedical/hmftools/tree/master/gripss. Alternatively cite GRIDSS (Cameron et al., Genome Res 27:2050-2060, 2017, doi 10.1101/gr.222109.117), as the Severus paper does.
- S152 / sniffles2: **OK**. Smolka et al. 2024 Nat Biotechnol is the Sniffles2 paper; Sniffles2 is one of the consensus callers.
- S152 / savana: **OK**. Elrick et al. 2025 Nat Methods is the SAVANA paper; SAVANA is one of the consensus callers.
- S153 / dysgu: **OK**. Cleal & Baird 2022 NAR is the dysgu paper. Together with S152 (severus, savana, nanomonsv, sniffles2) and the introduction, each of the ten HiFi callers is cited with its own paper when first named.
- S153 / svisionpro: **OK**. Wang et al. 2025 Nat Biotechnol, 'De novo and somatic structural variant discovery with SVision-pro', is the SVision-pro paper. Together with S152 (severus, savana, nanomonsv, sniffles2) and the introduction, each of the ten HiFi callers is cited with its own paper when first named.
- S153 / cutesv: **OK**. Jiang et al. 2020 Genome Biol is the cuteSV paper. Together with S152 (severus, savana, nanomonsv, sniffles2) and the introduction, each of the ten HiFi callers is cited with its own paper when first named.
- S153 / debreak: **OK**. Chen et al. 2023 Nat Commun is the DeBreak paper. Together with S152 (severus, savana, nanomonsv, sniffles2) and the introduction, each of the ten HiFi callers is cited with its own paper when first named.
- S153 / kled: **OK**. Zhang et al. 2024 Brief Bioinform is the Kled paper. Together with S152 (severus, savana, nanomonsv, sniffles2) and the introduction, each of the ten HiFi callers is cited with its own paper when first named.
- S153 / svdss: **OK**. Denti et al. 2022 Nat Methods is the SVDSS paper. Together with S152 (severus, savana, nanomonsv, sniffles2) and the introduction, each of the ten HiFi callers is cited with its own paper when first named.
- S163 / colo829: **OK**. Espejo Valle-Inclan et al. 2022 Cell Genomics is the COLO829 multi-platform validated somatic SV truth set (melanoma). The truth set holds 68 validated SVs (Severus paper: 'confident set of 68 somatic SV calls'); 62 is the authors' subset of at least 50 bp (supplementary.tex line 222).
- S163 / hg008: **WEAK**. McDaniel et al. 2025 Sci Data establishes HG008 (GIAB pancreatic ductal adenocarcinoma cell line with matched normal) and its sequencing data, so 'HG008 (pancreatic tumor)' is supported. It does not contain the somatic SV draft benchmark: the paper says benchmarks will be developed later and SV calls are out of its scope. 'GIAB draft v0.5' (and the 172 SVs) come from the NIST release NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5-20260318 on the GIAB FTP. Fix: Keep hg008 for the sample/data and add a data citation for the benchmark: NIST Genome in a Bottle, 'HG008-T somatic structural variant and CNV draft benchmark V0.5' (2026-03-18), https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NIST_HG008-T_somatic-stvar-CNV_DraftBenchmark_V0.5-20260318/ (file GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz). No peer-reviewed publication of this draft benchmark was found.
- S198 / clair3: **OK**. Zheng et al. 2022 Nat Comput Sci is the Clair3 paper; the supplement confirms Clair3 v2.0.3 was used to call small variants of the normal for phasing.
- S198 / longphase: **OK**. Lin et al. 2022 Bioinformatics is the LongPhase paper; used for phasing and haplotagging as stated in the supplement.

Uncited statements flagged:

- S147: Competing callers ran on minimap2 alignments (first mention of minimap2 in Results; repeated in S186, S196, S199 and the Fig. runtime caption) Proposed: Li_2018 (Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094-3100, 2018, doi 10.1093/bioinformatics/bty191). External tool named and relied on; it is cited in Methods but not here, and supplementary.tex line 200 also names minimap2 v2.31 without a citation. One citation at the first Results mention suffices.
- S149: 95% bootstrap confidence intervals over 2,000 resamples Proposed: Efron B. Bootstrap methods: another look at the jackknife. Annals of Statistics 7(1):1-26, 1979, doi 10.1214/aos/1176344552 (verified on Crossref; not in refs.bib). Named statistical method; low priority, as the bootstrap is standard, but a reference is conventional.
- S174: Simulated data with one normal genome carrying real germline variants and 20 lineage-related tumor clones Proposed: hgsvc3 (germline variants of HG00733 from the HGSVC3 call set), Battistella_2025 (insilicoSV), msprime (tree); all already in refs.bib and cited in the supplement. 'Real germline variants' relies on an external dataset (HGSVC3), and the simulation relies on external simulators; the Results point only to a Section.
- S176: We simulated 20x PacBio HiFi reads per genome Proposed: badread (Wick R. Badread: simulation of error-prone long reads. JOSS 4(36):1316, 2019, doi 10.21105/joss.01316). The read simulator is an external tool; it is cited in the supplement but not in the Results.
- S187: A neighbor-joining tree ... Robinson-Foulds distance 0 Proposed: Saitou N, Nei M. The neighbor-joining method: a new method for reconstructing phylogenetic trees. Mol Biol Evol 4(4):406-425, 1987, doi 10.1093/oxfordjournals.molbev.a040454; Robinson DF, Foulds LR. Comparison of phylogenetic trees. Mathematical Biosciences 53:131-147, 1981, doi 10.1016/0025-5564(81)90043-2 (both verified on Crossref, neither in refs.bib); optionally ape (existing key, the implementation used). Two named external methods with standard primary references; neither is cited anywhere in the paper (only ape in the supplement).

### discussion

- S215 / slidedna: **OK**. Zhao et al. 2022 (Nature 601:85-91) introduce Slide-DNA-seq, which captures spatially resolved DNA from intact tissue sections and discovers distinct tumour clones and their copy-number alterations in situ, showing that clones occupy distinct spatial regions (mouse metastasis model and a primary human cancer). This directly supports the outlook that spatially resolved DNA sequencing can locate subclones in tissue. Caveat (not a misattribution): the paper resolves copy-number alterations at low coverage, not breakpoint-level SVs such as balanced translocations or inversions. The sentence says 'could', so the citation fits as written.
- S216 / schram2017: **OK**. Schram et al. 2017 (Nat Rev Clin Oncol 14:735-748) review gene fusions arising from structural rearrangements in solid tumours, how they are detected in the clinic, kinase inhibitors that target them (ALK, ROS1, RET, NTRK, FGFR, BRAF) and the newer-generation inhibitors developed against resistance. This supports the claim that somatic SV detection guides the choice and development of targeted drugs. It covers only fusions, which are the clinically best-established class of actionable SVs, so it is adequate for this sentence. The 'pangenome-based methods scale' condition is the authors' own outlook and needs no citation.

Uncited statements flagged:

- S214: Reads of other long-read technologies such as ONT could also enter the graph. Proposed: Existing keys: Deamer_2016 (Three decades of nanopore sequencing, 2016) and/or logsdon (Long-read human genome sequencing and its applications, 2020); Delahaye_2021 is already cited in the Introduction for ONT error profiles and would also fit if the higher error rate of ONT is mentioned.. Low priority. The sentence is an outlook for janus, but it names an outside technology (Oxford Nanopore). The Discussion does not cite it here, though the Introduction cites logsdon and Delahaye_2021 for long-read platforms. A reader may expect a pointer at first mention of ONT in this section. Optional, not required.

## Changes made

- introduction.tex
  - was: PacBio HiFi reads achieve $\ge$99.8\% base accuracy~\cite{wenger}
  - now: PacBio HiFi reads achieve 99.8\% average base accuracy~\cite{wenger}
- introduction.tex
  - was: variant allele frequency (VAF) thresholds~\cite{severus,nanomonsv}
  - now: variant allele frequency (VAF) thresholds~\cite{severus,ith}
- introduction.tex
  - was: which capture sample differences directly~\cite{Pevzner_2001,debruijn}
  - now: which capture sample differences directly~\cite{Pevzner_2001,debruijn,cortex}
- introduction.tex
  - was: have notable drawbacks~\cite{roberts,rustmdbg,syncmers}
  - now: have notable drawbacks~\cite{roberts,syncmers}
- introduction.tex
  - was: An isolated error affects only 2.3 cores on average~\cite{lcpan}, and iterative parsing produces longer cores that capture structural context across multiple levels.
  - now: A single base substitution changes at most three cores at each level~\cite{gencore}, and iterative parsing produces longer cores that capture structural context across multiple levels~\cite{lcpan}.
- introduction.tex
  - was: On the CASTLE tumor--normal cell line pairs~\cite{severus}
  - now: On the CASTLE tumor--normal cell line pairs~\cite{deepsomatic,severus}
- methods.tex
  - was: After homopolymer compression, \tool partitions reads
  - now: After homopolymer compression~\cite{Li_2018}, \tool partitions reads
- methods.tex
  - was: constructs a colored de~Bruijn-like graph from these cores
  - now: constructs a colored de~Bruijn-like graph~\cite{cortex} from these cores
- methods.tex
  - was: After standard tip clipping and bubble popping~\cite{cortex}, \tool identifies
  - now: After standard tip clipping and bubble popping~\cite{velvet}, \tool identifies
- methods.tex
  - was: using gradient-boosted decision trees~\cite{sklearn}.
  - now: using gradient-boosted decision trees~\cite{friedman,sklearn}.
- methods.tex
  - was: \tool builds a neighbor-joining tree of the tumor samples and re-decides
  - now: \tool builds a neighbor-joining tree~\cite{nj} of the tumor samples and re-decides
- results.tex
  - was: we ran competing PacBio HiFi-based callers using minimap2 alignments generated
  - now: we ran competing PacBio HiFi-based callers using minimap2~\cite{Li_2018} alignments generated
- results.tex
  - was: we computed 95\% bootstrap confidence intervals for the F1 score differences
  - now: we computed 95\% bootstrap confidence intervals~\cite{efron} for the F1 score differences
- results.tex
  - was: cell line pairs from the CASTLE panel~\cite{severus} (Table~\ref{suptab:castle}). Lacking
  - now: cell line pairs from the CASTLE panel~\cite{deepsomatic,severus} (Table~\ref{suptab:castle}). Lacking
- results.tex
  - was: HG008 (pancreatic tumor, GIAB draft v0.5, 172 SVs)~\cite{hg008}. Trained on all five
  - now: HG008 (pancreatic tumor, GIAB draft v0.5, 172 SVs)~\cite{hg008,hg008bench}. Trained on all five
- results.tex
  - was: we evaluated joint calling on simulated data: one normal genome with real germline variants and 20 lineage-related tumor clones that inherit ancestral somatic SVs (Section
  - now: we evaluated joint calling on simulated data~\cite{Battistella_2025,msprime,badread}: one normal genome with real germline variants~\cite{hgsvc3} and 20 lineage-related tumor clones that inherit ancestral somatic SVs (Section
- results.tex
  - was: and a neighbor-joining tree built on their pairwise SV differences recovered all 17 simulated splits (Robinson--Foulds distance 0; Fig.
  - now: and a neighbor-joining tree~\cite{nj} built on their pairwise SV differences recovered all 17 simulated splits (Robinson--Foulds distance~\cite{rf} 0; Fig.
- refs.bib
  - was: @article{gripss, author={Shale, Charles and Cameron, Daniel L. and Baber, Jonathan and Wong, Marie and Cowley, Mark J. and Papenfuss, Anthony T. and Cuppen, Edwin and Priestley, Peter}, title={Unscrambling cancer genomes via integrated analysis of structural variation and copy number}, journal={Cell
  - now: @article{gripss, author={Cameron, Daniel L. and Baber, Jonathan and Shale, Charles and Espejo Valle-Inclan, Jose and Besselink, Nicolle and van Hoeck, Arne and Janssen, Roel and Cuppen, Edwin and Priestley, Peter and Papenfuss, Anthony T.}, title={{GRIDSS2}: comprehensive characterisation of somatic
- refs.bib
  - was: title={GenCore: Genomic distance estimation using Locally Consistent Parsing}
  - now: title={{GenCore}: Genomic distance estimation using Locally Consistent Parsing}
- refs.bib
  - was: author={Wang, Songbo and Lin, Jiadong and Jia, Peng and Xu, Tun and others}
  - now: author={Wang, Songbo and others}
- refs.bib
  - was: author={McDaniel, Jennifer H. and Patel, Vaidehi and Olson, Nathan D. and He, Hua-Jun and others}
  - now: author={McDaniel, Jennifer H. and others}
- refs.bib
  - was: author={Logsdon, Glennis A. and Ebert, Peter and Audano, Peter A. and Loftus, Mark and others}
  - now: author={Logsdon, Glennis A. and others}
- refs.bib
  - was: svdss year={2022}
  - now: year={2023} (volume 20 is the 2023 volume)
- refs.bib
  - was: (end of file)
  - now: 7 new entries: deepsomatic, velvet, friedman, nj, rf, efron, hg008bench
- introduction.tex, results.tex: order of the paired CASTLE citations set to severus, deepsomatic so the numbers print in order.
- results.tex: the Robinson--Foulds citation placed after the number ("distance 0~\cite{rf}").

## Flags left as they are

- COLO829 "62 validated SVs": the supplement (Section S2.3) already states that 62 of the 68 validated SVs are at least 50 bp; the main text is unchanged.
- Methods "the CASTLE runs" (parameter sweep), the 50 bp signature threshold and the 50 bp merging tolerance: own choices, CASTLE is cited at its first Methods mention.
- Introduction "Bulk sequencing therefore remains the standard approach" and Discussion "other long-read technologies such as ONT": transitions, both topics cited in the Introduction.
- Discussion / slidedna: the cited work resolves clone copy-number changes in tissue, not SV breakpoints; the sentence says "could", left as is.
- refs.bib keeps the unused entries (rustmdbg, Keskus_2025 duplicate of severus, the placeholder entries small and big); nothing removed.
