# Citations of the janus main text, in reading order

One line per sentence that cites something. Works are named by first author, year and title.


## abstract


## introduction

[Sentence 13: "Structural variants (SVs) are genomic alterations that span at least 50\,bp." -> cites: Chaisson 2019, "Multi-platform discovery of haplotype-resolved structural variation in"]

[Sentence 14: "Beyond inherited germline variants, replication and repair defects generate somatic SVs in cancer cells, ranging from focal deletions to chromosome-scale rearrangements." -> cites: Cosenza 2022, "Structural Variation in Cancer: Role, Prevalence, and Mechanisms"; Stephens 2011, "Massive Genomic Rearrangement Acquired in a Single Catastrophic Event "]

[Sentence 15: "SVs represent the primary class of driver mutations across most tumor types, and the resulting genomic instability fuels the acquisition of cancer hallmarks." -> cites: Cosenza 2022, "Structural Variation in Cancer: Role, Prevalence, and Mechanisms"; Hanahan 2011, "Hallmarks of cancer: the next generation"; Hanahan 2022, "Hallmarks of cancer: new dimensions"]

[Sentence 16: "Because cell lineages inherit these alterations, ongoing structural variation drives clonal diversification and tumor progression." -> cites: McGranahan 2017, "Clonal heterogeneity and tumor evolution: past, present, and the futur"; Dentro 2021, "Characterizing genetic intra-tumor heterogeneity across 2,658 human ca"]

[Sentence 17: "Consequently, rapid and accurate somatic SV identification is critical for clinical oncology, yet remains technically demanding." -> cites: Schram 2017, "Fusions in solid tumours: diagnostic strategies, targeted therapy, and"; Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"; Shiraishi 2023, "Precise characterization of somatic complex structural variations from"]

[Sentence 18: "Early approaches relied on short reads to infer SVs indirectly from discordant read pairs, split reads, depth fluctuations, and local assembly." -> cites: Chen 2016, "Manta: rapid detection of structural variants and indels for germline "; Wala 2018, "SvABA: genome-wide detection of structural variants and indels by loca"; Cameron 2021, "GRIDSS2: comprehensive characterisation of somatic structural variatio"]

[Sentence 19: "However, short reads map poorly to repetitive regions, limiting detection accuracy." -> cites: Chaisson 2019, "Multi-platform discovery of haplotype-resolved structural variation in"]

[Sentence 20: "Integrating long reads resolved this limitation, uncovering 3 to 7 times as many SVs per genome as short-read datasets alone." -> cites: Chaisson 2019, "Multi-platform discovery of haplotype-resolved structural variation in"]

[Sentence 21: "Long-read sequencing from Oxford Nanopore Technologies (ONT) and PacBio HiFi overcomes mapping ambiguity because individual reads can span repetitive elements or entire SVs." -> cites: Sedlazeck 2018, "Accurate detection of complex structural variations using single-molec"; Logsdon 2020, "Long-read human genome sequencing and its applications"; Ahsan 2023, "A survey of algorithms for the detection of genomic structural variant"]

[Sentence 22: "Between the two platforms, PacBio HiFi reads achieve 99.8\% average base accuracy, making them better suited than ONT reads for exact sequence matching and breakpoint resolution despite their shorter lengths." -> cites: Wenger 2019, "Accurate circular consensus long-read sequencing improves variant dete"; Logsdon 2020, "Long-read human genome sequencing and its applications"; Delahaye 2021, "Sequencing DNA with nanopores: Troubles and biases"]

[Sentence 23: "Most studies analyze somatic SVs in bulk samples that pool DNA from thousands of cells." -> cites: Li 2020, "Patterns of somatic structural variation in human cancer genomes"]

[Sentence 24: "Because ongoing mutations create distinct subclones, many variants occur only in small cell subpopulations." -> cites: Dentro 2021, "Characterizing genetic intra-tumor heterogeneity across 2,658 human ca"]

[Sentence 25: "Normal-cell contamination dilutes these mutations that often fall below standard variant allele frequency (VAF) thresholds." -> cites: Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"; Dentro 2021, "Characterizing genetic intra-tumor heterogeneity across 2,658 human ca"]

[Sentence 26: "Single-cell methods directly isolate distinct clones, but they introduce coverage gaps and allelic dropout that impede variant discovery." -> cites: Navin 2011, "Tumour evolution inferred by single-cell sequencing"; Gawad 2016, "Single-cell genome sequencing: current state of the science"]

[Sentence 28: "To improve clonal resolution, researchers sequence multiple bulk samples that capture tumor diversity." -> cites: McGranahan 2017, "Clonal heterogeneity and tumor evolution: past, present, and the futur"; Gerlinger 2012, "Intratumor heterogeneity and branched evolution revealed by multiregio"]

[Sentence 29: "Standard germline tools fail at low VAFs because they expect variants that appear in half or all of the reads." -> cites: Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"]

[Sentence 30: "Even specialized tools like Sniffles2 include mosaic modes that target rare alleles, yet they cannot distinguish somatic from germline events without normal data." -> cites: Smolka 2024, "Detection of mosaic and population-level structural variants with Snif"]

[Sentence 31: "Somatic discovery therefore requires dedicated tumor--normal callers that evaluate matched control samples to filter non-somatic signal." -> cites: Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"; Elrick 2025, "SAVANA: reliable analysis of somatic structural variants and copy numb"; Shiraishi 2023, "Precise characterization of somatic complex structural variations from"]

[Sentence 32: "All of the aforementioned tumor--normal callers detect SVs from reads that align to a single linear reference genome." -> cites: Qin 2025, "Improving long-read somatic structural variant calling with pangenome "]

[Sentence 33: "Reference bias causes alignment errors that often lead callers to misidentify germline variants as somatic SVs." -> cites: Ballouz 2019, "Is it time to change the reference genome?"; Qin 2025, "Improving long-read somatic structural variant calling with pangenome "; Zhang 2025, "Diploid donor-specific assembly enhances somatic structural variant de"]

[Sentence 34: "Personal assemblies and pangenome references mitigate this bias, but they still require full read alignments that demand heavy computation." -> cites: Cheng 2021, "Haplotype-resolved de novo assembly using phased assembly graphs with "; Liao 2023, "A draft human pangenome reference"; Garrison 2018, "Variation graph toolkit improves read mapping by representing genetic "; Alser 2020, "Accelerating Genome Analysis: A Primer on an Ongoing Journey"; Qin 2025, "Improving long-read somatic structural variant calling with pangenome "; Zhang 2025, "Diploid donor-specific assembly enhances somatic structural variant de"]

[Sentence 35: "Graph-based approaches instead compare sequences using $k$-mers and de Bruijn graphs, which capture sample differences directly." -> cites: Pevzner 2001, "An Eulerian path approach to DNA fragment assembly"; Compeau 2011, "How to apply de Bruijn graphs to genome assembly"; Iqbal 2012, "De novo assembly and genotyping of variants using colored de Bruijn gr"]

[Sentence 36: "Tools like novoBreak assemble reads that carry tumor-specific $k$-mers, while colored de Bruijn graphs represent multiple samples in a single model." -> cites: Chong 2017, "novoBreak: local assembly for breakpoint detection in cancer genomes"; Iqbal 2012, "De novo assembly and genotyping of variants using colored de Bruijn gr"]

[Sentence 37: "However, standard $k$-mers and sparse markers, such as minimizers and syncmers, have notable drawbacks." -> cites: Roberts 2004, "Reducing storage requirements for biological sequence comparison"; Edgar 2021, "Syncmers are more sensitive than minimizers for selecting conserved \e"]

[Sentence 38: "A single mutation disrupts up to $k$ overlapping $k$-mers, while subsampled markers leave sequence gaps that miss critical variation." -> cites: Roberts 2004, "Reducing storage requirements for biological sequence comparison"; Edgar 2021, "Syncmers are more sensitive than minimizers for selecting conserved \e"; Ashyralyyev 2026, "LCPan: efficient variation graph construction using locally consistent"]

[Sentence 39: "Locally consistent parsing (LCP) avoids these issues by generating variable-length cores that guarantee full sequence coverage." -> cites: Sahinalp 1994, "On a parallel-algorithms method for string matching problems (overview"; Sahinalp 1994, "Symmetry breaking for suffix tree construction"]

[Sentence 40: "LCP parses shared regions into identical cores that match consistently across reads." -> cites: Ashyralyyev 2026, "LCPan: efficient variation graph construction using locally consistent"; Ashyralyyev 2026, "GenCore: Genomic distance estimation using Locally Consistent Parsing"]

[Sentence 41: "A single base substitution changes at most three cores at each level, and iterative parsing produces longer cores that capture structural context across multiple levels." -> cites: Ashyralyyev 2026, "GenCore: Genomic distance estimation using Locally Consistent Parsing"; Ashyralyyev 2026, "LCPan: efficient variation graph construction using locally consistent"]

[Sentence 47: "On the CASTLE tumor--normal cell line pairs, janus predictions showed the highest average F1 score against a consensus of 18 call sets." -> cites: Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"; Park 2026, "Accurate somatic small variant discovery for multiple sequencing techn"]


## methods

[Sentence 52: "After homopolymer compression, janus partitions reads into variable-length cores across multiple granularity levels via locally consistent parsing (LCP)." -> cites: Li 2018, "Minimap2: pairwise alignment for nucleotide sequences"; Ashyralyyev 2026, "LCPan: efficient variation graph construction using locally consistent"; Sahinalp 1994, "On a parallel-algorithms method for string matching problems (overview"]

[Sentence 54: "janus then constructs a colored de Bruijn-like graph from these cores, tracking sample-specific read counts on each node (Fig. REFa)." -> cites: Iqbal 2012, "De novo assembly and genotyping of variants using colored de Bruijn gr"]

[Sentence 55: "After standard tip clipping and bubble popping, janus identifies edge pairs traversed exclusively by tumor reads as candidate somatic SV sites, detecting shared sites across multiple tumors only once (Fig. REFb,c)." -> cites: Zerbino 2008, "Velvet: Algorithms for de novo short read assembly using de Bruijn gra"]

[Sentence 57: "It maps snippets with minimap2 to short reference windows seeded by shared cores, and the alignments resolve breakpoints and SV types." -> cites: Li 2018, "Minimap2: pairwise alignment for nucleotide sequences"]

[Sentence 68: "janus homopolymer-compresses each read and partitions both strands into level-5 LCP cores using LCPtools, indexing each core by a 64-bit sequence hash that serves as its graph node identifier." -> cites: Ashyralyyev 2026, "LCPan: efficient variation graph construction using locally consistent"]

[Sentence 105: "janus aligns snippets with minimap2 to short candidate reference windows rather than the whole genome (Section REF)." -> cites: Li 2018, "Minimap2: pairwise alignment for nucleotide sequences"]

[Sentence 110: "Following Severus, each signature comprises two oriented breakpoint coordinates on the reference (Section REF)." -> cites: Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"]

[Sentence 126: "In a post-processing step, janus classifies candidate records from the VCF using gradient-boosted decision trees." -> cites: Friedman 2001, "Greedy function approximation: A gradient boosting machine"; Pedregosa 2011, "Scikit-learn: Machine Learning in Python"]

[Sentence 137: "From the carrier sets of these calls, janus builds a neighbor-joining tree of the tumor samples and re-decides the carriers of each SV with it: every side of every edge of the tree is a candidate carrier set, scored by the summed evidence of its members ($\operatorname{logit} p_k - \operatorname{logit} \tau_k$ for a tumor $k$ with a record of the SV, $-\mu$ with $\mu=0.5$ otherwise) plus a prior from the number of SVs whose carriers are exactly that set; if the best candidate scores above zero, its members without an accepted record become carriers of the SV, and carriers outside it are dropped (Section REF)." -> cites: Saitou 1987, "The neighbor-joining method: a new method for reconstructing phylogene"]


## results

[Sentence 147: "For comparison, we ran competing PacBio HiFi-based callers using minimap2 alignments generated from identical read sets (Section REF) and constructed consensus truth sets using Minda (Section REF)." -> cites: Li 2018, "Minimap2: pairwise alignment for nucleotide sequences"; Lab 2024, "Minda: a tool for benchmarking and merging structural variant calls"]

[Sentence 149: "To assess significance, we computed 95\% bootstrap confidence intervals for the F1 score differences between janus and the best other caller across 2,000 resamples (sampling consensus SVs with replacement for recall and candidate calls for precision)." -> cites: Efron 1979, "Bootstrap methods: Another look at the jackknife"]

[Sentence 151: "We first evaluated janus on five PacBio HiFi tumor--normal cell line pairs from the CASTLE panel (Table REF)." -> cites: Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"; Park 2026, "Accurate somatic small variant discovery for multiple sequencing techn"]

[Sentence 152: "Lacking a curated ground truth, the panel provides an 11-call-set Minda consensus requiring support from at least two technologies and four callers across PacBio HiFi, ONT, and Illumina." -> cites: Lab 2024, "Minda: a tool for benchmarking and merging structural variant calls"; Keskus 2026, "Severus detects somatic structural variation and complex rearrangement"; Shiraishi 2023, "Precise characterization of somatic complex structural variations from"; Chen 2016, "Manta: rapid detection of structural variants and indels for germline "; Wala 2018, "SvABA: genome-wide detection of structural variants and indels by loca"; Cameron 2021, "GRIDSS2: comprehensive characterisation of somatic structural variatio"; Smolka 2024, "Detection of mosaic and population-level structural variants with Snif"; Elrick 2025, "SAVANA: reliable analysis of somatic structural variants and copy numb"]

[Sentence 153: "Because that consensus includes only four of the ten HiFi callers evaluated here, we expanded the truth set to 18 call sets by incorporating janus, our runs of all ten HiFi callers (Severus, SAVANA, nanomonsv, Sniffles2, dysgu, SVision-pro, cuteSV, DeBreak, kled, and SVDSS; Section REF), and the panel's seven ONT and Illumina sets." -> cites: Cleal 2022, "Dysgu: efficient structural variant calling using short or long reads"; Wang 2025, "De novo and somatic structural variant discovery with SVision-pro"; Jiang 2020, "Long-read-based human genomic structural variation detection with cute"; Chen 2023, "Deciphering the exact breakpoints of structural variations using long "; Zhang 2024, "Kled: an ultra-fast and sensitive structural variant detection tool fo"; Denti 2023, "SVDSS: structural variation discovery in hard-to-call genomic regions "]

[Sentence 163: "To assess generalization, we evaluated two independent benchmarks outside the training set: COLO829 (melanoma, 62 validated SVs) and HG008 (pancreatic tumor, GIAB draft v0.5, 172 SVs)." -> cites: Espejo Valle-Inclan 2022, "A multi-platform reference for somatic structural variation detection"; McDaniel 2025, "Development and extensive sequencing of a broadly-consented Genome in "; Standards 2026, "HG008-T somatic structural variant and copy number draft benchmark V0."]

[Sentence 174: "Due to the lack of multi-sample long-read tumor--normal cohorts with established ground truth, we evaluated joint calling on simulated data: one normal genome with real germline variants and 20 lineage-related tumor clones that inherit ancestral somatic SVs (Section REF)." -> cites: Battistella 2025, "insilicoSV: a flexible grammar-based framework for structural variant "; Kelleher 2016, "Efficient Coalescent Simulation and Genealogical Analysis for Large Sa"; Wick 2019, "Badread: simulation of error-prone long reads"; Logsdon 2025, "Complex genetic variation in nearly complete human genomes"]

[Sentence 187: "The first-pass calls at 20/5$\times$ yielded 5,873 unique SVs (1,846 false positives, with 1,730 restricted to single tumors), and a neighbor-joining tree built on their pairwise SV differences recovered all 17 simulated splits (Robinson--Foulds distance 0; Fig. REFd), as did the trees reconstructed at the other four depth settings (Table REF)." -> cites: Saitou 1987, "The neighbor-joining method: a new method for reconstructing phylogene"; Robinson 1981, "Comparison of phylogenetic trees"]

[Sentence 198: "Severus and SAVANA each spent 56 min of their run time on phasing using Clair3 and LongPhase (Table REF, Fig. REFa)." -> cites: Zheng 2022, "Symphonizing pileup and full-alignment for deep learning-based long-re"; Lin 2022, "LongPhase: an ultra-fast chromosome-scale phasing algorithm for small "]


## discussion

[Sentence 215: "Spatially resolved DNA sequencing could locate the detected SVs and the subclones carrying them in the tissue." -> cites: Zhao 2022, "Spatial genomics enables multi-modal study of clonal heterogeneity in "]

[Sentence 216: "If pangenome-based methods scale to large cohorts, fast and accurate somatic SV calling could become a routine part of precision medicine, guiding the choice and development of targeted drugs." -> cites: Schram 2017, "Fusions in solid tumours: diagnostic strategies, targeted therapy, and"]

