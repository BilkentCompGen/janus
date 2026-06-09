<!-- Deep-research report (workflow wq9lq0ixj, 2026-06-08).
Recommendation: Adopt Sim-it as the primary HiFi read engine for somatic-SV benchmarking — but mixing real HiFi reads is what the somatic-SV literature actually does and is the stronger move if full in-silico isn't mandatory. -->

# Replacing pbsim3 for PacBio HiFi somatic-SV simulation

**Bottom line:** For benchmarking *somatic SV calling on PacBio HiFi*, switch the primary read engine from **pbsim3 (+ external `ccs`)** to **[Sim-it](https://github.com/ndierckx/Sim-it)** — the only actively-maintained tool that ships a dedicated *Sequel CCS/HiFi* error profile, simulates SVs natively (single-bp truth), and won the indel-size-distribution category in the one peer-reviewed long-read SV benchmark. Keep pbsim3+ccs as a read-level-marginals sanity check, and reserve **HI.SIM** / **Badread** as escalation options for context/error-structure stress tests. **Strategic caveat (read §6):** the modern somatic-SV literature overwhelmingly does *not* simulate reads de-novo — it mixes *real* HiFi reads — so if the project is not contractually fully in-silico, real-read mixing is the more defensible benchmark.

---

## 1. Why pbsim3 may be discouraged for HiFi somatic-SV work

The realism axis that actually gates your caller is **technology-specific error STRUCTURE** (indel-dominant, homopolymer-concentrated, context-dependent) plus an unambiguous normal-vs-tumor truth set. pbsim3 is weakest on exactly those axes for HiFi.

1. **Its error model is explicitly sequence-INDEPENDENT.** The PBSIM3 paper states verbatim: *"Our error model is sequence independent, whereas the k-mer profile is sequence dependent,"* and *"Both PBSIM3 models are sequence independent; therefore, homopolymer bias cannot be simulated"* (the lone `--hp-del-bias` exception is ONT-only and "only slightly" simulates the bias). Real HiFi residual error is the opposite — indel-dominant and homopolymer-concentrated: of CCS read discordances, **~92.0% are homopolymer indels**, 4.6% non-homopolymer indels, 3.4% mismatches; the paper's own real-HiFi (CHM13) breakdown is **0.13% del / 0.10% ins / 0.02% sub**. Any homopolymer realism in pbsim3's HiFi output comes from the *external ccs consensus step + the sequence-dependent k-mer profile*, **not** from the error model — the paper confirms simulated HiFi "showed the same characteristics as real HiFi reads, with longer homopolymer lengths resulting in higher error rates" precisely *because* the CLR output is fed to ccs. For a caller whose precision is gated by systematic/context-dependent residual error (your "phantom-bubble" residuals), depending on an incidental side-effect rather than a calibrated model is the single biggest reason to switch. ([paper](https://academic.oup.com/nargab/article/4/4/lqac092/6855700), [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9713900/), [CCS-accuracy/DeepConsensus](https://www.biorxiv.org/content/10.1101/2021.06.28.450238v1.full))

2. **"pbsim3 HiFi" is really pbsim3 + external `ccs` — a brittle two-tool handoff.** PBSIM3 emits only multi-pass CLR subread BAM; you must run PacBio's `ccs` (paper used **v6.0.0**, current **v6.4.0**) to get consensus reads. The README states it directly: *"PBSIM does not directly simulate HiFi reads, but only simulate generation of CLR by multi-pass sequencing; the output … is input to ccs software, which generates HiFi reads."* This handoff is documented as fragile — [GitHub #41](https://github.com/yukiteruono/pbsim3/issues/41) shows a parameter mismatch (an `ERRHMM-SEQUEL.model` passed to the `qshmm` method) driving ccs to *"ZMWs pass filters: 81 (0.103%)"* with 98% draft-generation error — i.e. near-zero HiFi reads. ([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9713900/), [#41](https://github.com/yukiteruono/pbsim3/issues/41))

3. **Maintainer-admitted model flaw: read length is decoupled from pass count / quality.** Real HiFi gives longer inserts fewer passes → lower accuracy; pbsim3 uses one fixed `--pass-num` for all reads. The author confirmed this in [GitHub #25](https://github.com/yukiteruono/pbsim3/issues/25) (author_association = OWNER): *"the longer the DNA fragment, the fewer the number of passes, so the longer the read, the lower the quality … However … it is currently not possible to implement this relationship in PBSIM3."* The 2025 **MHASS** tool ([Bioinformatics btaf656](https://academic.oup.com/bioinformatics/article/42/1/btaf656/8373413) / [arXiv 2509.12428](https://arxiv.org/abs/2509.12428)) independently confirms the artifact — it exists specifically to patch *"its restriction to a single pass number per run"* using empirical/lognormal pass-number distributions, citing the PBSIM3 *paper* (not the issue) for the limitation.

4. **Known deletion mis-shaping — bad for a DEL-sensitive caller.** The paper reports pbsim3 OVER-produces deletions across *all four* tested methods: *"Compared to the error rate of real HiFi reads, the deletion rate was higher for all four simulation methods … In particular, the deletion rate was high in the two error models,"* and the homopolymer slope is steeper than real: *"the rate of increase in the deletion rate of simulated HiFi reads was greater than that of real HiFi reads."* A simulator that systematically inflates short deletions in homopolymers biases the false-positive characterization of a DEL-sensitive somatic caller — and is consistent with this project's own finding that ccs/HiFi systematic del-errors drive phantom calls. ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9713900/) / [Oxford Academic](https://academic.oup.com/nargab/article/4/4/lqac092/6855700))

5. **Compute/scale/stability cost at whole-genome.** pbsim3 hard-rejects any sequence >1 Gbp — the source defines `#define REF_SEQ_LEN_MAX 1000000000` and errors *"Reference is too long. Acceptable length <= …"* ([src/pbsim.cpp](https://github.com/yukiteruono/pbsim3/blob/master/src/pbsim.cpp)) — a per-contig cap that forces splits interacting badly with ccs ZMW consensus. [GitHub #26](https://github.com/yukiteruono/pbsim3/issues/26) reports a ~1 Gbp genome at pass-num 20, 5–25X on 16 cores/128 GB *"taking days, if it works at all"* with crashes, and an open segfault remains at [#33](https://github.com/yukiteruono/pbsim3/issues/33). You run whole-genome 30× HiFi — exactly this regime. *(Caveats for accuracy: #26 is filed as a slowness, not segfault, report and is closed; #28's segfault is also closed — #33 is the open one. The frequently-cited "~5 TB intermediate SAM" figure is a plausible extrapolation but is **not** documented in any pbsim3 source.)*

**What pbsim3 IS good at (so the trade-off is honest):** the strongest third-party evidence — a May-2026 study of six simulators vs real ONT R10.4.1 — found pbsim3 most closely reproduced *read length, read accuracy, and FASTQ quality scores* (the best read-level all-rounder), but *"did not capture context-dependent substitution rates and homopolymer-length errors,"* and that study is **ONT, not HiFi**. So pbsim3 wins on the marginals you care about *least* for SV-precision work and loses on the error structure you care about *most*.

---

## 2. Recommended replacement: Sim-it

[Sim-it](https://github.com/ndierckx/Sim-it) is the best *maintained, HiFi-capable, SV-aware* single-tool option.

- **Native PacBio HiFi profile.** Ships `error_profile_PB_Sequel_CCS_hifi.txt` (plus ONT / RS II / Sequel II profiles) — one of the few tools covering the actual HiFi platform rather than ONT-only.
- **Built-in SV simulation with single-bp truth.** Takes a reference FASTA + config and an optional **SV-position VCF**, so the normal-vs-tumor pair *and* its ground truth come from one tool at single-bp resolution — no ccs ground-truth ambiguity, no separate spike-in step.
- **Best peer-reviewed indel realism for SV work.** In the dedicated long-read SV benchmark that introduced it (Dierckxsens et al., *Genome Biology* 2021), Sim-it produced *"the most accurate length distributions of deletion/insertion errors"* — the exact property that distinguishes a real SV from a clustered indel artifact.
- **Actively maintained.** v1.4.5 (May 2026), Apache-2.0.
- **Honest limitation:** in that same benchmark, *LongISLND* had the most accurate context-specific errors — but LongISLND is a pre-HiFi CLR tool (last commit 2017, no HiFi/CCS profile) and not a viable HiFi replacement. No simulator wins every metric; Sim-it is the best on the SV-relevant ones while remaining HiFi-aware and maintained.

---

## 3. Runner-ups (trade-offs)

| Tool | HiFi/CCS? | SV-aware? | Strength | Why it loses to Sim-it |
|---|---|---|---|---|
| **[HI.SIM](https://github.com/thegenemyers/HI.SIM)** (Gene Myers) | Yes (CCS-pass model) | No (pair with VISOR/insilicoSV) | **Most explicit error-structure architecture:** a 7-mer context table (target base ±3) for del/ins/sub, dedicated homopolymer/di-/tri-nucleotide microsatellite tables with length-dependent indel rates, and a length-stratified per-read error-rate distribution capturing the ccs pass-count accuracy spread. Use as an **error-structure stress test** if homopolymer/context realism turns out to dominate your caller's behavior. | No native SV-spiking; realism is architecture-based with **no independent head-to-head benchmark** proving it quantitatively closest. *(Correction vs the common claim: it is **not** dormant — last commit **2026-02-12**, [verified](https://api.github.com/repos/thegenemyers/HI.SIM/commits/main); and it is **not** "the only" tool to model HiFi error structure — pbsim3's ccs-pathway output also exhibits homopolymer-length-dependent error, just indirectly.)* |
| **[Badread](https://github.com/rrwick/Badread)** | Yes — real `pacbio2021` Sequel-IIe HiFi model | Indirectly (it's the read backend [VISOR](https://github.com/davidebolo1993/VISOR) pins) | Best for read **failure modes** (junk/glitches/chimeras/full identity spectrum) and fine user control; a strong *complement* for robustness testing. | README does not describe context-dependent error modeling, and the author states it *"does not try to be best at imitating real reads"* — weak as the error-structure-realism workhorse. |
| **pbsim3 + ccs** (incumbent) | Yes (via external ccs) | No (you pair with insilicoSV) | Best read-level marginals (length / accuracy / FASTQ Q); already wired into `simulation/scripts/sim_reads_wg_ccs.sh`. | All of §1 — keep only as a *marginals sanity check*, not the realism primary. |

**Explicitly unsuitable** (so the choice is defensible): **NanoSim / DeepSimulator** (ONT only; DeepSimulator academic-license-only), **PBSIM2** (removed CCS/HiFi; cannot distinguish error types), **PaSS / NPBSS / pbsim v1** (PacBio *CLR* only), **SimLoRD** (CCS-capable but dormant since 2016, simplified errors — no upgrade over pbsim3).

---

## 4. Practical integration — replacing the pbsim3→ccs pipeline

Your current path is `pbsim3` (multipass CLR → BAM) → `ccs` consensus, run per-FASTQ in `simulation/scripts/sim_reads_wg_ccs.sh`, producing `reads/normal/` and `reads/wgs/`, with `insilicoSV` supplying the tumor truth assembly/VCF.

**Drop-in with Sim-it (single tool, no ccs handoff):**

1. **Inputs:** germline reference FASTA (`simulation/ref.fasta`) + a Sim-it config + the **HiFi error profile** `error_profile_PB_Sequel_CCS_hifi.txt` + (optionally) an **SV-position VCF**.
2. **Normal arm:** run Sim-it on `ref.fasta` with no SVs (or germline-only VCF) at 30× → replaces `reads/normal/`.
3. **Tumor arm:** two equivalent options —
   - feed Sim-it your existing **insilicoSV tumor VCF** as its SV-position input so it implants the somatic SVs and emits reads in one pass, **or**
   - keep insilicoSV building the tumor *assembly* and run Sim-it over that assembly. Either way the truth is single-bp and the `ccs` step disappears entirely (no ZMW-consensus ground-truth ambiguity, no >1 Gbp split, no #41-style pass-rate collapse).
4. **Wiring:** mirror the per-arm structure of `sim_reads_wg_ccs.sh`; long WG runs still belong in a detached `screen -dmS … -L -Logfile …` per CLAUDE.md. Output FASTQ feeds the caller unchanged (`svp run reads/normal/… reads/wgs/… …`) — Sim-it changes only how the FASTQs are generated, not the caller contract (reads-only, whole-genome).
5. **Error model knobs:** the HiFi profile is the realism driver; verify coverage and indel-length spectra against your phantom-bubble QC before committing the WG run.

**If you escalate to HI.SIM** (error-structure stress test): generate the SV-bearing haplotypes with VISOR/insilicoSV first, then run HI.SIM over those haplotypes for the reads — it has no native SV-spiking, so it is strictly a two-tool path, but its explicit 7-mer/homopolymer tables make it the sharpest probe of whether your residual-error precision floor is context-driven.

---

## 5. Caveat to surface before committing — the field doesn't simulate reads de-novo

Worth flagging to your supervisor: modern somatic/germline long-read SV-caller benchmarks overwhelmingly do **not** use a read simulator. They use **real HG002 HiFi** against an assembly-derived truth set:

- **Germline:** Qin & Li 2025 (*GigaScience*) benchmarks cuteSV/pbsv/Sniffles2/SVDSS/longcallD/Sawfish/SVIM/DeBreak against **real HG002 HiFi** vs both **GIAB-SV v0.6** (9,705 SVs) and **HG002-Q100 v1.1** (28,188 SVs), no simulated reads ([GigaScience](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf154/8378076)). The "Hitchhiker's Guide" ([Biomedicines 2025](https://www.mdpi.com/2227-9059/13/8/1949)) and [kcleal's HG002 HiFi benchmark](https://github.com/kcleal/SV_benchmark_PacBio_HiFi) do likewise ("no simulated reads").
- **Somatic / low-VAF:** they MIX *real* reads. **Severus** ([Nat Biotechnol 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12483193/)) builds its synthetic tumor from **real CHM1 + CHM13 HiFi** (CHM1 = normal), truth from those assemblies. **Sniffles2** ([Nat Biotechnol 2023](https://www.nature.com/articles/s41587-023-02024-y), mosaic mode) titrates **real HG002 into real HG00733** at 7–28% VAF.
- **insilicoSV** — which your repo already uses for tumor truth — is one of the few SV simulators paired with a read simulator, and it pairs with **pbsim3**. *(Bounding the generalization honestly: de-novo simulation is not extinct — e.g. SVision-pro, [Nat Biotechnol 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11825360/), benchmarks somatic SVs with PBSIM-simulated reads + VISOR-implanted variants — but it is the minority practice among the major modern tools.)*

So two moves are strategically stronger than swapping simulators:
- **(a) Somatic axis:** generate the normal/tumor pair by **mixing real HiFi reads** (the field-standard "synthetic tumor"); and/or
- **(b) Germline validity:** **benchmark on real HG002 HiFi vs a T2T-Q100 truth set**, reserving a simulator (Sim-it) only for controlled SV-class/edge-case coverage where real truth is unavailable.

If the project must stay fully in-silico end-to-end, **Sim-it is the recommendation in §2**. But real-read mixing is what the somatic-SV literature actually does, and is the most defensible answer to *"why should I trust this benchmark."*

---

### Sources
- PBSIM3 paper: https://academic.oup.com/nargab/article/4/4/lqac092/6855700 · https://pmc.ncbi.nlm.nih.gov/articles/PMC9713900/
- pbsim3 source / issues: https://github.com/yukiteruono/pbsim3/blob/master/src/pbsim.cpp · [#25](https://github.com/yukiteruono/pbsim3/issues/25) · [#26](https://github.com/yukiteruono/pbsim3/issues/26) · [#41](https://github.com/yukiteruono/pbsim3/issues/41) · [#33](https://github.com/yukiteruono/pbsim3/issues/33)
- CCS/HiFi error structure: https://www.biorxiv.org/content/10.1101/2021.06.28.450238v1.full
- MHASS (pass-number patch): https://academic.oup.com/bioinformatics/article/42/1/btaf656/8373413 · https://arxiv.org/abs/2509.12428
- Sim-it: https://github.com/ndierckx/Sim-it (Genome Biology 2021 SV benchmark)
- HI.SIM: https://github.com/thegenemyers/HI.SIM (latest commit https://api.github.com/repos/thegenemyers/HI.SIM/commits/main)
- Badread: https://github.com/rrwick/Badread · VISOR: https://github.com/davidebolo1993/VISOR
- Real-read somatic/germline benchmarks: https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf154/8378076 · https://pmc.ncbi.nlm.nih.gov/articles/PMC12483193/ (Severus) · https://www.nature.com/articles/s41587-023-02024-y (Sniffles2) · https://github.com/kcleal/SV_benchmark_PacBio_HiFi · https://www.mdpi.com/2227-9059/13/8/1949 · https://pmc.ncbi.nlm.nih.gov/articles/PMC11825360/ (SVision-pro counterexample)
