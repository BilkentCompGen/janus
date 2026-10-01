# Competitor callers: wrappers

Every competitor runs through one interface:

```
bench/callers/<tool>.sh <tumor.bam> <normal.bam> <outdir> <threads>
  -> <outdir>/calls.vcf   the tool's somatic SV calls (PASS records as the tool defines them)
  -> <outdir>/log.txt     everything the tool printed
```

- Reference: `ref/hg38_parmask.fa` for every tool; VNTR bed `ref/vntr.trf.bed` where a tool takes one (Severus,
  Sniffles2, nanomonsv). Calls are in GRCh38 coordinates for every tool.
- BAMs: sorted + indexed, aligned as `hifi/align.sh` does (minimap2 `-ax map-hifi -Y --secondary=no`), with distinct
  `@RG SM` tags for tumor and normal (Delly needs them).
- `common.sh` (sourced by every wrapper): re-executes the wrapper under `taskset -c 16-31`, caps threads at 16, sends
  all output to `log.txt`, and provides `dk <image> <cmd>` = `docker run --rm --cpuset-cpus 16-31 --memory 40g
  -e OMP_NUM_THREADS=<threads> -u <uid>:<gid>` with the repo, the BAM directories and the outdir mounted at their own
  paths. Docker containers do not inherit the host affinity, hence `--cpuset-cpus`.
- `subtract.py`: for single-sample mosaic callers. It keeps the tumor's PASS SVs (>= 50 bp, optionally only
  `SOMATIC`-tagged ones) that have no normal SV of the same SVTYPE within 1 kb.
- `run_micro.sh <tool> ...`: runs wrappers on a validation pair with 16 threads, timed with `/usr/bin/time`, and
  samples `docker stats` for the containers' memory. It appends `tool exit wall_s peak_kb n_calls` to
  `micro/results.tsv`, or to `micro/demo/results.tsv` with `PAIR=demo`.
- Nothing is installed on the host. Tools come from the existing venvs in `lib/callers/` or from docker images.
  Our own images are built from `images/<tool>/Dockerfile`. Extra references are fetched by `ref/fetch.sh`.

## Validation pairs

- `micro/` (`micro/make.sh`): the first 131,000 FASTQ lines (32,750 reads, about 0.5 Gb, about 0.16x) of
  `hifi/hcc1954/{n,t}.fastq.gz`, aligned like `hifi/align.sh`. Alignment took 55 s for the normal and 66 s for the
  tumor, at 18 GB. At this depth almost nothing can be somatic, so this pair only shows that a tool runs end to end.
- `micro/demo/` (`micro/make_demo.sh`): SomaSV's bundled HG008 chr20 HiFi demo (2,283 tumor reads at 30x and 357
  normal reads at 10x over one region with somatic SVs), turned back into FASTQ and re-aligned the same way. It has
  real signal, so it also exercises the output stages.

## Status

Times are wall seconds with 16 threads on cores 16-31. "Peak" is the larger of GNU time's max RSS for the wrapper
tree and the sampled `docker stats` peak.

| Tool | Wrapper | Source | micro: status / wall / peak / calls | demo: status / wall / calls |
|---|---|---|---|---|
| Severus 1.7 (unphased) | `severus.sh` | venv `lib/callers/severus` | ok / 6.7 s / 0.6 GB / 0 | ok / 4.9 s / 4 |
| Severus 1.7 (phased: Clair3 + LongPhase) | `severus_phased.sh` | `hkubal/clair3:v2.0.3` + venv | ok / 283 s / 6.2 GB / 0 | ok / 46 s / 4 |
| Sniffles2 2.8.0 (multisample, SUPP_VEC=01) | `sniffles2.sh` | venv `lib/callers/sniffles` | ok / 9.8 s / 1.4 GB / 30 | ok / 6.9 s / 0 |
| Sniffles2 2.8.0 `--mosaic` minus normal | `sniffles2_mosaic.sh` | venv | ok / 8.6 s / 1.4 GB / 0 | ok / 5.0 s / 0 |
| SAVANA 1.3.8 (`--pb`) | `savana.sh` | venv `lib/callers/savana` | ok / 16.9 s / 14.8 GB / 2 | ok / 6.7 s / 4 |
| nanomonsv 0.9.0 (`--qv25`) | `nanomonsv.sh` | venv `lib/callers/nanomonsv` + host racon | ok / 20.2 s / 0.1 GB* / 0 | ok / 2.9 s / 4 |
| Delly 2.6.0 `lr -y pb` + `filter -f somatic` | `delly.sh` | `quay.io/biocontainers/delly:2.6.0--h3752d28_0` | ok / 38.9 s / 0.9 GB / 1 | ok / 2.6 s / 0 (7 joint calls, none pass the somatic filter: tumor coverage < 10) |
| SVision-pro 2.5 (somatic) | `svision-pro.sh` | `snape-bench/svision-pro:2.5-pysam0.21` | ok / 6.1 s / 0.8 GB / 3 | ok / 9.8 s / 2 |
| SomaSV 0.0.1 (long-read-only) | `somasv.sh` | `snape-bench/somasv:0.0.1-92f1560` | ok / 34.3 s / 8.4 GB / 0 | ok / 9.0 s / 0 (on the original demo BAMs: 4 calls, identical to its shipped `filtered.only_long.somatic.vcf`) |
| minisv r29 (hg38 + CHM13 + CHM13 graph, no normal assembly) | `minisv.sh` | `snape-bench/minisv:r29` | ok / 416 s (incl. 41 s CHM13 index build) / 28.6 GB (minigraph) / 0 | ok / 127 s / 3 |
| longcallD 0.0.11 `-s` minus normal | `longcalld.sh` | `quay.io/biocontainers/longcalld:0.0.11--hbc58adc_0` | ok / 41 s / 0.5 GB / 0 | ok / 3.3 s / 0 |
| SVScope @0ad3686 | `svscope.sh` | `snape-bench/svscope:0ad3686` | **fails** at the last stage / 292 s / 1.9 GB | **fails** / 192 s |

\* GNU time reports the largest single process, not the sum over nanomonsv's worker processes.

### Fixes and caveats found while validating

- **SVision-pro**: the biocontainer ships pysam 0.9.1. SVision-pro 2.5 then crashes while writing breakend records
  (`record.info.pop`: `AttributeError`, `src/output_vcf.py:194`). It hit this on the demo pair. Our image adds pysam
  0.21.0 (`images/svision-pro/Dockerfile`).
- **SomaSV** writes no VCF when it has 0 calls. The wrapper then writes a header-only `calls.vcf`, but only if the log
  says `Final somatic variants: 0`. On the original demo BAMs it reproduces its shipped long-read-only output exactly.
  On our re-aligned copy of the same reads it calls 0: 5 merged breakpoints instead of 7, from small alignment
  differences.
- **SVScope** runs window selection, split-alignment SV calling (`InterALNSVs.vcf`) and the local-graph step. It then
  crashes in its feature/random-forest stage whenever no inner-alignment candidate survives: a pandas `ValueError` on
  micro and an `EmptyDataError` on demo. Both test pairs are too small to produce one. Whether it completes on full
  data is **not yet verified**. Its random forest was trained on five of our benchmark cell lines, so any score of it
  on those pairs is in-sample (see COMPETITORS.md).
- **Delly**'s somatic filter needs tumor coverage >= 10 at the site (`-v 10`), so it reports 0 on both small pairs.
- **minisv**: minigraph on the 464-haplotype CHM13 graph peaks at 28.6 GB. On full data the realignment is three
  extra whole-genome alignments of the tumor and one of the normal, so expect hours. The hifiasm assembly of the
  normal from the README's full recipe is left out (> 40 GB).
- **Sniffles2** has no tumor/normal mode. `sniffles2.sh` is how the CASTLE panel ran it (multisample, tumor-only
  SUPP_VEC); `sniffles2_mosaic.sh` and `longcalld.sh` rely on our normal subtraction (`subtract.py`).
- **Severus**: `severus.sh` (unphased) is what our CASTLE reproduction used (`repro.txt`: 0.97 F1 against the
  published Severus calls). `severus_phased.sh` is the README's recommended phased setup.

## Exact commands (`$T` tumor, `$N` normal, `$O` outdir, `$P` threads, `$REF`, `$VNTR`)

- **severus.sh**:
  `severus --target-bam $T --control-bam $N --out-dir $O/severus --vntr-bed $VNTR -t $P`
  -> PASS of `severus/somatic_SVs/severus_somatic.vcf`
- **severus_phased.sh**:
  1. `run_clair3.sh --bam_fn=$N --ref_fn=$REF --threads=$P --platform=hifi --model_path=/opt/models/hifi_revio --output=$O/clair3 --enable_phasing --longphase_for_phasing`
  2. `longphase haplotag -s clair3/phased_merge_output.vcf.gz -b {$N,$T} -r $REF -t $P --tagSupplementary -o {n,t}.hp`, then `samtools index` both
  3. `severus --target-bam t.hp.bam --control-bam n.hp.bam --phasing-vcf clair3/phased_merge_output.vcf.gz --out-dir $O/severus --vntr-bed $VNTR -t $P`
  -> PASS of `severus_somatic.vcf`
- **sniffles2.sh**:
  `sniffles --input {$N,$T} --snf {n,t}.snf --threads $P --reference $REF --tandem-repeats $VNTR`, then
  `sniffles --input n.snf t.snf --vcf joint.vcf --threads $P --reference $REF`
  -> PASS records with `SUPP_VEC=01`
- **sniffles2_mosaic.sh**:
  `sniffles --input $T --vcf t.mosaic.vcf --mosaic ...` and `sniffles --input $N --vcf n.vcf --minsupport 1 ...` (both with the reference and tandem repeats)
  -> `subtract.py t.mosaic.vcf n.vcf --dist 1000`
- **savana.sh**:
  `savana --tumour $T --normal $N --ref $REF --outdir $O/savana --threads $P --pb --sample s --overwrite`
  -> PASS of `s.classified.somatic.vcf`
- **nanomonsv.sh**:
  `nanomonsv parse {$T,$N} {t/t,n/n} --reference_fasta $REF`, then
  `nanomonsv get t/t $T $REF --control_prefix n/n --control_bam $N --processes $P --qv25 --simple_repeat_bed repeat.bed.gz` (`repeat.bed.gz` = the VNTR bed, bgzipped and tabixed)
  -> PASS of `t/t.nanomonsv.result.vcf`
- **delly.sh**:
  `delly lr -y pb -g $REF -x human.hg38.excl.tsv -o joint.bcf $T $N`, then
  `delly filter -f somatic -o somatic.bcf -s samples.tsv joint.bcf` (`samples.tsv`: tumor SM -> tumor, normal SM -> control)
  -> `bcftools view -f PASS`
- **svision-pro.sh**:
  `SVision-pro --target_path $T --base_path $N --genome_path $REF --model_path model_liteunet_256_8_16_32_32_32.pth --access_path hg38.access.10M.bed --out_path $O/svp --sample_name tumor --preset hifi --detect_mode somatic --process_num $P/2 --unet_cpu_num 2 --device cpu`, then
  `python extract_op.py --input_vcf svp/tumor.svision_pro_v2.5.s5.vcf --extract somatic`
  -> PASS of `*.somatic_s1.vcf`
- **somasv.sh**:
  `python main.py --tumor-bam $T --normal-bam $N --reference $REF --output-dir $O/somasv --sample-platform HIFI --mode long-read-only --processes $P`
  -> PASS of `final_somatic_variants.vcf`
- **minisv.sh**:
  1. `samtools fastq -F 0x900 $T > t.fq.gz`
  2. `minimap2 -cx map-hifi -s50 --ds -t $P {hg38_parmask.mmi, hs1.map-hifi.mmi} t.fq.gz`
  3. `minigraph -cxlr -t $P CHM13-464.gfa.gz t.fq.gz`
  4. `samtools fastq $N | minimap2 -cx map-hifi -s50 --ds hg38_parmask.mmi -`
  5. `minisv.js e -n TUMOR -b hs38.cen-mask.bed t.hg38l t.chm13l t.chm13g | bash > t.rsv`, and `minisv.js extract -n NORMAL -b hs38.cen-mask.bed n.hg38l > n.rsv`
  6. `sort -k1,1 -k2,2n | minisv.js merge - | grep TUMOR | grep -v NORMAL | minisv.js genvcf`
  -> `calls.vcf`
- **longcalld.sh**:
  `longcallD call -s -t $P $REF $T > t.mosaic.vcf` and `longcallD call -t $P $REF $N > n.vcf`
  -> `subtract.py --somatic-tag --dist 1000`
- **svscope.sh**:
  `python src/SVScope.py --verbose DataPrepare -T $T -N $N -t tumor -n normal -r $REF -s $O/svscope -p $P --selectwindows --FullProcess -M spoa --platform PacBio`
  -> PASS of `tumor.mergedSomatic.vcf`

## Files

- `COMPETITORS.md`: the literature/tool survey (which tools exist, which cannot run and why).
- `images/{minisv,somasv,svscope,svision-pro}/Dockerfile`: our images. Pulled images: delly, svision-pro, longcalld
  (biocontainers) and `hkubal/clair3:v2.0.3`.
- `ref/`: `fetch.sh`; `hs1.fa.gz` + `hs1.map-hifi.mmi` (T2T-CHM13v2.0, UCSC); `CHM13-464.gfa.gz` (HPRC minigraph
  graph, zenodo 16728828); `human.hg38.excl.tsv` (Delly); `svision-pro/` (v2.5 repo: model + access bed); `somasv/`
  and `svscope/` (source checkouts used for the images and the demo data).
- Older files kept as they were: `run.sh`, `sniffles_somatic.sh`, `repro.sh`, `repro.txt`, `times.txt`, and the full
  CASTLE runs in `severus/` and `sniffles/`.
