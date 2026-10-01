#!/usr/bin/env bash
# svscope.sh <tumor.bam> <normal.bam> <outdir> <threads>: SVScope (Goatofmountain/SVscope @ 0ad3686, image
# snape-bench/svscope:0ad3686, images/svscope/Dockerfile), DataPrepare --selectwindows --FullProcess, platform PacBio,
# MSA spoa (the README warns that the abPOA default can hang), its bundled hg38 RepeatMasker / 10 kb window beds and
# its shipped random forest (src/Model_pb.joblib) -> calls.vcf = PASS records of tumor.mergedSomatic.vcf.
# NOTE: that random forest was trained on labelled somatic SVs of HCC1395, H2009, H1437, HCC1954, HCC1937 (leave one
# cell line out in its paper), i.e. on our benchmark pairs; scores on those pairs are in-sample for SVScope.
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/svscope:0ad3686
dk $IMG python /opt/SVScope/src/SVScope.py --verbose DataPrepare -T $T -N $N -t tumor -n normal -r $REF -s $O/svscope -p $P \
   --selectwindows --FullProcess -M spoa --platform PacBio
pass $O/svscope/tumor.mergedSomatic.vcf $O/calls.vcf
done_msg
