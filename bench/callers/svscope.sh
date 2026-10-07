#!/usr/bin/env bash
source /home/yuio/snape/bench/callers/common.sh
IMG=snape-bench/svscope:0ad3686
dk $IMG python /opt/SVScope/src/SVScope.py --verbose DataPrepare -T $T -N $N -t tumor -n normal -r $REF -s $O/svscope -p $P \
   --selectwindows --FullProcess -M spoa --platform PacBio
pass $O/svscope/tumor.mergedSomatic.vcf $O/calls.vcf
done_msg
