#!/usr/bin/env bash
O=/home/yuio/snape/bench/runs/load.tsv
[ -s $O ] || echo -e "time\tload1\tall_cores\trerun_cores\tother_cores\tother_top" > $O
cpu() { awk '/^cpu /{print $2 + $3 + $4 + $7 + $8}' /proc/stat; }
own() { c=$(systemctl --user show snape-rerun.service -p CPUUsageNSec --value 2>/dev/null); case $c in ''|*[!0-9]*) echo 0;; *) echo $c;; esac; }
H=$(getconf CLK_TCK); a=$(cpu); r=$(own); t=$(date +%s.%N)
while true; do
    sleep 30
    b=$(cpu); s=$(own); u=$(date +%s.%N)
    mine=" $(cat /sys/fs/cgroup$(systemctl --user show snape-rerun.service -p ControlGroup --value)/cgroup.procs 2>/dev/null | tr '\n' ' ')"
    top=$(top -b -n 2 -d 2 -o %CPU -w 200 | awk -v m="$mine" '/^top -/{k++} k == 2 && $1 ~ /^[0-9]+$/ && $9 > 50 && index(m, " " $1 " ") == 0 {printf "%s:%.0f ", $12, $9}' | cut -c1-80)
    awk -v a=$a -v b=$b -v r=$r -v s=$s -v t=$t -v u=$u -v h=$H -v l="$(cut -d' ' -f1 /proc/loadavg)" -v top="$top" -v d="$(date '+%F %T')" \
        'BEGIN { all = (b - a) / h / (u - t); own = (s - r) / 1e9 / (u - t); if (own < 0) own = 0; printf "%s\t%s\t%.1f\t%.1f\t%.1f\t%s\n", d, l, all, own, all - own, top }' >> $O
    a=$b; r=$s; t=$u
done
