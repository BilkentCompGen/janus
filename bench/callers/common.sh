set -euo pipefail
[ -z "${PINNED:-}" ] && PINNED=1 exec taskset -c ${CPUS:-16-31} "$0" "$@"
[ $# -eq 4 ] || { echo "usage: $(basename $0) <tumor.bam> <normal.bam> <outdir> <threads>" >&2; exit 2; }
T=$(readlink -f "$1"); N=$(readlink -f "$2"); mkdir -p "$3"; O=$(readlink -f "$3"); P=$4
[ "$P" -gt 16 ] && P=16
REF=/home/yuio/snape/ref/hg38_parmask.fa; VNTR=/home/yuio/snape/ref/vntr.trf.bed
L=/home/yuio/snape/lib/callers; K=/home/yuio/snape/bench/callers
CPUS=${CPUS:-16-31}
exec > "$O/log.txt" 2>&1
echo "[$(date '+%F %T')] $(basename $0) tumor=$T normal=$N out=$O threads=$P"
dk() {
  local img=$1; shift
  local m=() d
  mkdir -p "$O/.tmp"
  for d in /home/yuio/snape $(printf "%s\n" "$(dirname $T)" "$(dirname $N)" "$O" | grep -v '^/home/yuio/snape\(/\|$\)' | sort -u); do m+=(-v "$d:$d"); done
  docker run --rm --cpuset-cpus $CPUS --memory ${DMEM:-40g} --memory-swap ${DMEM:-40g} --oom-score-adj 1000 ${GPU:+--gpus all} ${SHM:+--shm-size $SHM} ${CGPARENT:+--cgroup-parent=$CGPARENT} -e OMP_NUM_THREADS=$P -u "$(id -u):$(id -g)" -e HOME=/tmp -v "$O/.tmp:/tmp" -w "$O" "${m[@]}" "$img" "$@"
}
pass() { awk -F'\t' '/^#/ || $7=="PASS" || $7=="."' "$1" > "$2"; }
done_msg() { echo "[$(date '+%F %T')] done: $(grep -vc '^#' $O/calls.vcf) calls in $O/calls.vcf"; }
