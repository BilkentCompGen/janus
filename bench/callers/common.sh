# common.sh: sourced by every competitor wrapper <tool>.sh <tumor.bam> <normal.bam> <outdir> <threads>.
# Sets T (tumor BAM), N (normal BAM), O (outdir), P (threads, capped at 16), REF, VNTR, L (venvs), and
# dk <image> <cmd...>: run a command in a docker image on cores $CPUS (16-31), <= $DMEM (40 GB), in cgroup slice
# $CGPARENT when set (so the benchmark can read the containers' CPU time and peak memory), as the calling user, with the
# inputs, the outdir and the repo mounted at their own paths. All output goes to $O/log.txt.
set -euo pipefail
[ -z "${PINNED:-}" ] && PINNED=1 exec taskset -c ${CPUS:-16-31} "$0" "$@"   # everything on cores $CPUS (default 16-31)
[ $# -eq 4 ] || { echo "usage: $(basename $0) <tumor.bam> <normal.bam> <outdir> <threads>" >&2; exit 2; }
T=$(readlink -f "$1"); N=$(readlink -f "$2"); mkdir -p "$3"; O=$(readlink -f "$3"); P=$4
[ "$P" -gt 16 ] && P=16
REF=/home/yuio/snape/ref/hg38_parmask.fa; VNTR=/home/yuio/snape/ref/vntr.trf.bed
L=/home/yuio/snape/lib/callers; K=/home/yuio/snape/bench/callers
CPUS=${CPUS:-16-31}     # the benchmark (bench/compete.sh) sets CPUS=0-15, DMEM and CGPARENT
exec > "$O/log.txt" 2>&1
echo "[$(date '+%F %T')] $(basename $0) tumor=$T normal=$N out=$O threads=$P"
dk() {
  local img=$1; shift
  local m=() d
  for d in /home/yuio/snape $(printf "%s\n" "$(dirname $T)" "$(dirname $N)" "$O" | grep -v '^/home/yuio/snape\(/\|$\)' | sort -u); do m+=(-v "$d:$d"); done
  docker run --rm --cpuset-cpus $CPUS --memory ${DMEM:-40g} ${CGPARENT:+--cgroup-parent=$CGPARENT} -e OMP_NUM_THREADS=$P -u "$(id -u):$(id -g)" -e HOME=/tmp -w "$O" "${m[@]}" "$img" "$@"
}
# pass <in.vcf> <out.vcf>: header + FILTER PASS (or '.') records
pass() { awk -F'\t' '/^#/ || $7=="PASS" || $7=="."' "$1" > "$2"; }
done_msg() { echo "[$(date '+%F %T')] done: $(grep -vc '^#' $O/calls.vcf) calls in $O/calls.vcf"; }
