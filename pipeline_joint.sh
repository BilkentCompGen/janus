#!/bin/bash
# Joint run: one normal and several tumor samples of a patient. usage: pipeline_joint.sh DIR normal.fq tumor1.fq [tumor2.fq ...]
D=$1; shift; N=$(realpath $1); shift; TS=(); for f in "$@"; do TS+=($(realpath $f)); done; K=${#TS[@]}; P=${P:-16}
J=/home/yuio/janus; B=$J/build/janus; R=$J/ref; PY=${PY:-$J/.venv-ml/bin/python}
mkdir -p $D; cd $D
step() { n=$1; shift; [ -e done.$n ] && return; /usr/bin/time -f "$n %e s %M KB %U user %S sys" "$@" > log.$n 2>&1 || exit 1; tail -1 log.$n >> steps.txt; touch done.$n; }
S=(n); for k in $(seq 1 $K); do S+=(t$k); done
step pack_n $B pack -i $N -o n.svb -t $P
for k in $(seq 1 $K); do step pack_t$k $B pack -i ${TS[$((k-1))]} -o t$k.svb -t $P; done
for s in ${S[@]}; do step extract_$s $B extract -i $s.svb -o $s.svb.exc -t $P; done
step count_n $B count -i n.svb -o n.cnt -t $P -m 0
for k in $(seq 1 $K); do step count_t$k $B count -i t$k.svb -o t$k.cnt -t $P -d n.cnt; done
MS=$(for k in $(seq 1 $K); do python3 -c "print(open('t$k.cnt', 'rb').read(27)[26])"; done | paste -sd,); M=$(echo $MS | tr , '\n' | sort -n | head -1)
SVB=$(printf "%s.svb," ${S[@]} | sed 's/,$//'); CNT=$(printf "%s.cnt," ${S[@]} | sed 's/,$//')
step graph $B graph -i $CNT -o graph.gph -m $M
step prune $B prune -i graph.gph -r $SVB -o prune -t $P -m $M
step select $B select -i graph.gph -d prune.sites -r $SVB -o ev -t $P
step place $B place -i ev.fa -r $R/${REF:-hg38_parmask}.cidx -d $R/${REF:-hg38_parmask}.fa -o ev -t $P
step local minimap2 -cx map-hifi --secondary=no -t $P -o local.paf ev.win.fa ev.on.fa
step lift $B lift -i local.paf -d $R/${REF:-hg38_parmask}.fa -o aln.paf
step clip $B clip -i ev -r aln.paf -o ev.gaps.fa
step gaps minimap2 -cx sr --secondary=no -t $P -o ev.gaps.raw.paf $R/${REF:-hg38_parmask}.sr.mmi ev.gaps.fa
step patch $B patch -i ev -r ev.gaps.raw.paf -d aln.paf -o ev.paf
step call $B call -i ev -r $SVB -d $R/${REF:-hg38_parmask}.fa -g graph.gph -o calls.vcf -v $R/$( [ -n "$REF" ] && echo $REF.trf.bed || echo vntr.trf.bed ) -t $P -m $M
HS=$(for k in $(seq 1 $K); do printf "t%s.cnt.hist," $k; done | sed 's/,$//')
step joint env PYTHONPATH=$J/python $PY -m janus.joint -i calls.vcf -o janus.vcf -m $J/models/janus.pkl -H $HS -s $MS
