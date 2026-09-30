#!/bin/bash
# convert_pause.sh stop|cont <sample...>: pauses (SIGSTOP) or resumes (SIGCONT) the conversion writers of the samples
# (fastq-dump for SRA runs, samtools fastq for HG008), so their output stops growing while the disk is short
sig=$1; shift
for s in "$@"; do
  for p in $(pgrep -f "fastq-dump --stdout --split-spot ${s}\.") $(pgrep -f "samtools fastq -@ 4 ${s}\."); do kill -$sig $p && echo "$sig $s pid $p"; done
done
