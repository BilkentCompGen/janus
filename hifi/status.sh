#!/usr/bin/env bash
# status.sh: download state of fetch_reads.sh. Files are preallocated, so their size says nothing; instead:
#   total speed = bytes received on the network interface over a 10 s sample,
#   per file    = aria2c's last progress line (raw/tmp/<file>.dl.log, written every 60 s),
# plus finished files, free disk, and whether the fetch and conversion processes are alive.
cd /home/yuio/snape/hifi/raw || exit 1
rx() { awk -F'[: ]+' 'NR>2 && $2!="lo"{s+=$3} END{print s}' /proc/net/dev; }
a=$(rx); sleep 10; b=$(rx); echo "network: $(( (b - a) / 10 / 1000000 )) MB/s (10 s sample)"
for l in tmp/*.dl.log; do
    [ -e "$l" ] || continue
    f=$(basename "$l" .dl.log)
    if [ -e "$f.aria2" ]; then
        p=$(grep -o '\[#[^]]*\]' "$l" | tail -1)
        printf "  %-50s %s\n" "$f" "${p:-starting}"
    elif [ -e "$f" ]; then printf "  %-50s downloaded\n" "$f"
    else printf "  %-50s converted, source in trash\n" "$f"; fi
done
for g in *.fastq.gz *.fastq.gz.tmp; do [ -e "$g" ] && printf "  %-50s %s %s\n" "$g" "$(du -BG --apparent-size "$g" | cut -f1)" "$([[ $g == *.tmp ]] && echo converting || echo done)"; done
echo "free $(df --output=avail -BG . | tail -1 | tr -d ' '); trash $(du -sBG /home/yuio/trash 2>/dev/null | cut -f1) (empty it to free space)"
echo "fetch_reads.sh procs: $(pgrep -f 'bash fetch_reads.sh' | wc -l); aria2c: $(pgrep -x aria2c | wc -l); converting: $(pgrep -f 'fastq-dump|samtools fastq' | wc -l)"
grep -v "downloading" fetch.log | tail -n 3
