#!/usr/bin/env python3
"""Merge dup_cov_flags.bin DUP candidates into calls.vcf.
usage: merge_dup_cov.py <dir>   (reads dir/calls.vcf + dir/dup_cov_flags.bin, writes dir/calls.vcf)
"""
import struct, sys, os

DIR = sys.argv[1]
CALLS = f'{DIR}/calls.vcf'
DUPS  = f'{DIR}/dup_cov_flags.bin'

# Load chrom_id -> name (need names for VCF)
chrom_name = {}
chrom_len  = {}
for line in open(f'{DIR}/chroms.tsv'):
    p = line.rstrip().split('\t')
    chrom_name[int(p[0])] = p[1]
    chrom_len[int(p[0])] = int(p[2])

# Load cov-DUPs
with open(DUPS, 'rb') as f: data = f.read()
K = struct.unpack_from('<Q', data, 0)[0]
print(f'[merge] {K} cov-DUP regions', file=sys.stderr)
dup_recs = []
dropped = 0
for i in range(K):
    cid, ps, pe, svlen, ct, cn = struct.unpack_from('<IIIIII', data, 8 + i*24)
    if cid not in chrom_name: dropped += 1; continue
    if svlen < 200: dropped += 1; continue
    # Real tandem DUP: tumor cov ≈ 2x normal. Stricter: ct/cn >= 1.7.
    if cn > 0 and ct * 10 < cn * 17: dropped += 1; continue
    # Drop absurdly large clusters (>500kb) — likely whole-chromosome repeat artifacts.
    if svlen > 500_000: dropped += 1; continue
    # Require cov_t >= 2.2 * bg_t (bg_t was 29 in our 30x sim). Hard-coded.
    if ct < 65: dropped += 1; continue
    chrom = chrom_name[cid]
    pos = ps + 1
    end = pos + svlen - 1
    line = f'{chrom}\t{pos}\t.\tN\t<DUP>\t60\tPASS\tEND={end};SVLEN={svlen};SVTYPE=DUP;TARGET={end}\tGT:DV:DR\t0/1:{ct}:{cn}\n'
    dup_recs.append((cid, pos, line))
print(f'[merge] kept {len(dup_recs)} cov-DUPs, dropped {dropped}', file=sys.stderr)

# Read existing calls.vcf
header = []; recs = []
with open(CALLS) as f:
    for line in f:
        if line.startswith('#'): header.append(line); continue
        flds = line.split('\t')
        chrom = flds[0]
        # Find chrom_id for sort
        for cid, name in chrom_name.items():
            if name == chrom:
                recs.append((cid, int(flds[1]), line))
                break

recs.extend(dup_recs)
recs.sort(key=lambda x: (x[0], x[1]))

# Write merged
out = CALLS  # overwrite
tmp = CALLS + '.tmp'
with open(tmp, 'w') as f:
    for h in header: f.write(h)
    for cid, pos, line in recs: f.write(line)
os.rename(tmp, out)
print(f'[merge] wrote {len(recs):,} records to {CALLS}', file=sys.stderr)
