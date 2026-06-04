#!/usr/bin/env python3
# Filter counts.bin in-place to a new file by keeping only sigs with max(n,t) >= min_mult.
# Preserves the schema: u64 K | K×SigCount | u64 total | u16 seqlen[K] | u8 payload[total].

import struct
import sys

if len(sys.argv) != 4:
    print('usage: filter_counts.py <in_counts.bin> <out_counts.bin> <min_mult>')
    sys.exit(2)

IN, OUT, MM = sys.argv[1], sys.argv[2], int(sys.argv[3])

with open(IN, 'rb') as f:
    K = struct.unpack('<Q', f.read(8))[0]
    sigs_raw = f.read(K * 16)
    total = struct.unpack('<Q', f.read(8))[0]
    seqlen_raw = f.read(K * 2)
    payload = f.read(total)
print(f'loaded K={K:,}  total payload={total:,}', file=sys.stderr)

seqlens = struct.unpack(f'<{K}H', seqlen_raw)
offs = [0] * (K + 1)
for i in range(K):
    offs[i+1] = offs[i] + seqlens[i]

keep = []
for i in range(K):
    sig, n, t = struct.unpack_from('<QII', sigs_raw, i * 16)
    if max(n, t) >= MM:
        keep.append(i)
K_new = len(keep)
print(f'kept: {K_new:,} ({100*K_new/K:.2f}%)  min_mult={MM}', file=sys.stderr)

new_total = sum(seqlens[i] for i in keep)
print(f'new payload: {new_total:,}', file=sys.stderr)

with open(OUT, 'wb') as fo:
    fo.write(struct.pack('<Q', K_new))
    bsigs = bytearray(K_new * 16)
    for j, i in enumerate(keep):
        bsigs[j*16:(j+1)*16] = sigs_raw[i*16:(i+1)*16]
    fo.write(bsigs)
    fo.write(struct.pack('<Q', new_total))
    bsl = bytearray(K_new * 2)
    for j, i in enumerate(keep):
        struct.pack_into('<H', bsl, j*2, seqlens[i])
    fo.write(bsl)
    # write payload by concatenating slices
    chunk = bytearray()
    for i in keep:
        chunk.extend(payload[offs[i]:offs[i+1]])
        if len(chunk) > (1<<26):
            fo.write(chunk); chunk = bytearray()
    if chunk: fo.write(chunk)
print(f'written {OUT}', file=sys.stderr)
