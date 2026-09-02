import numpy as np, struct, sys

gp, ip, dp = sys.argv[1], sys.argv[2], sys.argv[3]

with open(gp,'rb') as f:
    mg, ver, S, resv, ovs, ovc, K = struct.unpack('<4sIIIQQQ', f.read(40))
assert mg == b'SNPG' and ver == 6
print('graph.bin: S=%d K=%d' % (S, K))

plane0 = np.memmap(gp, dtype=np.uint8, mode='r', offset=40 + 8*K, shape=(K,))
sig    = np.memmap(gp, dtype=np.uint64, mode='r', offset=40,       shape=(K,))

with open(ip,'rb') as f:
    imag, iver, iS, ipad, iK, iN, iresv = struct.unpack('<4sIIIQII', f.read(32))
assert imag == b'SNPI' and iK == K
off = np.memmap(ip, dtype=np.uint32, mode='r', offset=32, shape=(K+1,))
E   = int(off[K])
w   = np.memmap(ip, dtype=np.uint64, mode='r', offset=32 + 4*(K+1), shape=(2*E,))
print('index.bin: %d kayit, N=%d' % (E, iN))

som = np.flatnonzero(plane0 == 0).astype(np.int64)
print('somatik dugum (normal duzlemi 0): %d  (%.2f%% of K)' % (som.size, 100.0*som.size/K))

drop = set()
try:
    with open(dp) as f:
        for l in f: drop.add(int(l))
except FileNotFoundError:
    pass
if drop:
    ds = np.fromiter(drop, dtype=np.uint64, count=len(drop))
    keep = ~np.isin(np.asarray(sig)[som], ds)
    som2 = som[keep]
    print('drop.sigs sonrasi somatik: %d' % som2.size)
else:
    som2 = som

for name, s in (('tum somatik', som), ('drop sonrasi', som2)):
    has = np.zeros(s.size, dtype=bool)
    B = 4_000_000
    for i in range(0, s.size, B):
        idx = s[i:i+B]
        n = (off[idx+1] - off[idx])
        first = np.asarray(w[(off[idx].astype(np.int64))*2])
        has[i:i+B] = (n > 0) & (((first >> 28) & 0x3f) == 0)
    k = int(has.sum())
    print('%-14s %d dugum -> %d tanesinde normal okuma kaydi VAR (%.2f%%)' % (name, s.size, k, 100.0*k/s.size))
