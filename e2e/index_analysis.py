"""index.bin N=8 junctions icin yetiyor mu -- tahmin degil sayim.

junctions bir chunk'a hattan en fazla per_chunk=5 contig kesiyor, ve bir pencere ancak
iki yanda da min_flank=300 bp capa varsa yaziliyor. Indeks dugum basina en kucuk N=8
okumayi tutuyor. Soru: o 8 aday, 5 kabul icin yetiyor mu.

Indekste her kaydin kendi dleft/dright'i var, yani tek cekirdeklik bir dizi icin
capa testi dogrudan uygulanabiliyor. Gercek dizi birden cok cekirdek olabilir ama
lf = dleft(ilk cekirdek), rf = dright(son cekirdek) -- ikisi de kayitli alanlar.
"""
import struct, sys, numpy as np

IDX = sys.argv[1] if len(sys.argv) > 1 else "e2e/hcc1954/index.bin"
GRA = sys.argv[2] if len(sys.argv) > 2 else "e2e/hcc1954/graph.bin"
MINFLANK = 300

f = open(IDX, "rb")
mg, ver, S, pad, K, N, _ = struct.unpack("<4sIIIQII", f.read(32))
assert mg == b"SNPI", mg
print(f"index.bin: K={K:,} S={S} N={N} pad={pad}")
off = np.fromfile(f, dtype="<u4", count=K + 1)
E = int(off[K])
post = np.fromfile(f, dtype="<u8", count=2 * E).reshape(E, 2)
f.close()
print(f"  {E:,} kayit, dugum basina ortalama {E/K:.2f}")

cnt = np.diff(off)
print(f"  dolu (cnt=={N}) dugum: {int((cnt==N).sum()):,} ({100*(cnt==N).mean():.1f}%)")
for c in range(0, N + 1):
    n = int((cnt == c).sum())
    if n: print(f"    cnt={c:2d}: {n:12,} ({100*n/K:5.2f}%)")

w0 = post[:, 0]; w1 = post[:, 1]
dl = ((w1 >> 32) & 0xFFFF).astype(np.int64)
dr = ((w1 >> 48) & 0xFFFF).astype(np.int64)
ok = (dl >= MINFLANK) & (dr >= MINFLANK)
print(f"\n  capa testi (her iki yanda >= {MINFLANK} bp), tek cekirdeklik dizi varsayimiyla:")
print(f"    gecen kayit: {int(ok.sum()):,} / {E:,} = {100*ok.mean():.1f}%")

# dugum basina kac kayit capa testini geciyor
node_of = np.repeat(np.arange(K, dtype=np.int64), cnt)
passing = np.bincount(node_of[ok], minlength=K)
print("\n  dugum basina capa testini gecen kayit sayisi:")
for c in list(range(0, min(N, 8) + 1)):
    n = int((passing == c).sum())
    print(f"    {c:2d} aday: {n:12,} ({100*n/K:5.2f}%)")
ge5 = int((passing >= 5).sum())
print(f"\n  >=5 aday olan dugum: {ge5:,} ({100*ge5/K:.2f}%)")

# somatik dugumler (junctions yalniz onlari soruyor)
g = open(GRA, "rb")
mg2, ver2, S2, rv, ovs, ovc, K2 = struct.unpack("<4sIIIQQQ", g.read(40))
assert K2 == K
g.seek(40 + 8 * K)
nc = np.frombuffer(g.read(K), dtype=np.uint8)
g.close()
som = nc == 0
print(f"\n  somatik dugum: {int(som.sum()):,} ({100*som.mean():.1f}%)")
print(f"  somatik ve >=5 aday: {int((som & (passing>=5)).sum()):,} "
      f"({100*(passing[som]>=5).mean():.2f}% somatiklerin)")
print(f"  somatik ve  0 aday: {int((som & (passing==0)).sum()):,} "
      f"({100*(passing[som]==0).mean():.2f}%)")
