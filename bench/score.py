import gzip, re, sys
TOL = 500
ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")

def recs(path):
    op = gzip.open if path.endswith(".gz") else open
    out, mates = [], set()
    for ln in op(path, "rt"):
        if ln[0] == "#": continue
        f = ln.rstrip("\n").split("\t")
        if f[6] not in ("PASS", "."): continue
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        t, c, p = info.get("SVTYPE"), f[0], int(f[1])
        m = ALT.search(f[4])
        if m:
            if info.get("MATEID", info.get("MATE_ID")) in mates: continue
            mates.add(f[2])
            out.append((t, c, p, m.group(1), int(m.group(2))))
        else:
            out.append((t, c, p, c, int(info.get("END", p))))
    return out

def hit(a, b):
    for x in ((b[1], b[2], b[3], b[4]), (b[3], b[4], b[1], b[2])):
        if a[1] == x[0] and a[3] == x[2] and abs(a[2] - x[1]) <= TOL and abs(a[4] - x[3]) <= TOL: return True
    if a[0] == "INS" or b[0] == "INS":
        return a[1] == b[1] and abs(a[2] - b[2]) <= TOL
    return False

T, C = recs(sys.argv[1]), recs(sys.argv[2])
tp_t = sum(any(hit(t, c) for c in C) for t in T)
tp_c = sum(any(hit(c, t) for t in T) for c in C)
r, p = tp_t / len(T), (tp_c / len(C) if C else 0)
print("truth %d calls %d recall %.3f precision %.3f f1 %.3f" % (len(T), len(C), r, p, 2 * r * p / (r + p) if r + p else 0))
