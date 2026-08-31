"""snape cagrilarini CASTLE minda_4_2 truth'una karsi skorla.

DIKKAT: bu Minda DEGIL. Minda bu kutuda kurulu degil. Ayni iki kurali uyguluyor
(min_size 50, tolerance 500) ama eslestirme kodu bizim. docs/RESULTS.md'deki
tarihsel sayilar Minda ile uretilmisti; bu betikle uretilenler onlarla bire bir
karsilastirilamaz, kendi icinde tutarlidir.

Kural: bir truth kaydi, ayni kromozom ciftinde, iki kirilma noktasi da <=tolerance
icinde olan bir cagri varsa bulunmus sayilir. INS icin yalniz POS bakilir (truth
INS kayitlarinda END yok).

kullanim: score.py <calls.vcf> <truth.vcf> [tolerance=500] [min_size=50] [--types]
"""
import sys, collections

args = [a for a in sys.argv[1:] if a != "--types"]
BYTYPE = "--types" in sys.argv
CALLS = args[0]
TRUTH = args[1]
TOL = int(args[2]) if len(args) > 2 else 500
MINSZ = int(args[3]) if len(args) > 3 else 50


def info(s):
    d = {}
    for f in s.split(';'):
        if '=' in f:
            k, v = f.split('=', 1)
            d[k] = v
    return d


def load(path, truth):
    out = []
    for line in open(path):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        i = info(f[7])
        pos = int(f[1])
        t = i.get('SVTYPE', '')
        ln = abs(int(i.get('SVLEN', 0))) if 'SVLEN' in i else 0
        if 'END' in i:
            end = int(i['END'])
        elif truth and t != 'INS':
            end = pos + ln
        else:
            end = pos
        if t != 'BND' and ln and ln < MINSZ:
            continue
        out.append([f[0], pos, end, t, ln, False])
    return out


truth = load(TRUTH, True)
calls = load(CALLS, False)
tby = collections.defaultdict(list)
for t in truth:
    tby[t[0]].append(t)

for c in calls:
    for t in tby.get(c[0], ()):
        if t[5]:
            continue
        if abs(t[1] - c[1]) <= TOL and (abs(t[2] - c[2]) <= TOL or t[3] == 'INS'):
            c[5] = True
            t[5] = True
            break

tp = sum(1 for c in calls if c[5])
fnd = sum(1 for t in truth if t[5])
P = tp / len(calls) if calls else 0.0
R = fnd / len(truth) if truth else 0.0
F = 2 * P * R / (P + R) if P + R else 0.0
name = CALLS.split('/')[-2] if '/' in CALLS else CALLS
print("%-10s cagri %5d  truth %5d  TP %5d  P %.3f  R %.3f  F1 %.3f"
      % (name, len(calls), len(truth), tp, P, R, F))

if BYTYPE:
    print("   %-6s %7s %7s %7s | %7s %8s %7s" % ("tip", "cagri", "dogru", "P", "truth", "bulunan", "R"))
    for ty in sorted(set([c[3] for c in calls]) | set([t[3] for t in truth])):
        cc = [c for c in calls if c[3] == ty]
        tt = [t for t in truth if t[3] == ty]
        a = sum(1 for c in cc if c[5])
        b = sum(1 for t in tt if t[5])
        ps = "%.3f" % (a / len(cc)) if cc else "-"
        rs = "%.3f" % (b / len(tt)) if tt else "-"
        print("   %-6s %7d %7d %7s | %7d %8d %7s" % (ty, len(cc), a, ps, len(tt), b, rs))
