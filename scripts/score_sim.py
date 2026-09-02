import sys, re, collections

args = [a for a in sys.argv[1:] if not a.startswith("--")]
BYTYPE = "--types" in sys.argv
CALLS, TRUTH = args[0], args[1]
TOL = int(args[2]) if len(args) > 2 else 500
MINSZ = int(args[3]) if len(args) > 3 else 50

ALT_MATE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def info(s):
    d = {}
    for f in s.split(';'):
        if '=' in f:
            k, v = f.split('=', 1)
            d[k] = v
    return d


def junctions(chrom, pos, i, alt, is_truth):
    t = i.get("SVTYPE", "")
    ln = abs(int(i["SVLEN"])) if "SVLEN" in i else 0
    end = int(i["END"]) if "END" in i else pos
    tc, tp = i.get("TARGET_CHROM"), i.get("TARGET")
    if is_truth and tc and tp:
        tp = int(tp)
        return t, ln, [((tc, tp), (chrom, pos)), ((tc, tp), (chrom, end))]
    if t == "BND" or (not is_truth and ("[" in alt or "]" in alt)):
        m = ALT_MATE.search(alt)
        if m:
            return "BND", 0, [((chrom, pos), (m.group(1), int(m.group(2))))]
        return "BND", 0, [((chrom, pos), None)]
    if t == "INS":
        return t, ln, [((chrom, pos), None)]
    return t, ln, [((chrom, pos), (chrom, end))]


def load(path, is_truth):
    out = []
    for line in open(path):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        i = info(f[7])
        t, ln, js = junctions(f[0], int(f[1]), i, f[4], is_truth)
        if t != "BND" and ln and ln < MINSZ:
            continue
        out.append([t, js, False])
    if is_truth:
        return out
    seen, uniq = set(), []
    for r in out:
        k = (r[0], tuple(tuple(sorted([e for e in j if e])) for j in r[1]))
        if k in seen:
            continue
        seen.add(k)
        uniq.append(r)
    return uniq


def near(a, b):
    return a[0] == b[0] and abs(a[1] - b[1]) <= TOL


def jmatch(x, y):
    if x[1] is None or y[1] is None:
        p = x[0] if x[1] is None else y[0]
        o = y if x[1] is None else x
        return near(p, o[0]) or (o[1] is not None and near(p, o[1]))
    return (near(x[0], y[0]) and near(x[1], y[1])) or \
           (near(x[0], y[1]) and near(x[1], y[0]))


truth = load(TRUTH, True)
calls = load(CALLS, False)
bucket = collections.defaultdict(list)
for k, t in enumerate(truth):
    for j in t[1]:
        for e in (j[0], j[1]):
            if e:
                bucket[(e[0], e[1] // 100000)].append(k)

for c in calls:
    cand = set()
    for j in c[1]:
        for e in (j[0], j[1]):
            if e:
                for d in (-1, 0, 1):
                    cand.update(bucket.get((e[0], e[1] // 100000 + d), ()))
    for k in sorted(cand):
        t = truth[k]
        if t[2]:
            continue
        if any(jmatch(a, b) for a in c[1] for b in t[1]):
            c[2] = True
            t[2] = True
            break

tp = sum(1 for c in calls if c[2])
fnd = sum(1 for t in truth if t[2])
P = tp / len(calls) if calls else 0.0
R = fnd / len(truth) if truth else 0.0
F = 2 * P * R / (P + R) if P + R else 0.0
name = CALLS.split('/')[-1]
print("%-28s cagri %5d  truth %5d  TP %5d  P %.3f  R %.3f  F1 %.3f"
      % (name, len(calls), len(truth), tp, P, R, F))

if BYTYPE:
    print("   %-6s %7s %7s %7s | %7s %8s %7s" %
          ("tip", "cagri", "dogru", "P", "truth", "bulunan", "R"))
    for ty in sorted({c[0] for c in calls} | {t[0] for t in truth}):
        cc = [c for c in calls if c[0] == ty]
        tt = [t for t in truth if t[0] == ty]
        a = sum(1 for c in cc if c[2])
        b = sum(1 for t in tt if t[2])
        print("   %-6s %7d %7d %7s | %7d %8d %7s" %
              (ty, len(cc), a, "%.3f" % (a / len(cc)) if cc else "-",
               len(tt), b, "%.3f" % (b / len(tt)) if tt else "-"))
