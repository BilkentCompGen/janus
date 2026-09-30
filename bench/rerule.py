# rerule.py <tag> <pairs> <rule,...> (exploration of rule variants, not the caller): apply hand-written rules to the recorded measurements in calls.<tag>.vcf (INFO fields incl.
# the per-read length differences D0 normal / D1 tumor), write calls.<tag>.<rule>.vcf (FILTER PASS/LOW, duplicates within the
# benchmark tolerance marked DUPL) and print a quick record-level check against the truth; the fair harness is run separately.
import sys, re, gzip, bisect, math, collections as C
exec(open("/home/yuio/snape/bench/score.py").read().split("T, C = recs")[0])
tag = sys.argv[1]
NIST = "/home/yuio/snape/hifi/truth/nist_hg008/GRCh38_HG008-T-V0.5_somatic-stvar_PASS.draftbenchmark.vcf.gz"
CB = "/home/yuio/tmp/snape-castle/accepted3/castle-bench"
V = C.defaultdict(list)
for l in open("/home/yuio/snape/ref/vntr.trf.bed"):
    f = l.split(); V[f[0]].append((int(f[1]), int(f[2])))
for k in V: V[k].sort()
def invntr(c, p):
    v = V.get(c, []); i = bisect.bisect_right(v, (p + 50, 1 << 40))
    return any(v[j][1] + 50 >= p for j in range(max(0, i - 3), i))
def clusters(d, w=20):
    # single-linkage groups of sorted deltas within w bp; returns sizes, largest first
    if not d: return []
    d = sorted(d); g = [[d[0]]]
    for x in d[1:]:
        if x - g[-1][-1] <= w: g[-1].append(x)
        else: g.append([x])
    return sorted((len(x) for x in g), reverse=True)
def parse(p):
    R = []
    for l in open(f"{CB}/{p}/calls.{tag}.vcf"):
        if l[0] == "#": continue
        f = l.rstrip("\n").split("\t"); I = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv); m = ALT.search(f[4])
        g = lambda k: [int(x) for x in I[k].split(",")]
        dn = [] if I["D0"] == "-" else [int(x) for x in I["D0"].split(",")]; dt = [] if I["D1"] == "-" else [int(x) for x in I["D1"].split(",")]
        r = dict(line=l, f=f, key=(I["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (I["SVTYPE"], f[0], int(f[1]), f[0], int(I["END"])),
                 ty=I["SVTYPE"], L=int(I["SVLEN"]), mq=int(I["MQ"]), tdp=int(I["TDP"]), ndp=int(I["NDP"]), pc=int(I["PIECES"]), sp=g("SPAN"), sm=g("SAME"), nr=g("NEAR"), df=g("DIFF"), sq=g("SEQ"),
                 gc=int(f[-1].split(":")[1]), vn=invntr(f[0], int(f[1])), dn=dn, dt=dt, cn=clusters(dn), ct=clusters(dt), q=int(I["QPOS"]), a=int(I["SITE"].split("-")[0]), b=int(I["SITE"].split("-")[1]))
        r["same_n"], r["same_t"] = sum(abs(x) < 20 for x in dn), sum(abs(x) < 20 for x in dt)
        if not dn and not dt: r["same_t"] = r["gc"]   # no anchor pair to measure (one-sided site): the graph counts stand
        r["near_n"] = sum(20 <= abs(x) < 50 for x in dn)
        r["nmin"] = min((abs(x) for x in dn), default=10**9)
        r["nstable"] = sum(r["cn"][:2]) / max(len(dn), 1)   # share of normal reads in its two largest length clusters
        r["tstable"] = sum(r["ct"][:2]) / max(len(dt), 1)
        r["nclu"] = sum(1 for x in r["cn"] if x >= 2)         # normal length clusters with >= 2 reads
        r["vaf"] = r["same_t"] / max(len(dt), 1)
        r["vafd"] = r["same_t"] / max(r["tdp"], r["same_t"], 1)      # VAF against the anchor depth (spanning reads cannot see long-range events' reference reads)
        r["site"] = (I["PIECE"], I["SITE"])
        r["sup"] = round(1.5 + 0.1 * r["tdp"])                    # Sniffles2's coverage-scaled minimum support (local coverage)
        r["rank"] = 3 if r["ty"] in ("BND", "INV") else 2 if r["L"] >= 500 else 1  # a junction dominates the small gaps next to it
        cl = clusts(dn)
        r["nmin2"] = min((abs(m) for n_, m in cl if n_ >= max(2, 0.03 * len(dn))), default=10**9)  # nearest normal allele (>= 2 reads) to the tumor allele
        R.append(r)
    # per site: how many records, and the rank of the largest
    bysite = C.defaultdict(list)
    for r in R: bysite[r["site"]].append(r)
    for r in R:
        g = bysite[r["site"]]; r["nsite"] = len(g); r["top"] = r["rank"] == max(x["rank"] for x in g) and r["L"] == max(x["L"] for x in g if x["rank"] == r["rank"])
    return R
def clusts(d, w=20):
    if not d: return []
    d = sorted(d); g = [[d[0]]]
    for x in d[1:]:
        if x - g[-1][-1] <= w: g[-1].append(x)
        else: g.append([x])
    return [(len(x), sum(x) / len(x)) for x in g]
def base(r):  # the common core: normal depth, normal does not carry the allele, tumor carries it on >= 3 reads at VAF >= 5%
    return r["ndp"] > 3 and r["same_n"] <= 1 and r["same_n"] <= 0.03 * len(r["dn"]) and r["same_t"] >= max(3, math.ceil(0.05 * len(r["dt"])))
RULES = {
    "S1": lambda r: base(r) and r["nmin2"] >= 50,
    "S2": lambda r: base(r) and r["top"],
    "S3": lambda r: base(r) and r["nsite"] < 3,
    "S4": lambda r: base(r) and r["nmin2"] >= 50 and r["top"],
    "S5": lambda r: base(r) and r["nmin2"] >= 50 and r["top"] and r["vafd"] >= 0.05,
    "S6": lambda r: base(r) and r["nmin2"] >= 50 and r["top"] and r["same_t"] >= r["sup"],
    "S7": lambda r: base(r) and r["nmin2"] >= 50 and r["top"] and r["nsite"] < 3,
    "S8": lambda r: base(r) and r["nmin2"] >= 50 and r["top"] and r["nsite"] < 3 and r["vafd"] >= 0.05,
    "S10": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05,
    "S11": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] <= 2 or r["top"]),
    "S12": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] <= 2 or r["top"]) and r["mq"] >= 40,
    "S13": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.10 and (r["nsite"] <= 2 or r["top"]),
    "V1": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and r["dt"] and r["nsite"] == 1,
    "V2": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and r["dt"] and r["nsite"] == 1 and not r["vn"],
    "V3": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and r["dt"] and r["nsite"] == 1 and (not r["vn"] or (r["vafd"] >= 0.10 and len(r["dn"]) >= 5 and r["same_t"] >= 5)),
    "V4": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and r["dt"] and r["nsite"] == 1 and not r["vn"] and r["mq"] >= 40,
    "V5": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and r["dt"] and (r["nsite"] == 1 or r["top"]) and not r["vn"],
    "V6": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and r["dt"] and not r["vn"],
    "V7": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (not r["dt"] or r["nsite"] == 1 or r["top"]) and (not r["vn"] or (r["dt"] and r["nstable"] >= 0.8 and sum(r["ct"][:3]) >= 0.8 * len(r["dt"]) and len(r["dn"]) >= 5 and r["vafd"] >= 0.10)),
    "V8": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] == 1 or r["top"]) and (r["dt"] or r["same_t"] >= r["sup"]) and (not r["vn"] or (r["dt"] and r["nstable"] >= 0.8 and sum(r["ct"][:3]) >= 0.8 * len(r["dt"]) and len(r["dn"]) >= 5 and r["vafd"] >= 0.10)),
    "V9": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] == 1 or r["top"]) and (r["dt"] or r["same_t"] >= r["sup"]) and not r["vn"],
    "V10": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] == 1 or r["top"]) and (r["dt"] or r["same_t"] >= r["sup"]) and (not r["vn"] or (r["dt"] and r["nclu"] <= 2 and r["nstable"] >= 0.9 and r["b"] - r["a"] <= 2000 and r["vafd"] >= 0.10 and r["mq"] >= 40)),
    "V11": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] == 1 or r["top"]) and (r["dt"] or r["same_t"] >= r["sup"]) and not r["vn"] and r["mq"] >= 40 and r["vafd"] >= 0.10,
    "V12": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] == 1 or r["top"]) and (r["dt"] or r["same_t"] >= r["sup"]) and not r["vn"] and r["b"] - r["a"] <= 4000,
    "V13": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["nsite"] == 1 or r["top"]) and (r["dt"] or r["same_t"] >= r["sup"]),
    "V14": lambda r: base(r) and r["nmin2"] >= 50 and r["vafd"] >= 0.05 and (r["dt"] or r["same_t"] >= r["sup"]),
    "S9": lambda r: base(r) and r["nmin2"] >= 50 and r["top"] and r["nsite"] < 3 and r["vafd"] >= 0.10,
    "R1": base,
    "R1mq": lambda r: base(r) and r["mq"] >= 20,
    "R1n0": lambda r: base(r) and r["same_n"] == 0,
    "R1near": lambda r: base(r) and r["near_n"] == 0,
    "R1stab": lambda r: base(r) and (not r["dn"] or r["nstable"] >= 0.8),
    "R1clu": lambda r: base(r) and r["nclu"] <= 2,
    "R1vaf10": lambda r: base(r) and r["vaf"] >= 0.10,
    "R1nmin": lambda r: base(r) and (not r["dn"] or r["nmin"] >= 50),
    "R1vn": lambda r: base(r) and (not r["vn"] or (r["nstable"] >= 0.8 and r["near_n"] == 0)),
    "R1vn2": lambda r: base(r) and (not r["vn"] or (r["nstable"] >= 0.8 and r["near_n"] == 0 and r["vaf"] >= 0.10)),
    "R1all": lambda r: base(r) and r["mq"] >= 20 and r["near_n"] == 0 and (not r["dn"] or r["nstable"] >= 0.8) and r["tstable"] >= 0.8,
}
def dedupe(R, keep):
    # records that pass, within the benchmark tolerance (500 bp; INS by position), the better one (more tumor reads) stays
    ix = [i for i in range(len(R)) if keep[i]]
    ix.sort(key=lambda i: (R[i]["key"][1], R[i]["key"][2]))
    dup = [False] * len(R)
    for a in range(len(ix)):
        x = R[ix[a]]
        for b in range(a + 1, len(ix)):
            y = R[ix[b]]
            if y["key"][1] != x["key"][1] or y["key"][2] > x["key"][2] + 500: break
            if dup[ix[a]] or dup[ix[b]]: continue
            if hit(x["key"], y["key"]):
                dup[ix[b] if (x["same_t"], x["mq"]) >= (y["same_t"], y["mq"]) else ix[a]] = True
    return dup
if __name__ == "__main__":
    pairs, rules = sys.argv[2].split(","), sys.argv[3].split(",")
    for p in pairs:
        R = parse(p)
        T = recs(NIST if p == "hg008" else f"/home/yuio/snape/hifi/truth/castle_without_snape/{p}/None_minda_ensemble.vcf")
        by = {}
        for t in T: by.setdefault(t[1], []).append(t); by.setdefault(t[3], []).append(t)
        for r in R: r["tp"] = any(hit(r["key"], t) for t in by.get(r["key"][1], []))
        for name in rules:
            keep = [RULES[name](r) for r in R]; dup = dedupe(R, keep)
            P = [r for r, k, d in zip(R, keep, dup) if k and not d]
            tp = sum(r["tp"] for r in P); hit_t = sum(any(hit(t, r["key"]) for r in P) for t in T)
            pr, rc = tp / max(len(P), 1), hit_t / max(len(T), 1)
            print("%-8s %-8s PASS %5d  true %4d  P %.3f  R %.3f  F1 %.3f  (record level, quick check)" % (p, name, len(P), tp, pr, rc, 2 * pr * rc / max(pr + rc, 1e-9)))
            with open(f"{CB}/{p}/calls.{tag}.{name}.vcf", "w") as o:
                o.write("".join(l for l in open(f"{CB}/{p}/calls.{tag}.vcf") if l[0] == "#"))
                for r, k, d in zip(R, keep, dup):
                    f = list(r["f"]); f[6] = "DUPL" if k and d else "PASS" if k else "LOW"; o.write("\t".join(f) + "\n")
