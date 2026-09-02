import csv, json, os, re, sys
import numpy as np

SIM = sys.argv[1] if len(sys.argv) > 1 else \
      "/home/samples/davp-samples/simulation/phylo"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/home/yuio/snape/truthsim"
TUM = os.path.join(SIM, "genomes", "tumor")
GRM = os.path.join(SIM, "genomes", "germline")


class Map:
    def __init__(self):
        self.raw = {}
        self.idx = {}

    def add(self, q, qs, qe, st, t, ts, te):
        self.raw.setdefault(q, []).append((qs, qe, 1 if st == "+" else -1, t, ts, te))

    def freeze(self):
        for c, v in self.raw.items():
            v.sort()
            tn = sorted({x[3] for x in v})
            ti = {n: i for i, n in enumerate(tn)}
            self.idx[c] = (np.array([x[0] for x in v], dtype=np.int64),
                           np.array([x[1] for x in v], dtype=np.int64),
                           np.array([x[2] for x in v], dtype=np.int8),
                           np.array([ti[x[3]] for x in v], dtype=np.int32),
                           np.array([x[4] for x in v], dtype=np.int64),
                           np.array([x[5] for x in v], dtype=np.int64), tn)
        self.raw = {}

    def lift(self, contig, pos):
        e = self.idx.get(contig)
        if e is None:
            return None
        qs, qe, st, tn, ts, te, names = e
        i = int(np.searchsorted(qs, pos, side="right")) - 1
        if i < 0 or pos >= qe[i]:
            return None
        off = pos - qs[i]
        if st[i] > 0:
            return names[tn[i]], int(ts[i] + off)
        return names[tn[i]], int(te[i] - 1 - off)


def read_paf(path, strip_query_suffix):
    m = Map()
    with open(path) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            q = f[0]
            if strip_query_suffix and q.endswith("_hapA"):
                q = q[:-5]
            m.add(q, int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]))
    m.freeze()
    return m


plan = json.load(open(os.path.join(TUM, "plan.json")))
parent, leaves = plan["parent"], plan["leaves"]


def path_to_root(n):
    p = []
    while n:
        p.append(n)
        n = parent.get(n)
    return p


node_map = {}
for n in plan["nodes"]:
    node_map[n] = read_paf(os.path.join(TUM, n, "sim.paf"), True)
sys.stderr.write("dugum haritalari: %d\n" % len(node_map))

germ = Map()
for d in sorted(os.listdir(GRM)):
    for h in ("hapA", "hapB"):
        p = os.path.join(GRM, d, "sim.%s.paf" % h)
        if os.path.exists(p):
            with open(p) as fh:
                for ln in fh:
                    f = ln.rstrip("\n").split("\t")
                    if len(f) >= 9:
                        germ.add(f[0], int(f[2]), int(f[3]), f[4], f[5],
                                 int(f[7]), int(f[8]))
germ.freeze()
sys.stderr.write("germline haritasi: %d contig\n" % len(germ.idx))

CN = re.compile(r"_c\d+$")


def lift_up(contig, pos, birth):
    anc = path_to_root(parent.get(birth)) if parent.get(birth) else []
    c, p = contig, pos
    for a in anc:
        m = node_map[a]
        if c not in m.idx:
            c2 = CN.sub("", c)
            if c2 in m.idx:
                c = c2
            else:
                return None
        r = m.lift(c, p)
        if r is None:
            return None
        c, p = r
    c = CN.sub("", c)
    return germ.lift(c, p)


os.makedirs(OUT, exist_ok=True)
muts = {}
with open(os.path.join(TUM, "mutations.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["class"] == "SV":
            muts.setdefault(r["birth_node"], []).append(r)

TGT = re.compile(r"TARGET=(\d+)")
TCH = re.compile(r"TARGET_CHROM=([^;]+)")
stat = []
for leaf in leaves:
    rows, miss, n = [], 0, 0
    for node in reversed(path_to_root(leaf)):
        for r in muts.get(node, []):
            n += 1
            a = lift_up(r["contig"], int(r["pos_parent"]) - 1, node)
            b = lift_up(r["contig"], int(r["end_parent"]) - 1, node)
            if a is None or b is None or a[0] != b[0]:
                miss += 1
                continue
            info = ["SVTYPE=%s" % r["svtype"], "SVLEN=%s" % r["svlen"],
                    "END=%d" % (b[1] + 1), "BIRTH=%s" % node,
                    "SRCCHROM=%s" % r["contig"]]
            mt, mc = TGT.search(r["info"]), TCH.search(r["info"])
            if mt and mc:
                t = lift_up(mc.group(1), int(mt.group(1)) - 1, node)
                if t is None:
                    miss += 1
                    continue
                info += ["TARGET=%d" % (t[1] + 1), "TARGET_CHROM=%s" % t[0]]
            rows.append((a[0], a[1] + 1, r["mut_id"], "N", "<%s>" % r["svtype"],
                         ";".join(info)))
    rows.sort(key=lambda x: (x[0], x[1]))
    p = os.path.join(OUT, "truth.%s.vcf" % leaf)
    with open(p, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">\n")
        fh.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"\">\n")
        fh.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"\">\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, pos, i, ref, alt, info in rows:
            fh.write("%s\t%d\t%s\t%s\t%s\t100\tPASS\t%s\n" % (c, pos, i, ref, alt, info))
    stat.append((leaf, n, len(rows), miss))
    sys.stderr.write("%s: yol uzeri %d SV, kaldirilan %d, kaldirilamayan %d\n"
                     % (leaf, n, len(rows), miss))

with open(os.path.join(OUT, "summary.tsv"), "w") as fh:
    fh.write("leaf\tpath_sv\tlifted\tunlifted\n")
    for s in stat:
        fh.write("%s\t%d\t%d\t%d\n" % s)
