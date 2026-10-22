import sys, os, re, json, csv, time, glob, collections, pickle, subprocess, statistics as st
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(H)
sys.path.insert(0, H)
sys.dont_write_bytecode = True
import paper
R = f"{H}/runs"
RC = f"{H}/runs_clean"
OUT = os.environ.get("FIGOUT", f"{ROOT}/paper/overleaf/figures")
TAB = os.environ.get("TABOUT", f"{H}/out")
os.makedirs(OUT, exist_ok=True); os.makedirs(TAB, exist_ok=True)
csv.field_size_limit(10**9)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

TEX = f"{ROOT}/paper/figures/system/style.tex"
_s = open(TEX).read()
_pal = re.search(r"\\providecommand\\palette\{(\w+)\}", _s).group(1)
_th = re.search(r"\\newcommand\\theme%s\s*\{\\settheme((?:\{[0-9A-Fa-f]{6}\}){6})\}" % _pal, _s).group(1)
_v = re.findall(r"\{([0-9A-Fa-f]{6})\}", _th)
C = dict(zip(["normal", "tumor", "ink", "ink2", "muted", "shared"], ["#" + h for h in _v]))
C.update({n: "#" + h for n, h in re.findall(r"\\definecolor\{(\w+)\}\{HTML\}\{([0-9A-Fa-f]{6})\}", _s)})
GEO = {n: float(v) for n, v in re.findall(r"\\newcommand\{\\(plot\w+)\}\{([\d.]+)\}", _s)}
MM = 1 / 25.4
FONT = next((f for f in ("Helvetica", "Arial", "Liberation Sans", "Nimbus Sans", "DejaVu Sans")
             if any(f == x.name for x in font_manager.fontManager.ttflist)), "DejaVu Sans")
plt.rcParams.update({
    "font.family": FONT, "font.size": GEO["plotticks"], "axes.labelsize": GEO["plotlabels"], "axes.titlesize": GEO["plotlabels"],
    "xtick.labelsize": GEO["plotticks"], "ytick.labelsize": GEO["plotticks"], "legend.fontsize": GEO["plotticks"],
    "axes.edgecolor": C["ink2"], "axes.labelcolor": C["ink"], "xtick.color": C["ink2"], "ytick.color": C["ink2"],
    "text.color": C["ink"], "axes.linewidth": GEO["plotaxislw"], "xtick.major.width": GEO["plotaxislw"],
    "ytick.major.width": GEO["plotaxislw"], "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": GEO["plotlw"],
    "lines.markersize": GEO["plotms"], "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 300,
})

from matplotlib.colors import LinearSegmentedColormap
PAIRNAME = {"hcc1954": "HCC1954", "h2009": "NCI-H2009", "hcc1937": "HCC1937", "h1437": "NCI-H1437", "hg008": "HG008",
            "colo829": "COLO829", "hcc1395": "HCC1395", "hg009": "HG009", "h209": "NCI-H209", "h2126": "NCI-H2126"}

JANUS = C["tumor"]
GROUP = [("Packing", ["pack_n", "pack_t"]), ("Core extraction", ["extract_n", "extract_t"]), ("Counting", ["count_n", "count_t"]),
         ("Graph and pruning", ["graph", "prune"]), ("Selection", ["select"]), ("Placement", ["place", "lift", "clip", "patch"]),
         ("minimap2", ["local", "gaps"]), ("Calling", ["call"])]
STAGE = dict(zip([g for g, _ in GROUP], ["#d8cfc0", "#b7c9dd", "#86a8cc", "#4f7cae", "#f0b9a8", "#d9735f", "#8c8c8c", "#5fa37a"]))
STEPCOL = {s: STAGE[g] for g, ss in GROUP for s in ss}
CALLERCOL = {"janus": JANUS, "Severus": "#3f6ea5", "SAVANA": "#4bb3b3", "nanomonsv": "#5fa37a", "Sniffles2": "#8c8c8c", "dysgu": "#7e6bb0",
             "SVision-pro": "#e0a030", "cuteSV": "#b98b5e", "DeBreak": "#c9a0c4", "kled": "#b0b0b0", "SVDSS": "#d0d0d0"}
P5 = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]

W = 160
FS = 7.0
plt.rcParams.update({"font.size": 6.0, "axes.labelsize": 7.0, "axes.titlesize": 7.0, "xtick.labelsize": 6.0, "ytick.labelsize": 6.0, "legend.fontsize": 6.3,
                     "axes.linewidth": 0.45, "lines.linewidth": 0.9, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
GRID = dict(color="#ebebeb", lw=0.4, zorder=0)
NEUTRAL = "#8c8c8c"
LIGHT = "#c9c9c9"
TECHCOL = {"HiFi": "#2f5d8a", "ONT": "#1f9e89", "Illumina": "#d0679d"}
PAIRCOL = dict(zip(P5 + ["hg008", "colo829"], ["#6a3d9a", "#e0a21a", "#8c5a2b", "#b39ddb", "#56616c", "#cfae7a", "#9a9a32"]))
MARK = {"janus": "o", "Severus": "s", "SAVANA": "^", "nanomonsv": "D", "dysgu": "P", "Sniffles2": "v", "SVision-pro": "X", "cuteSV": "<", "DeBreak": ">", "kled": "h", "SVDSS": "d"}
SVTYPE = {"DEL": "#3f6ea5", "INS": "#e0a030", "DUP": "#5fa37a", "INV": "#7e6bb0", "BND": "#8f5b34"}
FILTER = {"PASS": JANUS, "LOW": NEUTRAL, "DUPL": LIGHT}
MATCH, NOMATCH = "#3f6ea5", LIGHT
HEAT = LinearSegmentedColormap.from_list("f", ["#f7f7f7", "#a9c8e2", "#2f5d8a"])
DIV = LinearSegmentedColormap.from_list("d", ["#2f5d8a", "#f7f7f7", JANUS])


def tag(ax, s, x, dy):
    w, h = ax.figure.get_size_inches() / MM
    ax.figure.text(x / w, ax.get_position().y1 + dy / h, s, fontweight="bold", fontsize=9.5, va="baseline", ha="left")


def lum(c):
    r, g, b = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in matplotlib.colors.to_rgb(c)]
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def ink(c, cmap):
    L = [lum(cmap(v)) for v in np.linspace(0, 1, 101)]
    return "white" if lum(c) < min(L) + (max(L) - min(L)) / 3 else C["ink"]


def ffmt(v, _=None):
    s = "0" if abs(v) < 1e-9 else "1" if abs(v - 1) < 1e-9 else f"{abs(v):.3f}".rstrip("0").lstrip("0")
    return ("−" if v < -1e-9 else "") + s


def fticks(ax, axis="y", ticks=(0, 0.25, 0.5, 0.75, 1)):
    a = ax.yaxis if axis == "y" else ax.xaxis
    if ticks is not None: a.set_ticks(list(ticks))
    a.set_major_formatter(plt.FuncFormatter(ffmt))


def page(w, h):
    fig = plt.figure(figsize=(w * MM, h * MM))
    return fig, lambda x, y, pw, ph: fig.add_axes([x / w, y / h, pw / w, ph / h])


def title(ax, s, y=1.03, fs=0.8):
    ax.text(0, y, s, transform=ax.transAxes, fontsize=FS + fs, va="bottom", ha="left")
TECH = TECHCOL
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.legend_handler import HandlerTuple

ONT = ["Severus ONT", "SAVANA ONT", "nanomonsv ONT", "Sniffles2 ONT"]
ILL = ["SvABA", "GRIPSS", "Manta"]
EXT = {"snape": "janus", "severus_phased": "Severus", "savana_phased": "SAVANA", "nanomonsv": "nanomonsv", "dysgu": "dysgu", "svision-pro": "SVision-pro",
       "sniffles2": "Sniffles2", "cutesv": "cuteSV", "debreak": "DeBreak", "kled": "kled", "svdss": "SVDSS"}
COV = ["Severus", "SAVANA", "nanomonsv"]
TOOLKEY = {"janus": "snape", "Severus": "severus", "SAVANA": "savana", "nanomonsv": "nanomonsv"}
DEP, DL = ["full", "30x", "20x", "15x", "10x", "5x"], ["full", "30×", "20×", "15×", "10×", "5×"]
PUR, PL = ["100p", "80p", "60p", "40p", "20p"], ["100%", "80%", "60%", "40%", "20%"]
TYPES = ["BND", "DEL", "DUP", "INS", "INV"]
SIZES = [("BND", "no length"), ("50-100", "50–100 bp"), ("100-1k", "0.1–1 kb"), ("1k-10k", "1–10 kb"), ("10k-100k", "10–100 kb"), (">100k", ">100 kb")]
CAL5 = ["janus", "Severus", "SAVANA", "nanomonsv", "Sniffles2"]
LBL = {"colo829": {"janus": (0.86, 0.96, "center"), "SAVANA": (0.66, 0.86, "center"), "Severus": (0.46, 0.84, "right"), "nanomonsv": (0.78, 0.71, "left"),
                   "dysgu": (0.775, 0.635, "left"), "SVision-pro": (0.17, 0.70, "left"), "DeBreak": (0.17, 0.635, "left"), "Sniffles2": (0.17, 0.57, "left"),
                   "kled": (0.17, 0.50, "left"), "cuteSV": (0.17, 0.43, "left"), "SVDSS": (0.075, 0.226, "left")},
       "hg008": {"janus": (0.99, 0.96, "right"), "Severus": (0.74, 0.905, "right"), "nanomonsv": (0.84, 0.76, "right"), "SAVANA": (0.99, 0.545, "right"),
                 "dysgu": (0.815, 0.69, "right"), "Sniffles2": (0.17, 0.83, "left"), "DeBreak": (0.16, 0.655, "left"), "kled": (0.315, 0.59, "left"),
                 "SVision-pro": (0.29, 0.475, "left"), "cuteSV": (0.13, 0.405, "left"), "SVDSS": (0.15, 0.337, "left")}}
ROWS = P5 + ["hg008", "colo829"]


def save(fig, key):
    fig.savefig(f"{OUT}/fig_{key}.pdf", facecolor="white", dpi=200)
    plt.close(fig)


SRC = os.environ.get("SRCDATA", f"{os.path.dirname(OUT)}/sourcedata")


def src(name, header, rows):
    os.makedirs(SRC, exist_ok=True)
    with open(f"{SRC}/{name}.csv", "w", newline="") as o:
        w = csv.writer(o); w.writerow(header)
        w.writerows([round(float(v), 10) if isinstance(v, (float, np.floating)) else v for v in r] for r in rows)


def prf(tp, fp, fn):
    return tp / max(1, tp + fp), tp / max(1, tp + fn), 2 * tp / max(1, 2 * tp + fp + fn)


def dotrows(pn, s, order):
    return [[pn, n, tech(n), PAIRNAME.get(p, p)] + [D[s][p][n].get(k, "") for k in ("P", "R", "F1", "tp", "fp", "fn", "truth")] for n in order for p in P5] + \
           [[pn, n, tech(n), "mean of five pairs"] + [mean(s, n, k) for k in ("P", "R", "F1")] + [""] * 4 for n in order]


def poolrows(pn, s, kind, cls, labs, cal):
    out = []
    for n in cal:
        for c, l in zip(cls, labs):
            tp, fp, fn = [sum(D[s][p][n][kind][c][x] for p in P5) for x in ("TP", "FP", "FN")]
            out.append([pn, n, tech(n), l, *prf(tp, fp, fn), tp, fp, fn, tp + fn])
    return out


ACCH = ["panel", "call_set", "technology", "item", "precision", "recall", "f1", "tp", "fp", "fn", "consensus_svs"]


E = f"{H}/eval_sets"
OURS = ["Severus", "SAVANA", "nanomonsv", "Sniffles2", "dysgu", "SVision-pro", "cuteSV", "DeBreak", "kled", "SVDSS"]
REST = ["Severus ONT", "SAVANA ONT", "nanomonsv ONT", "Sniffles2 ONT", "SvABA", "GRIPSS", "Manta"]
SETNAMES = {"all_with": OURS + ["janus"] + REST, "all_without": OURS + REST}
SIZEKEYS = ["BND", "50-100", "100-1k", "1k-10k", "10k-100k", ">100k"]


def sizeclass(t, n):
    try:
        n = abs(float(n))
    except ValueError:
        n = 0
    if t == "BND" or n <= 0:
        return "BND"
    return "50-100" if n <= 100 else "100-1k" if n <= 1000 else "1k-10k" if n <= 10000 else "10k-100k" if n <= 100000 else ">100k"


def breakdown(d, col, truth):
    ty = {t: dict(TP=0, FP=0, FN=0) for t in TYPES}
    sz = {s: dict(TP=0, FP=0, FN=0) for s in SIZEKEYS}
    rd = csv.reader(open(f"{d}/None_support.tsv"), delimiter="\t")
    h = next(rd)
    c, e = h.index("Minda_IDs") + 1 + col, h.index(truth) if truth else h.index("ensemble")
    for r in rd:
        k = {("True", "True"): "TP", ("True", "False"): "FN", ("False", "True"): "FP"}.get((r[e], r[c]))
        if k:
            t = r[h.index("SVTYPE")]
            st_ = r[h.index("STRANDS")] if "STRANDS" in h else (re.search(r"STRANDS=([+-]{2})", r[h.index("INFO_x")]) or [None, ""])[1]
            inv = t == "INV" or (r[0] == r[4] and st_ in ("++", "--"))
            if inv or t in ty:
                ty["INV" if inv else t][k] += 1
            sz[sizeclass(t, r[h.index("SVLEN")])][k] += 1
    if truth:
        for r in csv.DictReader(open(f"{d}/None_Unknown source_fp.tsv"), delimiter="\t"):
            t = r["SVTYPE_x"]
            inv = t == "INV" or (t == "BND" and r["#CHROM_x"] == r["#CHROM_y"] and re.search(r"^[A-Za-z.]+\]|\[[A-Za-z.]+$", r["ALT_x"]))
            if inv or t in ty:
                ty["INV" if inv else t]["FP"] += 1
            sz[sizeclass(t, r["SVLEN"])]["FP"] += 1
    return ty, sz


def overall(d):
    return [r for r in csv.reader(open(f"{d}/results/overall.tsv"), delimiter="\t")][1:]


def record(r, d, col, truth=None):
    ty, sz = breakdown(d, col, truth)
    return dict(tp=int(r[1]), fn=int(r[2]), fp=int(r[3]), P=float(r[4]), R=float(r[5]), F1=float(r[6]), calls=int(r[7]), truth=int(r[8]), type=ty, size=sz)


def accuracy():
    out, extra = {}, {}
    for s, names in SETNAMES.items():
        out[s] = {}
        for p in P5:
            d = f"{E}/{p}/{s}"; rows = overall(d)
            assert len(rows) == len(names), (s, p, len(rows))
            out[s][p] = {n: record(r, d, i) for i, (n, r) in enumerate(zip(names, rows))}
            if s == "all_without":
                out[s][p]["janus"] = record(overall(f"{d}_snape")[0], f"{d}_snape", 1, "MindaV0.0.2")
    for s in SETNAMES:
        extra[s] = {}
        for p in P5:
            extra[s][p] = {}
            for n, t in {"Sniffles2 (mosaic)": "sniffles2_mosaic", "Delly": "delly"}.items():
                if os.path.exists(f"{E}/{p}/{s}_{t}/results/overall.tsv"):
                    r = overall(f"{E}/{p}/{s}_{t}")[0]
                    extra[s][p][n] = dict(P=float(r[4]), R=float(r[5]), F1=float(r[6]), calls=int(r[7]), truth=int(r[8]))
    return json.loads(json.dumps(dict(sets=out, extra=extra)))


SIGTOOLS = ["snape", "severus", "savana", "nanomonsv"]
ENS = {"snape": "K_Unknown source", "severus": "A_Severus_v1.7", "savana": "B_SAVANAv1.3.8", "nanomonsv": "C_nanomonsv-0.9.0"}
B = 2000


def sigunit(d, name):
    for x in csv.DictReader(open(f"{d}/results/overall.tsv"), delimiter="\t"):
        if x["Caller"] == name:
            tp, fp, n = int(x["True Positives"]), int(x["False Positives"]), int(x["Ensemble Records"])
    with open(f"{d}/None_{name}_fn.tsv") as h:
        rd = csv.reader(h, delimiter="\t"); next(rd)
        miss = {tuple(x[:3]) for x in rd if x}
    return dict(tp=tp, fp=fp, n=n, miss=miss)


def sigunits(p, k):
    if k == "full":
        d = f"{E}/{p}/all_with"
        return {t: sigunit(d, ENS[t]) for t in SIGTOOLS if os.path.exists(f"{d}/None_{ENS[t]}_fn.tsv")}
    u = {}
    for t in SIGTOOLS:
        d = f"{R}/coverage_eval/{p}_{k}/{t}/minda"
        f = [x for x in glob.glob(f"{d}/None_*_fn.tsv") if "MindaV" not in x]
        if f:
            u[t] = sigunit(d, os.path.basename(f[0])[5:-7])
    return u


def draws(u, rng):
    keys = sorted(set().union(*(v["miss"] for v in u.values())))
    n = next(iter(u.values()))["n"]
    M = np.array([[k in u[t]["miss"] for k in keys] for t in u], bool).reshape(len(u), len(keys))
    idx = rng.integers(0, n, size=(B, n))
    hit, obs = {}, {}
    for i, t in enumerate(u):
        miss = np.zeros(n, bool); miss[:len(keys)] = M[i]
        R_ = 1 - miss[idx].mean(1)
        c = u[t]["tp"] + u[t]["fp"]
        P_ = rng.binomial(c, u[t]["tp"] / max(c, 1), size=B) / max(c, 1)
        hit[t] = np.where(P_ + R_ > 0, 2 * P_ * R_ / np.maximum(P_ + R_, 1e-12), 0)
    for t in u:
        c = u[t]["tp"] + u[t]["fp"]; P_ = u[t]["tp"] / max(c, 1); R_ = 1 - len(u[t]["miss"]) / n
        obs[t] = 2 * P_ * R_ / max(P_ + R_, 1e-12)
    return hit, obs


def pval(d):
    return float(min(1.0, 2 * min((d <= 0).mean(), (d >= 0).mean())) if len(d) else 1.0)


def significance():
    rng = np.random.default_rng(0); out = {}
    for k in ["full", "30x", "20x", "15x", "10x", "5x", "100p", "80p", "60p", "40p", "20p"]:
        per = {}
        for p in P5:
            u = sigunits(p, k)
            if "snape" in u and len(u) == 4:
                per[p] = draws(u, rng)
        if not per:
            continue
        pairs = sorted(per)
        mean_obs = {t: float(np.mean([per[p][1][t] for p in pairs])) for t in SIGTOOLS}
        best = max((t for t in SIGTOOLS if t != "snape"), key=lambda t: mean_obs[t])
        diff = np.mean([per[p][0]["snape"] - per[p][0][best] for p in pairs], 0)
        snape_mean = np.mean([per[p][0]["snape"] for p in pairs], 0)
        pp = {}
        for p in pairs:
            bp = max((t for t in SIGTOOLS if t != "snape"), key=lambda t: per[p][1][t])
            d = per[p][0]["snape"] - per[p][0][bp]
            pp[p] = dict(best=bp, diff=float(per[p][1]["snape"] - per[p][1][bp]), lo=float(np.percentile(d, 2.5)), hi=float(np.percentile(d, 97.5)))
        out[k] = dict(pairs=pairs, mean=mean_obs, best=best, diff=float(mean_obs["snape"] - mean_obs[best]), lo=float(np.percentile(diff, 2.5)), hi=float(np.percentile(diff, 97.5)),
                      snape_lo=float(np.percentile(snape_mean, 2.5)), snape_hi=float(np.percentile(snape_mean, 97.5)), per_pair=pp, p=pval(diff))
    ext = {}
    for p in ("hg008", "colo829"):
        u = {}
        for t, d in {"snape": "snape", "severus": "severus_phased", "savana": "savana_phased", "nanomonsv": "nanomonsv"}.items():
            d = f"{R}/external_eval/{p}/{d}/minda"
            rs = list(csv.DictReader(open(f"{d}/results/overall.tsv"), delimiter="\t"))
            u[t] = sigunit(d, [x for x in rs if not (x["False Positives"] == "0" and x["True Positives"] == x["Ensemble Records"])][0]["Caller"])
        hit, obs = draws(u, rng)
        bp = max((t for t in SIGTOOLS if t != "snape"), key=lambda t: obs[t])
        d = hit["snape"] - hit[bp]
        ext[p] = dict(f1=obs, best=bp, diff=float(obs["snape"] - obs[bp]), lo=float(np.percentile(d, 2.5)), hi=float(np.percentile(d, 97.5)), p=pval(d))
    return json.loads(json.dumps(dict(external=ext, settings=out)))


def steps(d, p):
    return {x["step"]: x for x in json.load(open(f"{d}/full/{p}/run.json"))["steps"]}


FASTQPACK = {"hcc1954": 235, "h2009": 284, "hcc1937": 246, "h1437": 236, "hcc1395": 140}


def cost():
    out = {"time": {}, "mem": {}}
    for p in P5:
        s = {k: v for k, v in steps(R, p).items() if k.startswith("pack")}
        s.update({k: v for k, v in steps(RC, p).items() if not k.startswith("pack")})
        out["time"][p] = {g: sum(s[k]["wall_s"] for k in ks) for g, ks in GROUP}
        out["time"][p]["Packing"] = FASTQPACK[p]
        out["mem"][p] = {g: max(s[k]["peak_rss_kib"] for k in ks) / 2**20 for g, ks in GROUP}
    return out


def simdata():
    S = json.load(open(f"{R}/sim/scores.json"))
    for k in S:
        d = f"{RC}/sim/timing/{k.replace('-', '_')}"
        t = [open(f"{d}/pair/T{i:02d}/time.txt").read().split()[:2] for i in range(1, 21)]
        j = open(f"{d}/joint/time.txt").read().split()[:2]
        S[k]["time"] = dict(pair_wall=[float(x[0]) for x in t], joint_wall=float(j[0]), pair_peak_kb=[float(x[1]) for x in t], joint_peak_kb=float(j[1]))
    return S


def mm2():
    out = {"5-5": dict(wall_min=117.6, peak_gb=19.0)}
    for k in ("5-5", "20-20"):
        f = sorted(glob.glob(f"{R}/sim/mm2/{k.replace('-', '_')}/*.time"))
        if len(f) == 21:
            v = [open(x).read().split() for x in f]
            out[k] = dict(wall_min=sum(float(x[0]) for x in v) / 60, peak_gb=max(int(x[1]) for x in v) / 2**20)
    return out


SIMSETS = ["20-20", "20-10", "20-5", "10-10", "5-5"]
SIMT = ["T%02d" % i for i in range(1, 21)]
ALL20 = frozenset(range(20))


def carriers(d):
    pr, du, own, NB, SM, SP, cut = [d[k] for k in ("pr", "du", "own", "NB", "SM", "SP", "cut")]
    N = len(pr)
    cut = np.broadcast_to(np.asarray(cut, float), (N,))
    conf = np.array([len({own[j] for j in NB.get(i, ()) if own[j] != own[i] and pr[j] >= cut[j]}) for i in range(N)])
    nown = np.array([len({own[j] for j in NB.get(i, ())} | {own[i]}) for i in range(N)])
    anyown = [set() for _ in range(20)]
    for i in np.where(~du)[0]:
        anyown[own[i]].update(NB[i])
    acc = ~du & ((pr >= cut + 0.2) | ((pr >= 0.1) & (conf >= 1)))
    root = list(range(N))

    def find(i):
        while root[i] != i:
            root[i] = root[root[i]]
            i = root[i]
        return i

    for i in range(N):
        for j in NB.get(i, ()):
            root[find(i)] = find(j)
    carr = collections.defaultdict(set)
    for i in np.where(acc)[0]:
        carr[find(i)].add(int(own[i]))
    for k in range(20):
        mine = np.where(acc & (own == k))[0]
        have = set(mine)
        for i in mine:
            have.update(NB[i])
        seen = set()
        ok = acc & (own != k) & (SM[:, k] >= 1) & (0.25 * SP[:, k] <= SM[:, k]) & (nown >= 3)
        for i in np.where(ok)[0]:
            if i in anyown[k] or i in have or i in seen:
                continue
            carr[find(i)].add(k)
            seen.update(NB[i])
    return collections.Counter(frozenset(v) for v in carr.values())


def splitkey(s):
    return min(s, ALL20 - s, key=lambda x: (len(x), sorted(x)))


def fits(a, b):
    return not (a & b) or a <= b or b <= a or not ((ALL20 - a) & (ALL20 - b))


def compatible(cnt):
    sup = collections.Counter()
    for s, n in cnt.items():
        if 1 < len(s) < 19:
            sup[splitkey(s)] += n
    chosen = []
    for s, n in sorted(sup.items(), key=lambda x: (-x[1], sorted(x[0]))):
        if all(fits(s, c) for c in chosen):
            chosen.append(s)
    return chosen


def nwksplits(nwk):
    s, stk, out, i = re.sub(r":[-0-9.eE]+", "", nwk.strip().rstrip(";")), [], set(), 0
    while i < len(s):
        if s[i] == "(":
            stk.append(set())
        elif s[i] == ")":
            c = stk.pop()
            if 1 < len(c) < 19:
                out.add(splitkey(frozenset(c)))
            if stk:
                stk[-1] |= c
        elif s[i] == "T":
            stk[-1].add(int(s[i + 1:i + 3]) - 1)
            i += 2
        i += 1
    return out


NJR = 'library(ape); for (f in commandArgs(trailingOnly = TRUE)) { d <- read.table(f, colClasses = c("character", "integer")); m <- do.call(cbind, strsplit(rep(d$V1, d$V2), "")) == "1"; rownames(m) <- sprintf("T%02d", 1:20); write.tree(nj(dist(m, method = "manhattan")), sub("tsv$", "nj.nwk", f)) }'


def simtrees():
    PAR = json.load(open(f"{R}/sim/truth/tree.json"))["parent"]
    below = collections.defaultdict(set)
    for t in SIMT:
        n = t
        while n:
            below[n].add(SIMT.index(t))
            n = PAR[n]
    TRUE = {frozenset(v): k for k, v in below.items()}
    SPLITS = {s for s in TRUE if 1 < len(s) < 20}
    ROOT = splitkey(frozenset(below["n18"])); TS = {splitkey(s) for s in SPLITS}
    W = f"{TAB}/nj"; os.makedirs(W, exist_ok=True)

    def rooted(sp):
        c = set()
        if ROOT in sp:
            for s in sp:
                for side in (s, ALL20 - s):
                    if side <= ROOT or side <= ALL20 - ROOT:
                        c.add(side)
        return c

    tr, cnts = {}, {}
    for name in SIMSETS:
        cnt = carriers(pickle.load(open(f"{R}/sim/settings/{name.replace('-', '_')}/joint/carry.pkl", "rb"))); cnts[name] = cnt
        rec = compatible(cnt); tot = sum(cnt.values()); cl = rooted(set(rec))
        tr[name] = dict(events=tot, sets=len(cnt), exact=sum(n for s, n in cnt.items() if s in TRUE), shared=len(set(rec) & TS), rf=len(set(rec) ^ TS),
                        clades=[sorted(c) for c in cl], wrong=[sorted(c) for c in cl if c not in SPLITS])
        with open(f"{W}/{name}.tsv", "w") as o:
            for s, n in sorted(cnt.items(), key=lambda x: sorted(x[0])):
                o.write("".join("1" if i in s else "0" for i in range(20)) + f"\t{n}\n")
    subprocess.run(["docker", "run", "--rm", "--log-driver", "none", "-v", f"{W}:/w", "-u", f"{os.getuid()}:{os.getgid()}", "--entrypoint", "/usr/local/bin/_entrypoint.sh",
                    "janus-bench/phangorn", "Rscript", "-e", NJR] + [f"/w/{n}.tsv" for n in SIMSETS], check=True, stdout=subprocess.DEVNULL)
    nj = {}
    for name in SIMSETS:
        sp = nwksplits(open(f"{W}/{name}.nj.nwk").read()); cl = rooted(sp)
        nj[name] = dict(splits=len(sp), shared=len(sp & TS), true=len(TS), rf=len(sp ^ TS), clades=[sorted(c) for c in cl], wrong=[sorted(c) for c in cl - SPLITS], missing=[TRUE[c] for c in SPLITS - cl])
    return json.loads(json.dumps(dict(parent=PAR, true=[sorted(s) for s in SPLITS], sets=tr))), json.loads(json.dumps(nj))


ACC = accuracy()
D = ACC["sets"]
SIG = significance()

def base(n):
    return n.replace(" ONT", "")


def tech(n):
    return "ONT" if n in ONT else "Illumina" if n in ILL else "HiFi"


def color(n):
    return JANUS if n == "janus" else CALLERCOL.get(base(n), TECHCOL[tech(n)])


def label(n):
    return n.replace(" ONT", " (ONT)")


def mean(s, n, k):
    return st.mean(D[s][p][n][k] for p in P5)


def cov():
    S = {}
    for l in open(f"{R}/coverage_eval/scores.tsv"):
        s, t, n, pp, r, f = l.split()
        p, k = s.rsplit("_", 1)
        S[(p, k, t)] = (float(pp), float(r), float(f))
    return S


def external():
    E = collections.defaultdict(dict)
    for l in open(f"{R}/external_eval/scores.tsv"):
        f = l.split("\t")
        if f[1] in EXT: E[f[0]][EXT[f[1]]] = (float(f[5]), float(f[6]), float(f[7]))
    return E


def pooled(s, n, kind, c):
    tp, fp, fn = [sum(D[s][p][n][kind][c][x] for p in P5) for x in ("TP", "FP", "FN")]
    return 2 * tp / max(1, 2 * tp + fp + fn), tp + fn


def f1_dots(ax, s, order, ttl, show_y, techs=False):
    for i, n in enumerate(order):
        v = [D[s][p][n]["F1"] for p in P5]; m = st.mean(v)
        ax.plot([min(v), max(v)], [i, i], color="#d0d0d0", lw=1.2, zorder=1)
        for p, x in zip(P5, v): ax.scatter([x], [i], s=7, color=PAIRCOL[p], zorder=2, lw=0)
        ax.scatter([m], [i], s=44, marker="|", color="black", lw=0.9, zorder=1.5)
        ax.text(1.03, i, f"{m:.3f}", fontsize=FS, va="center", ha="left", fontweight="bold" if m >= max(mean(s, x, "F1") for x in order) - 1e-9 else "normal", transform=ax.get_yaxis_transform())
        ax.text(1.47, i, f"{mean(s, n, 'P'):.2f}", fontsize=FS - 0.5, va="center", ha="right", color="#555555", transform=ax.get_yaxis_transform())
        ax.text(1.69, i, f"{mean(s, n, 'R'):.2f}", fontsize=FS - 0.5, va="center", ha="right", color="#555555", transform=ax.get_yaxis_transform())
    for x, hd, al in ((1.03, "F1", "left"), (1.47, "P", "right"), (1.69, "R", "right")): ax.text(x, -1.05, hd, fontsize=FS - 0.3, va="center", ha=al, color="#555555", transform=ax.get_yaxis_transform())
    ax.set_ylim(len(order) - 0.4, -1.6); ax.set_yticks(range(len(order)))
    if show_y:
        ax.set_yticklabels([label(n) for n in order], fontsize=FS)
        mixed = techs or any(tech(n) != "HiFi" for n in order)
        for t, n in zip(ax.get_yticklabels(), order): t.set_color(JANUS if n == "janus" else TECHCOL[tech(n)] if mixed else "#333333"); t.set_fontweight("bold" if n == "janus" else "normal")
    else:
        ax.set_yticklabels([])
    ax.set_xlim(0, 1); fticks(ax, "x")
    ax.set_xlabel("F1 score"); title(ax, ttl, 1.01, 0.0)
    ax.grid(axis="x", **GRID); ax.tick_params(axis="y", length=0)


def ext_scatter(ax, E, sample, ttl):
    for f in (0.2, 0.4, 0.6, 0.8):
        x = np.linspace(f / (2 - f) + 1e-3, 1.0, 300); y = f * x / (2 * x - f); m = y <= 1.0
        ax.plot(x[m], y[m], color="#d8d8d8", lw=0.5, zorder=0)
        if f < 0.7: ax.text(0.995, f / (2 - f) - 0.012, f"F1 {f:.1f}", fontsize=FS - 1.4, color="#9a9a9a", ha="right", va="top")
    for n, (p, r, f1) in sorted(E.items(), key=lambda x: x[1][2]):
        c = color(n); top = n == "janus"
        ax.scatter([p], [r], s=46 if top else 24, marker=MARK.get(n, "o"), facecolor=c, edgecolor="#333333" if c in ("#d0d0d0", "#b0b0b0") else c, lw=0.5, zorder=5 if top else 4)
        lx, ly, ha = LBL[sample][n]
        ax.annotate(f"{n} {f1:.2f}", (p, r), xytext=(lx, ly), fontsize=FS - 0.5, ha=ha, va="center", color=JANUS if top else "#333333", fontweight="bold" if top else "normal",
                    bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none"), arrowprops=dict(arrowstyle="-", color="#9a9a9a", lw=0.4, shrinkA=0, shrinkB=2.5), zorder=3)
    ax.set_xlim(0, 1.0); ax.set_ylim(0, 1.02); fticks(ax, "x"); fticks(ax, "y"); ax.set_xlabel("Precision"); ax.set_ylabel("Recall")
    g = SIG["external"][sample]; best = [k for k, v in TOOLKEY.items() if v == g["best"]][0]
    ax.text(0, 1.03, f"janus − {best}, F1: {g['diff']:+.3f} [{g['lo']:+.3f}, {g['hi']:+.3f}]", fontsize=FS - 0.8, color="#555555", va="bottom", transform=ax.transAxes)
    title(ax, ttl, 1.11)


def heat(ax, s, kind, cls, labs, cal, ttl):
    M = np.array([[pooled(s, n, kind, c)[0] for c in cls] for n in cal])
    ax.imshow(M, cmap=HEAT, vmin=0.4, vmax=1.0, aspect="auto")
    for i, n in enumerate(cal):
        for j, c in enumerate(cls):
            best = M[i, j] >= M[:, j].max() - 1e-9
            ax.text(j, i, f"{M[i, j]:.3f}", ha="center", va="center", fontsize=FS - 0.2, color=ink(HEAT((M[i, j] - 0.4) / 0.6), HEAT), fontweight="bold" if best else "normal")
    ax.set_xticks(range(len(cls))); ax.set_xticklabels([f"{l}\n({pooled(s, 'janus', kind, c)[1]:,})" for l, c in zip(labs, cls)], fontsize=FS - 1.4)
    ax.set_yticks(range(len(cal))); ax.set_yticklabels(cal, fontsize=FS)
    for t, n in zip(ax.get_yticklabels(), cal): t.set_color(JANUS if n == "janus" else "#333333"); t.set_fontweight("bold" if n == "janus" else "normal")
    ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    title(ax, ttl, 1.04)


def fbars(ax, s, kind, cls, labs, ttl):
    n = len(CAL5); w = 0.82 / n
    for k, nm in enumerate(CAL5):
        v = [pooled(s, nm, kind, c)[0] for c in cls]
        ax.bar(np.arange(len(cls)) + (k - (n - 1) / 2) * w, v, w * 0.9, color=color(nm), lw=0, zorder=2)
    ax.set_xticks(range(len(cls))); ax.set_xticklabels([f"{l}\n({pooled(s, 'janus', kind, c)[1]:,})" for l, c in zip(labs, cls)], fontsize=FS - 1)
    ax.set_ylim(0, 1); fticks(ax); ax.set_ylabel("F1 score")
    ax.set_xlim(-0.55, len(cls) - 0.45); ax.grid(axis="y", **GRID); ax.tick_params(axis="x", length=0); title(ax, ttl, 1.03)


def pfmt(p):
    return "<.001" if p < 0.001 else f"{p:.3f}"[1:] if p < 0.01 else f"{p:.2f}"[1:] if p < 0.995 else "1"


def full(n, k):
    return st.mean(D["all_with"][p][n][k] for p in P5)


def series(S, n, keys, pairs, j):
    t = TOOLKEY[n]
    return [full(n, ["P", "R", "F1"][j]) if k == "full" else st.mean(S[(p, k, t)][j] for p in pairs) for k in keys]


def f1_lines(ax, S, keys, labs, pairs, ttl, five=False, pvals=True, band=True):
    x = range(len(keys)); sg = SIG["settings"]
    if band: ax.fill_between(x, [sg[k]["snape_lo"] for k in keys], [sg[k]["snape_hi"] for k in keys], color=JANUS, alpha=0.25, lw=0, zorder=1)
    for n in ["janus"] + COV:
        ax.plot(x, series(S, n, keys, pairs, 2), color=color(n), lw=1.1 if n == "janus" else 1.0, marker=MARK[n], ms=3.6 if n == "janus" else 3.0, mfc=color(n), mec="white", mew=0.4, zorder=3 if n == "janus" else 2, label=n)
    for i, k in enumerate(keys if pvals else []):
        pv = sg[k]["p"]; ax.text(i, 1.07, pfmt(pv), ha="center", va="bottom", fontsize=FS - 1.6, color=JANUS if pv < 0.05 else "#777777", transform=ax.get_xaxis_transform())
    if pvals: ax.text(-0.04, 1.07, "p", ha="right", va="bottom", fontsize=FS - 1.4, style="italic", color="#555555", transform=ax.transAxes)
    ax.set_xticks(list(x)); ax.set_xticklabels(labs); ax.set_ylim(0, 1.0); fticks(ax)
    ax.set_ylabel("F1 score, mean"); title(ax, ttl, 1.19 if pvals else 1.03, 0.2); ax.grid(axis="y", **GRID)


def pr_lines(axp, axr, S, keys, labs, pairs):
    x = range(len(keys))
    for n in ["janus"] + COV:
        for a, j in ((axp, 0), (axr, 1)):
            a.plot(x, series(S, n, keys, pairs, j), color=color(n), lw=1.5 if n == "janus" else 0.9, marker=MARK[n], ms=2.8 if n == "janus" else 2.3, mfc=color(n), mec="white", mew=0.3, zorder=3 if n == "janus" else 2)
    for a, lab in ((axp, "Precision"), (axr, "Recall")):
        a.set_ylim(0, 1.05); a.set_yticks([0, 0.5, 1]); a.set_yticklabels(["0", ".5", "1"]); a.set_ylabel(lab); a.set_xticks(list(x)); a.grid(axis="y", **GRID)
    axp.set_xticklabels([]); axr.set_xticklabels(labs, fontsize=FS - 1)


def fig2():
    S = cov(); E = external()
    H = 150; T1 = H - 7; T2 = T1 - 34 - 17; T3 = T2 - 40 - 13
    fig, ax = page(160, H)
    order = [n for n in sorted(D["all_with"]["hcc1954"], key=lambda n: -mean("all_with", n, "F1")) if tech(n) == "HiFi"]
    a1 = ax(22, T1 - 34, 32, 34); f1_dots(a1, "all_with", order, "", True); title(a1, "Consensus of all 18 call sets", 1.03, 0.0); tag(a1, "a", 2, 4.8)
    ph = [Line2D([], [], ls="", marker="o", color=PAIRCOL[p], ms=3.4, label=PAIRNAME[p]) for p in P5] + [Line2D([], [], ls="", marker="|", color="black", ms=6, mew=0.9, label="mean of five pairs")]
    f = ax(88, T1 - 34, 34, 34); f1_lines(f, S, DEP, ["full", "30", "20", "15", "10", "5"], P5, "Depth, five pairs", pvals=False, band=False); f.set_xlabel("Depth of both samples (×)"); tag(f, "b", 78, 4.8)
    h = ax(124, T1 - 34, 34, 34); f1_lines(h, S, PUR, ["100", "80", "60", "40", "20"], P5, "Purity at 30×", pvals=False, band=False); h.set_ylabel(""); h.set_yticklabels([]); h.set_xlabel("Tumor purity (%)"); tag(h, "c", 117, 4.8)
    b = ax(23, T2 - 40, 40, 40); ext_scatter(b, E["colo829"], "colo829", "COLO829, 62 truth SVs"); tag(b, "d", 2, 9.6)
    c = ax(101, T2 - 40, 40, 40); ext_scatter(c, E["hg008"], "hg008", "HG008, 172 truth SVs"); tag(c, "e", 84, 9.6)
    d = ax(13, T3 - 20, 58, 20); fbars(d, "all_with", "type", TYPES, TYPES, "F1 by SV type, five CASTLE pairs pooled"); tag(d, "f", 2, 4.6)
    e = ax(88, T3 - 20, 70, 20); fbars(e, "all_with", "size", [k for k, _ in SIZES], [l for _, l in SIZES], "F1 by SV size (consensus SVs in parentheses)"); tag(e, "g", 87, 4.6)
    lh = [Line2D([], [], color=color(n), marker=MARK[n], ms=3.4, mec="white", lw=1.1 if n == "janus" else 1.0, label=n) for n in ["janus"] + COV]
    lh += [Patch(color=color("Sniffles2"), label="Sniffles2 (f, g)")]
    fig.legend(handles=lh, loc="lower center", bbox_to_anchor=(0.5, 5.0 / H), ncol=5, fontsize=FS - 0.7, handletextpad=0.4, columnspacing=1.2, frameon=False, borderaxespad=0)
    fig.legend(handles=[Patch(color="none", label="Pairs (a):")] + ph, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=7, fontsize=FS - 0.7, handletextpad=0.4, columnspacing=1.2, frameon=False, borderaxespad=0.3)
    fig.savefig(f"{OUT}/fig_accuracy.pdf", facecolor="white"); plt.close(fig)
    sg = SIG["settings"]; rows = [r + ["", ""] for r in dotrows("a", "all_with", order)]
    for pn, keys, labs in (("b", DEP, DL), ("c", PUR, PL)):
        for n in ["janus"] + COV:
            v = [series(S, n, keys, P5, j) for j in range(3)]
            rows += [[pn, n, "HiFi", l, v[0][i], v[1][i], v[2][i]] + [""] * 4 + ([sg[k]["snape_lo"], sg[k]["snape_hi"]] if n == "janus" else ["", ""]) for i, (k, l) in enumerate(zip(keys, labs))]
    for pn, q in (("d", "colo829"), ("e", "hg008")):
        rows += [[pn, n, "HiFi", PAIRNAME[q], *E[q][n]] + [""] * 6 for n in sorted(E[q], key=lambda n: -E[q][n][2])]
        g = SIG["external"][q]; rows.append([pn, f"janus minus {[k for k, v in TOOLKEY.items() if v == g['best']][0]}", "", "F1 difference", "", "", g["diff"]] + [""] * 4 + [g["lo"], g["hi"]])
    rows += [r + ["", ""] for r in poolrows("f", "all_with", "type", TYPES, TYPES, CAL5) + poolrows("g", "all_with", "size", [k for k, _ in SIZES], [l for _, l in SIZES], CAL5)]
    src("fig_accuracy", ACCH + ["ci_low", "ci_high"], rows)


def loo():
    L = {}
    for p in P5:
        L[p] = {}
        for k, n in EXT.items():
            r = list(csv.reader(open(f"{H}/eval_loo/{p}/{k}/score/results/overall.tsv"), delimiter="\t"))[1]
            L[p][n] = dict(P=float(r[4]), R=float(r[5]), F1=float(r[6]))
    return L


def supp_allsets():
    D["loo"] = loo()
    fig, ax = page(160, 123.5)
    order = sorted(D["all_with"]["hcc1954"], key=lambda n: -mean("all_with", n, "F1"))
    a1 = ax(29, 61, 38, 54); f1_dots(a1, "all_with", order, "Consensus of all 18 call sets", True); tag(a1, "a", 2, 4.5)
    a2 = ax(98, 61, 38, 54); f1_dots(a2, "all_without", order, "Consensus of 17, janus excluded", False); tag(a2, "b", 88, 4.5)
    a3 = ax(29, 9, 38, 33); f1_dots(a3, "loo", [n for n in order if tech(n) == "HiFi"], "Each call set against the consensus of the other 17", True, techs=True); tag(a3, "c", 2, 4.5)
    ph = [Line2D([], [], ls="", marker="o", color=PAIRCOL[p], ms=3.4, label=PAIRNAME[p]) for p in P5] + [Line2D([], [], ls="", marker="|", color="black", ms=6, mew=0.9, label="mean")]
    ph += [Patch(color=TECHCOL[t], label=l) for t, l in (("HiFi", "PacBio HiFi"), ("ONT", "ONT"), ("Illumina", "Illumina"))]
    fig.legend(handles=ph, loc="upper left", bbox_to_anchor=(98 / 160, 42 / 123.5), ncol=2, fontsize=FS - 0.5, handletextpad=0.4, columnspacing=1.6, handlelength=1.0, labelspacing=0.7)
    fig.savefig(f"{OUT}/fig_allsets.pdf", facecolor="white"); plt.close(fig)
    src("fig_allsets", ACCH, dotrows("a", "all_with", order) + dotrows("b", "all_without", order) + dotrows("c", "loo", [n for n in order if tech(n) == "HiFi"]))


def supp_without():
    fig, ax = page(160, 41.5)
    a = ax(17, 8, 55, 24); heat(a, "all_without", "type", TYPES, TYPES, CAL5, "F1 by SV type, consensus without janus"); tag(a, "a", 2, 5.2)
    b = ax(92, 8, 66, 24); heat(b, "all_without", "size", [k for k, _ in SIZES], [l for _, l in SIZES], CAL5, "F1 by SV size, consensus without janus"); tag(b, "b", 77, 5.2)
    fig.savefig(f"{OUT}/fig_accuracy_without.pdf", facecolor="white"); plt.close(fig)
    src("fig_accuracy_without", ACCH, poolrows("a", "all_without", "type", TYPES, TYPES, CAL5) + poolrows("b", "all_without", "size", [k for k, _ in SIZES], [l for _, l in SIZES], CAL5))


def supp_depth():
    S = cov(); xs = [5, 10, 15, 20, 30]; qs = [20, 40, 60, 80, 100]
    fig, ax = page(160, 118)
    ax = [ax(14 + 38.5 * (i % 4), 67 - 56 * (i // 4), 28, 38) for i in range(7)]
    for p in P5 + ["hg008", "colo829"]:
        ax[0].plot(xs, [S[(p, f"{x}x", "snape")][2] for x in xs], "-o", ms=2.2, lw=0.9, color=PAIRCOL[p], label=PAIRNAME[p])
    ax[0].set_xticks(xs); ax[0].set_xlabel("Depth of both samples (×)"); ax[0].set_ylabel("F1 score"); title(ax[0], "janus, depth", fs=0.0)
    for i, p in enumerate(P5):
        for n in ["janus"] + COV:
            ax[1 + i].plot(qs, [S[(p, f"{q}p", TOOLKEY[n])][2] for q in qs], "-", marker=MARK[n], ms=2.4, lw=1.2 if n == "janus" else 0.9, color=color(n), mec="white", mew=0.3, label=n)
        ax[1 + i].set_xticks(qs); ax[1 + i].set_xlabel("Tumor purity (%)"); ax[1 + i].set_ylabel("F1 score"); title(ax[1 + i], f"{PAIRNAME[p]}, purity", fs=0.0)
    for p in ["hg008", "colo829"]:
        ax[6].plot(qs, [S[(p, f"{q}p", "snape")][2] for q in qs], "-o", ms=2.2, lw=0.9, color=PAIRCOL[p])
    ax[6].set_xticks(qs); ax[6].set_xlabel("Tumor purity (%)"); ax[6].set_ylabel("F1 score"); title(ax[6], "janus, HG008 and COLO829", fs=0.0)
    for i, a in enumerate(ax):
        a.set_ylim(0, 1); fticks(a); a.grid(axis="y", **GRID); tag(a, "abcdefg"[i], 4 + 38.5 * (i % 4), 2.5)
        if i == 1: a.legend(loc="lower right", fontsize=FS - 1.2, handlelength=1.4)
    fig.legend(handles=[Line2D([], [], color=PAIRCOL[p], marker="o", ms=2.6, lw=0.9, label=PAIRNAME[p]) for p in P5 + ["hg008", "colo829"]], loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=7, fontsize=FS - 0.7, handlelength=1.4, columnspacing=1.2)
    fig.savefig(f"{OUT}/fig_depth.pdf", facecolor="white"); plt.close(fig)


def trace(p):
    def load(d, keep):
        steps = {s["step"]: (float(s["start"]), float(s["end"])) for s in csv.DictReader(open(f"{d}/full/{p}/steps.tsv"), delimiter="\t") if keep(s["step"])}
        t0 = min(v[0] for v in steps.values()); t1 = max(v[1] for v in steps.values())
        day = time.mktime(time.strptime(time.strftime("%Y-%m-%d", time.localtime(t0)), "%Y-%m-%d"))
        rss = collections.defaultdict(float)
        for l in open(f"{d}/full/{p}/monitor/pidstat.txt"):
            x = l.split()
            if len(x) < 19 or x[0] == "#" or x[-1] not in ("snape", "janus", "minimap2"): continue
            hh, mm, ss = [int(v) for v in x[0].split(":")]
            t = day + hh * 3600 + mm * 60 + ss
            if t < t0 - 60: t += 86400
            if t0 - 60 <= t <= t1 + 60: rss[t] += float(x[12]) / 2**20
        return steps, rss
    s1, r1 = load(R, lambda k: k.startswith("pack")); s2, r2 = load(RC, lambda k: not k.startswith("pack") and k != "total")
    steps = {**s1, **s2}; rss = {**r1, **r2}
    ts = np.array(sorted(rss))
    return steps, ts, np.array([rss[t] for t in ts])


def fig3(main=False):
    CJ = cost(); T, M = CJ["time"], CJ["mem"]
    ORDER = ["pack_n", "pack_t", "extract_n", "extract_t", "count_n", "count_t", "graph", "prune", "select", "place", "local", "gaps", "call"]
    LOG = lambda p, step: open(f"{RC}/full/{p}/logs/log.{step}").read()
    SD = []
    if not main:
        fig, ax = page(80, 52)
        g = ax(16, 11, 60, 36)
        STEPS = [("Tumor reads", lambda p: json.load(open(f"{R}/full/{p}/run.json"))["reads"]["t"]["reads"]), ("Sites", lambda p: int(re.search(r"dropped nodes, (\d+) sites", LOG(p, "prune")).group(1))),
                 ("Snippets", lambda p: int(re.search(r"reads; (\d+) pieces", LOG(p, "select")).group(1))), ("Records", lambda p: int(re.search(r"at sites, (\d+) records", LOG(p, "call")).group(1))),
                 ("PASS calls", lambda p: sum(1 for l in open(f"{H}/verdict/{p}.vcf") if "\tPASS\t" in l))]
        SD += [["", PAIRNAME[p], n, fn(p), "count"] for p in P5 for n, fn in STEPS]
        for p in P5: g.plot(range(len(STEPS)), [fn(p) for _, fn in STEPS], marker="o", ms=2.6, lw=0.9, color=PAIRCOL[p], mec="white", mew=0.3, label=PAIRNAME[p])
        g.set_yscale("log"); g.set_ylim(100, 3e7); g.set_xticks(range(len(STEPS))); g.set_xticklabels([n for n, _ in STEPS], rotation=30, ha="right", rotation_mode="anchor"); g.set_xlim(-0.3, len(STEPS) - 0.7)
        g.yaxis.set_minor_locator(matplotlib.ticker.NullLocator()); g.set_ylabel("Number per pair"); g.grid(axis="y", **GRID); g.legend(loc="lower left", fontsize=FS - 1.2, handlelength=1.2, labelspacing=0.3)
        fig.savefig(f"{OUT}/fig_cost.pdf", facecolor="white"); plt.close(fig)
        src("fig_cost", ["panel", "series", "item", "value", "unit"], SD); return
    ROWS = [("Severus", "severus_phased", 1), ("SAVANA", "savana_phased", 1), ("nanomonsv", "nanomonsv", 0), ("dysgu", "dysgu", 0), ("SVision-pro", "svision-pro", 0),
            ("cuteSV", "cutesv", 0), ("DeBreak", "debreak", 0), ("Sniffles2", "sniffles2", 0), ("kled", "kled", 0), ("SVDSS", "svdss", 0)]
    Y1, Y2 = 58, 12
    fig, ax = page(160, 102)
    # a: pipelines on HCC1954
    ALN = sum(json.load(open(f"{R}/align/hcc1954/{x}/run.json"))["total"]["wall_s"] for x in "nt")
    PH = {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(f"{R}/phase/hcc1954/times.tsv") if not l.startswith("step")}; phase = sum(PH.values())
    d = ax(19, Y1, 34, 36); sn = sum(T["hcc1954"].values()); pad = dict(xytext=(2.8, 0), textcoords="offset points", va="center")
    SD.append(["a", "janus", "total", sn / 3600, "h"]); d.barh(0, sn / 3600, color=JANUS, height=0.66); d.annotate(f"{sn / 60:.0f} min", (sn / 3600, 0), fontsize=FS, fontweight="bold", **pad)
    for k, (n, t, ph) in enumerate(ROWS, 1):
        w = json.load(open(f"{R}/competitors/{t}/hcc1954/run.json"))["total"]["wall_s"]; p_ = phase if ph else 0
        d.barh(k, ALN / 3600, color="#b9b9b9", height=0.66); d.barh(k, p_ / 3600, left=ALN / 3600, color="#e0a030", height=0.66); d.barh(k, w / 3600, left=(ALN + p_) / 3600, color="#3f6ea5", height=0.66)
        tot = (ALN + p_ + w) / 3600; SD += [["a", n, x_, v, "h"] for x_, v in (("mapping", ALN / 3600), ("phasing", p_ / 3600), ("caller", w / 3600), ("total", tot))] + [["a", n, "factor over janus", tot * 3600 / sn, "x"]]; d.annotate(f"{int(tot)}:{int(round((tot - int(tot)) * 60)):02d} h", (tot, k), fontsize=FS - 0.4, **pad); d.text(0.07, k, f"{tot * 3600 / sn:.0f}×", va="center", ha="left", fontsize=FS - 0.8, color=C["ink"], fontweight="bold")
    d.legend(handles=[Patch(color="#b9b9b9", label="mapping"), Patch(color="#e0a030", label="phasing"), Patch(color="#3f6ea5", label="caller")],
             loc="lower left", bbox_to_anchor=(-0.05, 1.0), ncol=3, fontsize=FS - 0.6, handlelength=1, handletextpad=0.4, columnspacing=0.8, frameon=False, borderaxespad=0.3)
    d.set_xlim(0, 7); d.set_xticks(range(8)); d.set_yticks(range(len(ROWS) + 1)); d.set_yticklabels(["janus"] + [r[0] for r in ROWS], fontsize=FS); d.get_yticklabels()[0].set_fontweight("bold"); d.get_yticklabels()[0].set_color(JANUS)
    d.set_ylim(len(ROWS) + 0.6, -0.6); d.set_xlabel("Wall clock on HCC1954, 16 threads (h)"); d.tick_params(axis="y", length=0); tag(d, "a", 2, 6.2)
    # b: memory trace
    TR = {p: trace(p) for p in P5}
    width = {k: st.median(TR[p][0][k][1] - TR[p][0][k][0] for p in P5) / 60 for k in ORDER}
    for g, ks in GROUP:
        ks = [k for k in ks if k in width]; tot = sum(width[k] for k in ks); tgt = st.median(T[p][g] for p in P5) / 60
        for k in ks: width[k] *= tgt / tot
    b = ax(74.5, Y1, 35, 36); x0 = 0
    for k in ORDER:
        w = width[k]; n = max(4, int(round(w * 30))); u = np.linspace(0.02, 0.98, n); Y = []
        for p in P5:
            steps, ts, y = TR[p]; s0, s1 = steps[k]; sel = (ts >= s0) & (ts <= s1)
            Y.append(np.interp(u, (ts[sel] - s0) / (s1 - s0), y[sel]) if sel.sum() > 1 else np.full(n, y[sel].max() if sel.any() else np.nan))
        Y = np.array(Y); xs = x0 + u * w
        SD += [["b", "resident memory, five pairs", k, x_, np.nanmedian(Y[:, i]), np.nanmin(Y[:, i]), np.nanmax(Y[:, i])] for i, x_ in enumerate(xs)]
        b.axvspan(x0, x0 + w, ymin=0.955, ymax=1, color=STEPCOL[k], lw=0)
        b.fill_between(xs, np.nanmin(Y, 0), np.nanmax(Y, 0), color="#9a9a9a", alpha=0.3, lw=0); b.plot(xs, np.nanmedian(Y, 0), color=C["ink"], lw=0.9)
        x0 += w; b.axvline(x0, color="#dddddd", lw=0.4, zorder=0)
    MAP = max(json.load(open(f"{R}/align/hcc1954/{x}/run.json"))["total"]["peak_rss_kib"] for x in "nt") / 2**20
    SD.append(["b", "minimap2 mapping HCC1954, peak", "", "", MAP, "", ""]); b.axhline(MAP, color="#444444", ls=(0, (4, 3)), lw=0.9); b.text(x0 * 0.985, MAP - 0.6, f"minimap2,\nHCC1954: {MAP:.1f} GB", ha="right", va="top", fontsize=FS - 1.2, zorder=6)
    b.set_xlim(0, x0); b.set_ylim(0, 22); b.set_yticks([0, 5, 10, 15, 20]); b.set_ylabel("Resident memory (GB)"); b.set_xlabel("Wall clock of the run (min)")
    tag(b, "b", 63.5, 6.2)
    fig.legend(handles=[Patch(color="none", label="Stages (strip in b, bars in d):")] + [Patch(color=STAGE[g], label=g) for g, _ in GROUP], loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=9, fontsize=FS - 1.2,
               handlelength=0.9, handletextpad=0.35, columnspacing=0.9, borderaxespad=0.3)
    # c: disk
    GBY = 2.0**30
    KIND = [("Packed reads", "#3f6ea5"), ("Core positions", "#86a8cc"), ("Counts", "#e0a030"), ("Graph", "#5fa37a"), ("Snippets, other", "#7e6bb0")]
    SZ, BAMS = [], []
    for p in P5:
        f = {x["file"]: x["bytes"] for x in json.load(open(f"{RC}/full/{p}/run.json"))["files"] if not x["file"].startswith(("log.", "steps", "calls.vcf"))}
        pk = sum(v for k, v in f.items() if k in ("n.svb", "t.svb", "n.svb.idx", "t.svb.idx")); ex = sum(v for k, v in f.items() if k.endswith(".exc"))
        cn = sum(v for k, v in f.items() if k.endswith(".cnt")); gr = f.get("graph.gph", 0)
        SZ.append([x / GBY for x in (pk, ex, cn, gr, sum(f.values()) - pk - ex - cn - gr)])
        BAMS.append([json.load(open(f"{R}/full/{p}/run.json"))["reads"][x]["bytes"] / GBY for x in "nt"])
    SZ = np.array(SZ).mean(0); BAMS = np.array(BAMS).mean(0)
    c = ax(125.5, Y1, 13, 36); bot = 0
    for (lab, col), v in zip(KIND, SZ): c.bar(0, v, bottom=bot, color=col, width=0.62, label=lab); bot += v
    c.text(0, bot + 3, f"{bot:.0f}", ha="center", fontsize=FS - 0.4)
    SD += [["c", "janus", lab, v, "GB"] for (lab, _), v in zip(KIND, SZ)] + [["c", "janus", "Total", bot, "GB"], ["c", "BAM", "normal", BAMS[0], "GB"], ["c", "BAM", "tumor", BAMS[1], "GB"], ["c", "BAM", "Total", BAMS.sum(), "GB"]]
    c.bar(1, BAMS[0], color="#b9b9b9", width=0.62, label="BAM, normal"); c.bar(1, BAMS[1], bottom=BAMS[0], color="#d6d6d6", width=0.62, label="BAM, tumor"); c.text(1, BAMS.sum() + 3, f"{BAMS.sum():.0f}", ha="center", fontsize=FS - 0.4)
    c.set_xticks([0, 1]); c.set_xticklabels(["janus", "BAM"]); c.set_xlim(-0.6, 1.6); c.set_ylim(0, 180); c.set_ylabel("Disk space, mean of pairs (GB)"); c.grid(axis="y", **GRID); c.set_axisbelow(True)
    c.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=FS - 1.3, handlelength=1, handletextpad=0.4, labelspacing=0.5, frameon=False, borderaxespad=0); c.tick_params(axis="x", length=0); tag(c, "c", 114.5, 6.2)
    # d: stage times per pair
    a = ax(19, Y2, 34, 30); rows = P5 + ["mean"]; tot = []
    for i, p in enumerate(rows):
        left = 0
        for g, _ in GROUP:
            v = (st.mean(T[q][g] for q in P5) if p == "mean" else T[p][g]) / 60
            a.barh(i, v, left=left, height=0.66, color=STAGE[g], lw=0); left += v; SD.append(["d", PAIRNAME.get(p, "mean of five pairs"), g, v, "min"])
        a.text(left + 0.15, i, f"{left:.1f}", va="center", fontsize=FS - 0.4, fontweight="bold" if p == "mean" else "normal"); tot.append(left); SD.append(["d", PAIRNAME.get(p, "mean of five pairs"), "Total", left, "min"])
    a.set_yticks(range(len(rows))); a.set_yticklabels([PAIRNAME.get(p, "mean") for p in rows], fontsize=FS)
    a.get_yticklabels()[-1].set_fontweight("bold"); a.set_ylim(len(rows) - 0.4, -0.6); a.set_xlim(0, 16); a.set_xticks([0, 5, 10, 15]); a.set_xlabel("Wall clock, 16 threads (min)")
    a.tick_params(axis="y", length=0); tag(a, "d", 2, 6.2)
    # e: sequence aligned
    e = ax(74.5, Y2, 35, 30); X = np.arange(5)
    rb = [sum(json.load(open(f"{R}/full/{p}/run.json"))["reads"][x]["bases"] for x in "nt") / 1e9 for p in P5]
    pb = [int(re.search(r"pieces, (\d+) bp", LOG(p, "select")).group(1)) / 1e9 for p in P5]
    e.bar(X - 0.19, rb, width=0.36, color="#b9b9b9", lw=0, zorder=2, label="minimap2, both read sets")
    e.bar(X + 0.19, pb, width=0.36, color=JANUS, lw=0, zorder=2, label="janus, snippets")
    SD += [["e", PAIRNAME[p], x_, v, u_] for i, p in enumerate(P5) for x_, v, u_ in (("minimap2, both read sets", rb[i], "Gbp"), ("janus, snippets", pb[i], "Gbp"), ("ratio", rb[i] / pb[i], "x"))]
    for i in range(5): e.text(i + 0.19, pb[i] * 1.5, f"{rb[i] / pb[i]:,.0f}×", ha="center", va="bottom", rotation=90, fontsize=FS - 1.2, color=JANUS)
    e.set_yscale("log"); e.set_ylim(0.1, 30000); e.set_yticks([0.1, 1, 10, 100, 1000]); e.set_yticklabels(["0.1", "1", "10", "100", "1,000"]); e.set_ylabel("Sequence aligned (Gbp)")
    e.set_xticks(X); e.set_xticklabels([PAIRNAME[p] for p in P5], rotation=30, ha="right", rotation_mode="anchor", fontsize=FS - 0.6); e.tick_params(axis="x", length=0); e.grid(axis="y", **GRID)
    e.legend(loc="upper left", fontsize=FS - 1.0, handlelength=1.0, handletextpad=0.4, borderaxespad=0.3, frameon=False); e.yaxis.set_minor_locator(matplotlib.ticker.NullLocator()); tag(e, "e", 63.5, 6.2)
    # f: CPU time
    f = ax(127, Y2, 29, 30); ACPU = sum(json.load(open(f"{R}/align/hcc1954/{x}/run.json"))["total"]["cpu_s"] for x in "nt") / 3600
    PCPU = sum(float(l.split("\t")[2]) for l in open(f"{R}/phase/hcc1954/times.tsv") if not l.startswith("step")) / 3600
    sc = (sum(x["cpu_s"] for x in json.load(open(f"{R}/full/hcc1954/run.json"))["steps"] if x["step"].startswith("pack")) + json.load(open(f"{RC}/full/hcc1954/run.json"))["total"]["cpu_s"]) / 3600
    SD.append(["f", "janus", "total", sc, "h"]); f.barh(0, sc, color=JANUS, height=0.66); f.text(sc + 1.5, 0, f"{sc:.1f} h", va="center", fontsize=FS - 0.4, fontweight="bold")
    for k, (n, t, ph) in enumerate(ROWS, 1):
        cc = json.load(open(f"{R}/competitors/{t}/hcc1954/run.json"))["total"]["cpu_s"] / 3600; pc = PCPU if ph else 0
        f.barh(k, ACPU, color="#b9b9b9", height=0.66); f.barh(k, pc, left=ACPU, color="#e0a030", height=0.66); f.barh(k, cc, left=ACPU + pc, color="#3f6ea5", height=0.66)
        SD += [["f", n, x_, v, "h"] for x_, v in (("mapping", ACPU), ("phasing", pc), ("caller", cc), ("total", ACPU + pc + cc))]
        f.text(ACPU + pc + cc + 1.5, k, f"{ACPU + pc + cc:.0f} h", va="center", fontsize=FS - 0.4)
    f.set_xlim(0, 105); f.set_xticks([0, 25, 50, 75, 100]); f.set_yticks(range(len(ROWS) + 1)); f.set_yticklabels(["janus"] + [r[0] for r in ROWS], fontsize=FS); f.get_yticklabels()[0].set_fontweight("bold"); f.get_yticklabels()[0].set_color(JANUS)
    f.set_ylim(len(ROWS) + 0.6, -0.6); f.set_xlabel("CPU time on HCC1954 (h)"); f.tick_params(axis="y", length=0); f.grid(axis="x", **GRID); f.set_axisbelow(True); tag(f, "f", 114.5, 6.2)
    fig.savefig(f"{OUT}/fig_runtime.pdf", facecolor="white"); plt.close(fig)
    src("fig_runtime", ["panel", "series", "item", "x_min", "value", "low", "high", "unit"], [r[:3] + ["", r[3], "", "", r[4]] if r[0] != "b" else r + ["GB"] for r in SD])


def tree(ax, clades, wrong, ttl, mark=()):
    clades = [frozenset(c) for c in clades] + [frozenset([i]) for i in range(20)] + [frozenset(range(20))]
    clades = sorted(set(clades), key=len)
    par = {c: next(d for d in clades if c < d) for c in clades if len(c) < 20}
    kids = {}
    for c, p in par.items(): kids.setdefault(p, []).append(c)
    order, depth, y = [], {}, {}

    def walk(c, d):
        depth[c] = d
        if len(c) == 1:
            y[c] = len(order); order.append(c); return
        for k in sorted(kids[c], key=min): walk(k, d + 1)
        y[c] = sum(y[k] for k in kids[c]) / len(kids[c])

    root = clades[-1]; walk(root, 0); deep = max(depth.values())
    for c, p in par.items():
        x0, x1 = depth[p], deep if len(c) == 1 else depth[c]
        col = "#e0a030" if c in wrong else JANUS if c in mark else C["ink"]
        ax.plot([x0, x0, x1], [y[p], y[c], y[c]], color=col, lw=1.3 if (c in wrong or c in mark) else 0.6, solid_joinstyle="miter", zorder=3 if (c in wrong or c in mark) else 2)
    ax.plot([-0.6, 0], [y[root]] * 2, color=C["ink"], lw=0.6)
    for c in order: ax.text(deep + 0.25, y[c], "T%02d" % (min(c) + 1), fontsize=FS - 1, va="center")
    ax.set_ylim(19.7, -0.7); ax.set_xlim(-0.7, deep + 1.8); ax.axis("off"); ax.text(0.08, 0.99, ttl, transform=ax.transAxes, fontsize=FS + 0.8, va="bottom")


def fig4():
    S = simdata(); TR, MP = simtrees(); M2 = mm2()
    SETS = ["20-20", "20-10", "20-5", "10-10", "5-5"]; LAB = [s.replace("-", "/") for s in SETS]; BASE = "20-5"
    ONE, ALL = "#8c8c8c", JANUS
    fig, ax = page(160, 125)
    true = [frozenset(c) for c in TR["true"]]; rec = MP[BASE]
    kids = collections.defaultdict(list)
    for n, q in TR["parent"].items(): kids[q].append(n)

    def leaves(n): return {int(n[1:]) - 1} if n.startswith("T") else set().union(*(leaves(k) for k in kids[n]))

    miss55 = {frozenset(leaves(n)) for n in MP["5-5"].get("missing", [])}
    at = ax(70, 4, 40, 42); tree(at, true, set(), "Simulated tree", mark=miss55); tag(at, "d", 64, 2.5)
    ar = ax(116, 4, 40, 42); tree(ar, rec["clades"], {frozenset(c) for c in rec["wrong"]}, f"Reconstructed at {BASE.replace('-', '/')}")
    if miss55: at.plot([], [], color=JANUS, lw=1.3, label="the one split missed at 5/5"); at.legend(loc="upper left", bbox_to_anchor=(0.0, 0.0), fontsize=FS - 0.8, handlelength=1.5, frameon=False)
    ab = ax(14, 10, 44, 36); PO = json.load(open(f"{R}/sim/pooled/scores.json"))[BASE.replace("-", "/").replace("/", "-")]; CL = list(PO["joint"]["classes"]); w = 0.26
    POOLED = "#e0a030"
    for j, (key, col, lab) in enumerate((("pooled", POOLED, "20 tumors pooled into one sample"), ("pairwise", ONE, "20 pairwise runs"), ("joint", ALL, "joint run"))):
        ab.bar(np.arange(len(CL)) + (j - 1) * w, [PO[key]["classes"][c][1] for c in CL], w * 0.92, color=col, lw=0, zorder=2, label=lab)
    ab.set_ylim(0, 1); fticks(ab); ab.set_xticks(range(len(CL))); ab.set_xticklabels([f"{c}\n({PO['joint']['classes'][c][0]:,})" for c in CL], fontsize=FS - 0.6); ab.set_xlabel("Tumors carrying the SV (SVs)"); ab.set_ylabel("Recall of the cohort's SVs"); ab.grid(axis="y", **GRID); ab.tick_params(axis="x", length=0)
    ab.legend(loc="upper left", fontsize=FS - 1.0, handlelength=1.0, handletextpad=0.4, borderaxespad=0.3, frameon=False, labelspacing=0.3); title(ab, f"At {BASE.replace('-', '/')}", 1.03, 0.0); tag(ab, "c", 2, 2.5)

    def f1(d): return [2 * p * r / (p + r) for p, r in zip(d["P"], d["R"])]

    ax = [ax(12, 64, 60, 52)]
    for a, get, lab in zip(ax, (f1,), ("F1 score",)):
        for i, s in enumerate(SETS):
            u, v = get(S[s]["pair"]), get(S[s]["joint"])
            for x, yv in zip(u, v): a.plot([i - 0.18, i + 0.18], [x, yv], color="#cfcfcf", lw=0.35, zorder=1)
            a.plot([i - 0.18] * 20, u, "o", color=ONE, ms=1.9, mew=0, zorder=2, label="20 pairwise runs" if i == 0 else None)
            a.plot([i + 0.18] * 20, v, "o", color=ALL, ms=1.9, mew=0, zorder=2, label="joint run" if i == 0 else None)
            a.plot([i - 0.36, i - 0.0], [st.median(u)] * 2, color="#3a3a3a", lw=2.0, solid_capstyle="butt", zorder=4); a.plot([i + 0.0, i + 0.36], [st.median(v)] * 2, color="#7d1631", lw=2.0, solid_capstyle="butt", zorder=4)
            a.text(i, 0.315, f"{st.median(v):.2f}", fontsize=FS - 1, color=ALL, va="bottom", ha="center")
        a.set_ylim(0.3, 1.0); fticks(a, "y", np.arange(0.3, 1.01, 0.1)); a.set_xlim(-0.5, 4.5); a.set_xticks(range(5)); a.set_xticklabels(LAB); a.set_xlabel("Depth, normal/tumor (×)"); a.set_ylabel(lab + ", per tumor"); a.grid(axis="y", **GRID)
    ax[0].legend(handles=[Line2D([], [], ls="", marker="o", color=ONE, ms=3.8, mew=0, label="20 pairwise runs"), Line2D([], [], ls="", marker="o", color=ALL, ms=3.8, mew=0, label="joint run"),
                          Line2D([], [], color="#3a3a3a", lw=2.0, label="median, pairwise"), Line2D([], [], color="#7d1631", lw=2.0, label="median, joint")],
                 loc="lower left", bbox_to_anchor=(0.0, 0.1), handlelength=1.2, borderaxespad=0.2, fontsize=FS - 0.6, ncol=2, columnspacing=1.0)
    tag(ax[0], "a", 2, 2.5)
    aw = fig.add_axes([92 / 160, 64 / 125, 66 / 160, 52 / 125]); y = np.arange(5)[::-1]; KM = 3.2
    for i, (s_, l) in enumerate(zip(SETS, LAB)):
        t = S[s_]["time"]; x0, y0, x1, y1 = sum(t["pair_wall"]) / 60, max(t["pair_peak_kb"]) / 2**20, t["joint_wall"] / 60, t["joint_peak_kb"] / 2**20
        aw.barh(y[i] + 0.19, x0, 0.34, color="#c4c4c4", lw=0); aw.barh(y[i] - 0.19, x1, 0.34, color="#f2b8c1", lw=0)
        aw.barh(y[i] + 0.19, -y0 * KM, 0.34, color=ONE, lw=0); aw.barh(y[i] - 0.19, -y1 * KM, 0.34, color=ALL, lw=0)
        aw.text(x0 + 1.2, y[i] + 0.19, f"{x0:.0f}", va="center", fontsize=FS - 1, color="#555555"); aw.text(x1 + 1.2, y[i] - 0.19, f"{x1:.0f}", va="center", fontsize=FS - 1, color=ALL)
        aw.text(-y0 * KM - 1.2, y[i] + 0.19, f"{y0:.1f}", va="center", ha="right", fontsize=FS - 1, color="#555555"); aw.text(-y1 * KM - 1.2, y[i] - 0.19, f"{y1:.1f}", va="center", ha="right", fontsize=FS - 1, color=ALL)
    aw.set_yticks(y); aw.set_yticklabels(LAB); aw.axvline(0, color=C["ink"], lw=0.6); aw.tick_params(axis="y", length=0)
    aw.set_xticks([-20 * KM, -10 * KM, 0, 20, 40, 60]); aw.set_xticklabels(["20", "10", "0", "20", "40", "60"])
    aw.set_ylim(-0.75, 4.75); aw.set_xlim(-24 * KM, 76); aw.set_ylabel("Depth, normal/tumor (×)"); aw.grid(axis="x", **GRID)
    aw.text(-12 * KM, 5.0, "Peak memory (GB)", ha="center", va="bottom", fontsize=FS); aw.text(38, 5.0, "Wall clock, 16 threads (min)", ha="center", va="bottom", fontsize=FS)
    if "5-5" in M2: aw.text(76, -0.72, f"minimap2, 21 read sets at 5/5: {M2['5-5']['wall_min']:.0f} min, {M2['5-5']['peak_gb']:.1f} GB", ha="right", va="bottom", fontsize=FS - 1.2, color="#555555")
    aw.legend(handles=[Patch(color=ONE, label="20 pairwise runs (sum of times, largest memory)"), Patch(color=ALL, label="joint run")], loc="upper center", bbox_to_anchor=(0.45, -0.1), ncol=2, fontsize=FS - 0.8, handlelength=1.0, handletextpad=0.4, columnspacing=1.0, frameon=False, borderaxespad=0)
    for sp in ("left",): aw.spines[sp].set_visible(False)
    tag(aw, "b", 80, 6.0)
    fig.savefig(f"{OUT}/fig_sim.pdf", facecolor="white"); plt.close(fig)
    SD = []
    for s, l in zip(SETS, LAB):
        for m, ser in (("pair", "pairwise run"), ("joint", "joint run")):
            q = S[s][m]; v = f1(q); SD += [["a", l, ser, t, q["P"][i], q["R"][i], v[i], "", ""] for i, t in enumerate(SIMT)] + [["a", l, ser, "median", st.median(q["P"]), st.median(q["R"]), st.median(v), "", ""]]
        t = S[s]["time"]; d = TR["sets"][s]
        SD += [["b", l, "janus, 20 pairwise runs (sum)", "wall clock (min)", "", "", "", sum(t["pair_wall"]) / 60, ""], ["b", l, "janus, joint run", "wall clock (min)", "", "", "", t["joint_wall"] / 60, ""],
               ["b", l, "janus, pairwise runs (largest)", "peak memory (GB)", "", "", "", max(t["pair_peak_kb"]) / 2**20, ""], ["b", l, "janus, joint run", "peak memory (GB)", "", "", "", t["joint_peak_kb"] / 2**20, ""],
               ["d", l, "joint run", "splits recovered", "", "", "", MP[s]["shared"], f"of {MP[s]['true']}"]]
    for key, ser in (("pooled", "20 tumors pooled into one sample"), ("pairwise", "20 pairwise runs"), ("joint", "joint run")):
        o = PO[key]; SD += [["c", BASE.replace("-", "/"), ser, f"SVs carried by {c} tumors", "", "", "", v[1], f"{v[0]} SVs"] for c, v in o["classes"].items()]
        SD.append(["c", BASE.replace("-", "/"), ser, "all distinct SVs", o["P"], o["R"], o["F1"], o["found"], f"{o['calls']} calls, {o['tp']} true, {PO['truth']} SVs"])
    for k, v in M2.items():
        SD += [["b", k.replace("-", "/"), "minimap2, 21 read sets", "wall clock (min)", "", "", "", v["wall_min"], ""], ["b", k.replace("-", "/"), "minimap2, 21 read sets", "peak memory (GB)", "", "", "", v["peak_gb"], ""]]
    SD += [["d", "", "simulated tree", ";".join("T%02d" % (i + 1) for i in sorted(c)), "", "", "", "", "split missed at 5/5" if c in miss55 else ""] for c in sorted(true, key=lambda c: (len(c), sorted(c)))]
    SD += [["d", BASE.replace("-", "/"), "reconstructed tree", ";".join("T%02d" % (i + 1) for i in c), "", "", "", "", "not in the simulated tree" if c in rec["wrong"] else ""] for c in sorted(rec["clades"], key=lambda c: (len(c), c))]
    src("fig_sim", ["panel", "setting", "series", "item", "precision", "recall", "f1", "value", "note"], SD)



def gridvals():
    S, O = {}, {}
    for l in open(f"{H}/runs/coverage_eval/scores.tsv"):
        s, t, n, pp, r, f = l.split()
        p, k = s.rsplit("_", 1)
        if t == "snape":
            S[(p, k)] = float(f)
        else:
            O[(p, k)] = max(O.get((p, k), 0), float(f))
    E = external()
    for p in P5:
        S[(p, "full")] = D["all_with"][p]["janus"]["F1"]
        O[(p, "full")] = max(D["all_with"][p][n]["F1"] for n in ("Severus", "SAVANA", "nanomonsv"))
    for p in ("hg008", "colo829"):
        S[(p, "full")] = E[p]["janus"][2]
        O[(p, "full")] = max(E[p][n][2] for n in ("Severus", "SAVANA", "nanomonsv"))
    return S, O


def grid():
    S, O = gridvals()
    J = SIG
    SIGP = {(p, k): (v['lo'], v['hi']) for k, x in J['settings'].items() for p, v in x['per_pair'].items()}
    SIGP.update({(p, 'full'): (v['lo'], v['hi']) for p, v in J['external'].items()})
    cols = DEP + PUR
    head = ["full"] + ["30×", "20×", "15×", "10×", "5×"] + ["100%", "80%", "60%", "40%", "20%"]
    fig, ax = page(W, 63.5)
    a = ax(15, 1.5, 124, 62)
    lim = 0.25; GC = LinearSegmentedColormap.from_list("g", [(0, DIV(0.4)), (1 / 6, DIV(0.5)), (1, DIV(1.0))]); NG = plt.Normalize(-0.05, lim, clip=True)
    for i, p in enumerate(ROWS):
        for j, k in enumerate(cols):
            x = j + (j > 0) * 0.25 + (j > 5) * 0.25
            s = S.get((p, k))
            if s is None:
                continue
            o = O.get((p, k))
            if o is None:
                a.add_patch(plt.Rectangle((x, i), 0.94, 0.9, facecolor="white", edgecolor=LIGHT, lw=0.4, hatch="////"))
                a.text(x + 0.47, i + 0.45, f"{s:.2f}", ha="center", va="center", fontsize=FS - 0.6, color=C["ink"], bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"))
                continue
            d = s - o
            fc = GC(NG(d)); a.add_patch(plt.Rectangle((x, i), 0.94, 0.9, facecolor=fc, edgecolor="none"))
            q = SIGP.get((p, k))
            sym = "=" if q is None else "▲" if q[0] > 0 else "▼" if q[1] < 0 else "="
            a.text(x + 0.47, i + 0.29, f"{s:.2f}", ha="center", va="center", fontsize=FS - 0.6, color=ink(fc, GC), fontweight="bold")
            a.text(x + 0.47, i + 0.66, f"{sym} {d:+.2f}".replace("-", "\u2212"), ha="center", va="center", fontsize=FS - 1, color=ink(fc, GC) if ink(fc, GC) == "white" else C["ink2"])
    for j, h in enumerate(head):
        x = j + (j > 0) * 0.25 + (j > 5) * 0.25
        a.text(x + 0.47, -0.12, h, ha="center", va="bottom", fontsize=FS - 1)
    a.text(1.25 + 2.5, -0.62, "Depth (tumor and normal)", ha="center", va="bottom", fontsize=FS)
    a.text(6.5 + 2.5, -0.62, "Tumor purity (30×)", ha="center", va="bottom", fontsize=FS)
    for i, p in enumerate(ROWS):
        a.text(-0.12, i + 0.45, PAIRNAME[p], ha="right", va="center", fontsize=FS - 0.6)
    a.set_xlim(-0.05, 11.5)
    a.set_ylim(len(ROWS), -1.2)
    a.axis("off")
    c = ax(143, 6, 2.4, 46)
    sm = plt.cm.ScalarMappable(cmap=GC, norm=NG)
    cb = fig.colorbar(sm, cax=c, ticks=[-0.05, 0, 0.1, 0.2]); cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(ffmt))
    cb.ax.tick_params(length=1.5, pad=1)
    cb.outline.set_linewidth(0.4)
    cb.set_label("janus F1 minus best of\nSeverus, SAVANA, nanomonsv", fontsize=FS - 1.4, labelpad=2)
    fig.savefig(f"{OUT}/fig_grid.pdf", facecolor="white")
    plt.close(fig)
    C_, E = cov(), external(); SD = []
    pr_ = lambda p, k, n: (tuple(D["all_with"][p][n][x] for x in ("P", "R", "F1")) if p in P5 else E[p][n]) if k == "full" else C_.get((p, k, TOOLKEY[n]))
    for p in ROWS:
        for k, h in zip(cols, head):
            if (p, k) not in S: continue
            o = [n for n in COV if pr_(p, k, n)]; b = max(o, key=lambda n: pr_(p, k, n)[2]) if o else None; q = SIGP.get((p, k))
            SD.append([PAIRNAME[p], h, *pr_(p, k, "janus"), b or "", *(pr_(p, k, b) if b else ("", "", "")), S[(p, k)] - O[(p, k)] if b else "", *(q or ("", "")),
                       ("above" if q[0] > 0 else "below" if q[1] < 0 else "contains 0") if q and b else ""])
    src("fig_grid", ["pair", "setting", "janus_precision", "janus_recall", "janus_f1", "best_other", "other_precision", "other_recall", "other_f1", "f1_difference", "ci_low", "ci_high", "interval_vs_0"], SD)


def featfig():
    J = json.load(open(f"{H}/supp/features.json"))
    G = [(n, ix) for n, ix in J["groups"] if n != "Minimum and count-stage context"]
    P7 = list(J["pairs"])
    M = {n: [float(np.mean(J["pairs"][p]["drop"][n])) for p in P7] for n, _ in G}
    order = sorted(M, key=lambda n: np.mean(M[n]))
    SHOW = {"Tumor reads, single anchor pair": "Tumor reads, anchors of the site", "Normal reads, single anchor pair": "Normal reads, anchors of the site",
            "Tumor reads, anchor grid and paths": "Tumor reads, anchors next to the event", "Normal reads, anchor grid and paths": "Normal reads, anchors next to the event",
            "Graph context": "Graph", "Alignment quality": "Alignment (mapping quality, divergence)", "Reference context and depth": "Reference context and depth", "Type and size": "Type and size"}
    fig, ax = page(W, 63)
    a = ax(51, 9.5, 105, 47)
    for i, n in enumerate(order):
        a.barh(i, np.mean(M[n]), height=0.62, color=LIGHT, lw=0, zorder=2)
        for k, p in enumerate(P7):
            a.plot(M[n][k], i + (k - 3) * 0.09, "o", color=PAIRCOL[p], ms=3.0, mec="white", mew=0.3, zorder=3)
    a.axvline(0, color=C["ink2"], lw=0.5); a.spines["left"].set_visible(False); a.set_xlim(-0.035, max(max(v) for v in M.values()) + 0.03); fticks(a, "x", np.arange(0, 0.71, 0.1))
    a.set_yticks(range(len(order)))
    a.set_yticklabels([f"{SHOW[n]} ({len(dict(G)[n])})" for n in order])
    a.tick_params(axis="y", length=0)
    a.set_xlabel("Drop in F1 of the pair when the group is permuted")
    a.set_ylim(-0.6, len(order) - 0.4)
    a.grid(axis="x", **GRID)
    hand = [Patch(color=LIGHT, label="mean of seven pairs")] + [Line2D([], [], ls="", marker="o", color=PAIRCOL[p], mec="white", mew=0.3, ms=3.4, label=PAIRNAME[p]) for p in P7]
    fig.legend(handles=hand, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=8, handlelength=1.0, handletextpad=0.3, columnspacing=1.0)
    fig.savefig(f"{OUT}/fig_features.pdf", facecolor="white")
    plt.close(fig)
    src("fig_features", ["feature_group", "features", "pair", "f1_drop"], [[SHOW[n], len(dict(G)[n]), PAIRNAME[p], M[n][k]] for n in order for k, p in enumerate(P7)] + [[SHOW[n], len(dict(G)[n]), "mean of seven pairs", np.mean(M[n])] for n in order])


def clsfig():
    F = json.load(open(f"{H}/verdict/folds.json"))
    P7 = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395", "hg008", "colo829"]
    CUTS = np.arange(0.2, 0.81, 0.05)
    def f1(p, y, c):
        return 2 * ((p >= c) & (y > 0)).sum() / max((p >= c).sum() + (y > 0).sum(), 1)
    def load(p):
        d = np.array(F[p]["dupl"]) == 0
        return np.array(F[p]["prob"])[d], np.array(F[p]["label"])[d]
    fig, ax = page(W, 64)
    a = ax(14, 15, 64, 42)
    pr = np.concatenate([load(p)[0] for p in P5]); y = np.concatenate([load(p)[1] for p in P5])
    B = np.linspace(0, 1, 21); SD = []
    for v, c, n in ((1, MATCH, "Match the training labels"), (0, NOMATCH, "Do not match")):
        h = np.histogram(pr[y == v], B)[0]; SD += [["a", n, B[i], B[i + 1], int(h[i]), "", "", "", "histogram"] for i in range(len(h))]
        a.stairs(h, B, baseline=None, color=c, lw=1.3, label=f"{n} ({(y == v).sum():,})")
        print(n, (y == v).sum(), "p<0.1", (pr[y == v] < 0.1).mean().round(3), "p>=0.9", (pr[y == v] >= 0.9).mean().round(3), "0.2<=p<0.8", ((pr[y == v] >= 0.2) & (pr[y == v] < 0.8)).mean().round(3))
    a.set_yscale("log"); a.set_xlim(0, 1); a.set_ylim(8, 2e4); fticks(a, "x", np.arange(0, 1.01, 0.2))
    a.set_xlabel("Held-out probability"); a.set_ylabel("Records, five pairs pooled"); a.legend(loc="upper center", handlelength=1.4)
    tag(a, "a", 2, 2.5)
    a = ax(94, 15, 64, 42)
    for p in P7:
        pp, yy = load(p); cut = F[p]["cut"]
        v = [f1(pp, yy, t) for t in CUTS]
        SD += [["b", PAIRNAME[p], t, "", int((pp >= t).sum()), ((pp >= t) & (yy > 0)).sum() / max((pp >= t).sum(), 1), ((pp >= t) & (yy > 0)).sum() / max((yy > 0).sum(), 1), f1(pp, yy, t), "threshold used" if m_ else "curve"] for t, m_ in [(t, 0) for t in CUTS] + [(cut, 1)]]
        a.plot(CUTS, v, color=PAIRCOL[p], lw=1.0, label=PAIRNAME[p])
        a.plot([cut], [f1(pp, yy, cut)], "o", color=PAIRCOL[p], ms=3.6, mec="white", mew=0.5, zorder=3)
        b = int(np.argmax(v))
        print(f"{p} cut {cut} PASS {int((pp >= cut).sum())} F1 at cut {f1(pp, yy, cut):.3f} best {v[b]:.3f} at {CUTS[b]:.2f} range {min(v):.3f}-{max(v):.3f}")
    a.set_xlim(0.18, 0.82); a.set_ylim(0.8, 0.96); fticks(a, "x", [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]); fticks(a, "y", [0.8, 0.85, 0.9, 0.95])
    a.set_xlabel("Probability threshold"); a.set_ylabel("F1 of the held-out pair"); a.grid(axis="y", **GRID)
    fig.legend(handles=a.get_legend_handles_labels()[0], loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=7, handlelength=1.4, columnspacing=1.2, borderaxespad=0.3)
    tag(a, "b", 82, 2.5)
    fig.savefig(f"{OUT}/fig_classifier.pdf", facecolor="white")
    plt.close(fig)
    src("fig_classifier", ["panel", "series", "x", "x_high", "records", "precision", "recall", "f1", "kind"], SD)


def missedfig():
    ns = {}
    exec(paper.SCRIPTS["score"].split("T, C = recs")[0], ns)
    hit, ALT = ns["hit"], ns["ALT"]
    KIND = {}; miss = {}; bt = collections.defaultdict(lambda: [0, 0, 0])
    for p in P5:
        T = []
        for l in open(f"{H}/eval_sets/{p}/all_with/None_minda_ensemble.vcf"):
            if l[0] == "#":
                continue
            f = l.rstrip("\n").split("\t")
            i = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            m = ALT.search(f[4])
            r = (i.get("SVTYPE"), f[0], int(float(f[1])), m.group(1), int(m.group(2))) if m else (i.get("SVTYPE"), f[0], int(float(f[1])), f[0], int(float(i.get("END", f[1]))))
            T.append(r)
            KIND[r] = "INV" if r[0] == "INV" or (r[1] == i.get("CHR2", r[3]) and i.get("STRANDS") in ("++", "--")) else r[0]
        PS, LO = [], []
        for l in open(f"{H}/verdict/{p}.vcf"):
            if l[0] == "#":
                continue
            f = l.split("\t")
            i = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            m = ALT.search(f[4])
            r = (i["SVTYPE"], f[0], int(f[1]), m.group(1), int(m.group(2))) if m else (i["SVTYPE"], f[0], int(f[1]), f[0], int(i["END"]))
            (PS if f[6] == "PASS" else LO).append(r)
        fn = [t for t in T if not any(hit(t, c) for c in PS)]
        lo = 0
        for t in T:
            bt[KIND[t]][2] += 1
        for t in fn:
            h = any(hit(t, c) for c in LO)
            lo += h; bt[KIND[t]][0 if h else 1] += 1
        miss[p] = (len(T), len(fn), lo, len(fn) - lo)
        print(PAIRNAME[p], miss[p])
    VAL["missed"] = miss
    print("total", [sum(miss[p][i] for p in P5) for i in range(4)])
    print("by type (not PASS, no record, consensus SVs)", dict(bt))
    LA, LB = "Record that is not PASS", "No record"
    fig, ax = page(W, 65.5)


    def panel(ax, keys, labs, A, B):
        ax.bar(range(len(keys)), A, width=0.62, color=FILTER["LOW"], lw=0, zorder=2, label=LA)
        ax.bar(range(len(keys)), B, bottom=A, width=0.62, facecolor="white", edgecolor=NEUTRAL, lw=0.5, zorder=2, label=LB)
        top = max(x + y for x, y in zip(A, B))
        for i, (x, y) in enumerate(zip(A, B)):
            if x > top * 0.06:
                ax.text(i, x / 2, str(x), ha="center", va="center", fontsize=FS - 1.4, color="white")
            if y > top * 0.06:
                ax.text(i, x + y / 2, str(y), ha="center", va="center", fontsize=FS - 1.4, color=C["ink"])
            ax.text(i, x + y + top * 0.015, str(x + y), ha="center", va="bottom", fontsize=FS - 1.4, color=C["ink2"])
        ax.set_xticks(range(len(keys))); ax.set_xticklabels(labs); ax.set_ylim(0, top * 1.12)
        ax.set_ylabel("Consensus SVs without a PASS call"); ax.tick_params(axis="x", length=0); ax.grid(axis="y", **GRID)


    a = ax(14, 11.5, 62, 47)
    panel(a, P5, [f"{PAIRNAME[p]}\n({miss[p][0]:,})" for p in P5], [miss[p][2] for p in P5], [miss[p][3] for p in P5])
    a.set_xlabel("Pair (consensus SVs)"); a.legend(loc="upper left", handlelength=1.2, frameon=True, facecolor="white", edgecolor="none", framealpha=1, borderpad=0.3)
    tag(a, "a", 2, 2.5)
    TY = sorted(bt, key=lambda t: (t not in ("BND", "DEL", "DUP", "INS", "INV"), t))
    a = ax(96, 11.5, 62, 47)
    panel(a, TY, [f"{t}\n({bt[t][2]:,})" for t in TY], [bt[t][0] for t in TY], [bt[t][1] for t in TY])
    a.set_xlabel("SV type (consensus SVs, five pairs pooled)")
    tag(a, "b", 84, 2.5)
    fig.savefig(f"{OUT}/fig_missed.pdf", facecolor="white")
    plt.close(fig)
    src("fig_missed", ["panel", "item", "consensus_svs", "record_not_pass", "no_record", "without_pass_call"], [["a", PAIRNAME[p], miss[p][0], miss[p][2], miss[p][3], miss[p][1]] for p in P5] + [["b", t, bt[t][2], bt[t][0], bt[t][1], bt[t][0] + bt[t][1]] for t in TY])


def evidencefig():
    F = json.load(open(f"{H}/verdict/folds.json"))
    N, N26, T, Y = [], [], [], []
    for p in P5:
        k = 0
        for l in open(f"{H}/verdict/{p}.vcf"):
            if l[0] == "#": continue
            f = l.rstrip().split("\t"); i = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if not F[p]["dupl"][k]:
                x = [float(v) for v in i["F"].split(",")]
                N.append(int(i["SAME"].split(",")[0])); N26.append(round(2 ** x[26] - 1)); T.append(round(2 ** x[1] - 1)); Y.append(F[p]["label"][k])
            k += 1
    N, N26, T, Y = [np.array(v) for v in (N, N26, T, Y)]
    fig, ax = page(W, 60)
    a = ax(14, 9.5, 62, 44)
    for j, (v, c, n) in enumerate(((1, MATCH, "Match the training labels"), (0, NOMATCH, "Do not match"))):
        h = [((np.minimum(N[Y == v], 3) == k).mean()) for k in range(4)]; SD = (SD if j else []) + [["a", n, ["0", "1", "2", "3+"][k], "", "", int((np.minimum(N[Y == v], 3) == k).sum()), h[k]] for k in range(4)]
        a.bar(np.arange(4) + (j - 0.5) * 0.44, h, width=0.42, color=c, lw=0, zorder=2, label=f"{n} ({(Y == v).sum():,})")
        for k in range(4):
            a.text(k + (j - 0.5) * 0.44, h[k] + 0.015, f"{h[k]:.2f}" if h[k] >= 0.005 else f"{h[k]:.3f}", ha="center", va="bottom", fontsize=FS - 1.4, color=C["ink2"])
        print(n, (Y == v).sum(), "exact count 0/1/2/3+", np.round(h, 4), ">=1", (N[Y == v] >= 1).mean().round(4), "feature 26 >=1", (N26[Y == v] >= 1).mean().round(4),
              "f26 0/1/2/3+", np.round([(np.minimum(N26[Y == v], 3) == k).mean() for k in range(4)], 4), "tumor support median", np.median(T[Y == v]), "quartiles", np.percentile(T[Y == v], [25, 75]), "<=3", (T[Y == v] <= 3).mean().round(3))
    a.set_xticks(range(4)); a.set_xticklabels(["0", "1", "2", "3+"]); a.set_ylim(0, 1.08); fticks(a, "y", [0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_xlabel("Normal reads with |D| < 20 bp"); a.set_ylabel("Fraction of records"); a.legend(loc="upper right", handlelength=1.2); a.tick_params(axis="x", length=0); a.grid(axis="y", **GRID)
    tag(a, "a", 2, 2.5)
    a = ax(96, 9.5, 62, 44)
    B = np.arange(0.5, 63, 2)
    for v, c in ((1, MATCH), (0, NOMATCH)):
        a.stairs(np.histogram(np.minimum(T[Y == v], 61), B)[0] / (Y == v).sum(), B, baseline=None, color=c, lw=1.3)
        hh = np.histogram(np.minimum(T[Y == v], 61), B)[0]; SD += [["b", "Match the training labels" if v else "Do not match", f"{B[i] + 0.5:.0f}-{B[i + 1] - 0.5:.0f}" if i < len(hh) - 1 else ">=61", B[i], B[i + 1], int(hh[i]), hh[i] / (Y == v).sum()] for i in range(len(hh))]
    a.set_xlim(0, 63); a.set_ylim(0, 0.31); fticks(a, "y", [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]); a.set_xticks([1, 10, 20, 30, 40, 50, 61.5]); a.set_xticklabels(["1", "10", "20", "30", "40", "50", "≥61"])
    a.set_xlabel("Tumor reads supporting the record"); a.set_ylabel("Fraction of records"); a.grid(axis="y", **GRID)
    tag(a, "b", 84, 2.5)
    print("tumor support 1-2 / 3-4 fraction", [[round(float(((T[Y == v] >= a) & (T[Y == v] <= a + 1)).mean()), 3) for a in (1, 3)] for v in (1, 0)], "min", T.min())
    print("agreement feature 26 vs exact", (N == N26).mean().round(3), "records", len(Y))
    VAL["evidence"] = {v: (int((Y == v).sum()), float((N[Y == v] >= 1).mean())) for v in (1, 0)}
    fig.savefig(f"{OUT}/fig_evidence.pdf", facecolor="white")
    plt.close(fig)
    src("fig_evidence", ["panel", "series", "bin", "bin_low", "bin_high", "records", "fraction"], SD)


def filesfig():
    GB = 2.0 ** 30
    NAME = PAIRNAME
    bam, tot = {}, {}
    for p in P5:
        J = json.load(open(f"{H}/runs_clean/full/{p}/run.json"))
        bam[p] = [os.path.getsize(os.path.realpath(f"{ROOT}/hifi/{p}/{s}.bam")) / GB for s in "nt"]
        f = {x["file"]: x["bytes"] for x in J["files"] if not x["file"].startswith(("log.", "steps", "calls.vcf"))}
        pk = sum(v for k, v in f.items() if k in ("n.svb", "t.svb", "n.svb.idx", "t.svb.idx"))
        ex = sum(v for k, v in f.items() if k.endswith(".exc"))
        cn = sum(v for k, v in f.items() if k.endswith(".cnt"))
        gr = f.get("graph.gph", 0)
        al = sum(f.values())
        tot[p] = [x / GB for x in (pk, ex, cn, gr, al - pk - ex - cn - gr, al)]
        print(NAME[p], "files", ["%.1f" % x for x in tot[p]], "BAM n+t (2^30) %.1f" % sum(bam[p]), "BAM n+t (1e9) %.1f" % (sum(bam[p]) * GB / 1e9), "ratio %.3f" % (tot[p][5] / sum(bam[p])))


    def grab(p, step, pat):
        m = re.search(pat, open(f"{H}/runs_clean/full/{p}/logs/log.{step}").read())
        return [int(x) for x in m.groups()]


    ROWS = [("Sites", "prune", r"dropped nodes, (\d+) sites"),
            ("Reads cut", "select", r"on (\d+) reads"), ("Snippets", "select", r"reads; (\d+) pieces"),
            ("Snippets placed", "place", r"pieces, (\d+) placed"), ("Sites with an event", "call", r"\[call\] (\d+) of \d+ sites carry"),
            ("Records", "call", r"at sites, (\d+) records")]
    stat = {lab: [grab(p, step, pat)[0] for p in P5] for lab, step, pat in ROWS}
    stat["PASS calls"] = [sum(1 for l in open(f"{H}/verdict/{p}.vcf") if "\tPASS\t" in l) for p in P5]
    for k, v in stat.items():
        print(k, v)
    VAL["files"] = (tot, bam, stat)
    KIND = [("Packed reads", STAGE["Packing"]), ("Core positions", STAGE["Core extraction"]), ("Counts", STAGE["Counting"]), ("Graph", STAGE["Graph and pruning"]), ("Snippets, alignments, other", STAGE["Placement"])]
    fig, ax = page(W, 76)
    a = ax(14, 19, 62, 40)
    for i, p in enumerate(P5):
        b = 0
        for j, (lab, c) in enumerate(KIND):
            a.bar(i - 0.2, tot[p][j], bottom=b, width=0.36, color=c, lw=0, zorder=2, label=lab if i == 0 else None)
            b += tot[p][j]
        a.text(i - 0.2, b + 2, "%.0f" % b, ha="center", va="bottom", fontsize=FS - 1.4, color=C["ink2"])
        a.bar(i + 0.2, sum(bam[p]), width=0.36, color=LIGHT, lw=0, zorder=2, label="Two sorted BAM files" if i == 0 else None)
        a.text(i + 0.2, sum(bam[p]) + 2, "%.0f" % sum(bam[p]), ha="center", va="bottom", fontsize=FS - 1.4, color=C["ink2"])
    a.set_xticks(range(5)); a.set_xticklabels([NAME[p] for p in P5], rotation=35, ha="right", rotation_mode="anchor"); a.set_ylabel("Size (GB)"); a.set_ylim(0, 185)
    a.tick_params(axis="x", length=0); a.grid(axis="y", **GRID)
    a.legend(ncol=2, loc="lower left", bbox_to_anchor=(-0.02, 1.0), handlelength=1.1, columnspacing=1.0, handletextpad=0.4)
    tag(a, "a", 2, 11.6)
    a = ax(96, 19, 62, 40)
    L = list(stat)
    for i, p in enumerate(P5):
        a.plot(range(len(L)), [stat[l][i] for l in L], marker="o", ms=2.8, lw=1.0, color=PAIRCOL[p], mec="white", mew=0.3, label=NAME[p])
    a.set_yscale("log"); a.set_ylim(100, 1e6); a.set_xticks(range(len(L))); a.set_xticklabels(L, rotation=35, ha="right", rotation_mode="anchor"); a.set_ylabel("Number of objects")
    a.grid(axis="y", **GRID)
    a.legend(ncol=3, loc="lower left", bbox_to_anchor=(-0.02, 1.0), handlelength=1.4, columnspacing=1.0, handletextpad=0.4)
    tag(a, "b", 84, 11.6)
    fig.savefig(f"{OUT}/fig_files.pdf", facecolor="white")
    plt.close(fig)


def jointfig():
    from matplotlib.patches import Rectangle
    S = json.load(open(f"{H}/runs/sim/scores.json"))
    fig, ax = page(W, 60.5)
    for n, (name, lab, x0) in enumerate((("20-20", "Normal 20×, tumors 20×", 16), ("20-5", "Normal 20×, tumors 5×", 62))):
        w = S[name]["sweep"]; ra, fl = w["raise_"], w["floor"]; g = np.array(w["PRF"])[:, :, 2]
        a = ax(x0, 12.5, 40, 38)
        lo, hi = g.min() - 0.002, g.max() + 0.002
        a.imshow(g, cmap=HEAT, vmin=lo, vmax=hi, aspect="auto", origin="lower")
        for i in range(len(ra)):
            for j in range(len(fl)):
                a.text(j, i, f"{g[i, j]:.3f}", ha="center", va="center", fontsize=FS - 0.5, color=ink(HEAT((g[i, j] - lo) / (hi - lo)), HEAT), fontweight="bold" if (ra[i], fl[j]) == (0.2, 0.1) else "normal")
        a.add_patch(Rectangle((fl.index(0.1) - 0.5, ra.index(0.2) - 0.5), 1, 1, fill=False, ec=JANUS, lw=1.4))
        a.set_xticks(range(len(fl))); a.set_xticklabels([f"{v:g}" for v in fl]); a.set_yticks(range(len(ra))); a.set_yticklabels([f"+{v:g}" for v in ra])
        a.set_xlabel("Lowest probability accepted\nwith support from another tumor"); a.set_ylabel("Raise of the threshold\nwithout such support") if n == 0 else None
        for s in a.spines.values(): s.set_visible(False)
        a.tick_params(length=0)
        title(a, lab, 1.03)
        tag(a, "ab"[n], 2 if n == 0 else 54, 5.3)
        VAL.setdefault("joint", {})[name] = (float(g.min()), float(g.max()))
        SD = (SD if n else []) + [["ab"[n], name.replace("-", "/"), ra[i], fl[j], "", *w["PRF"][i][j], int((ra[i], fl[j]) == (0.2, 0.1))] for i in range(len(ra)) for j in range(len(fl))]
        print(name, "F1 range %.3f-%.3f" % (g.min(), g.max()))
    c = ax(118, 12.5, 40, 38); MUC = ["#6a3d9a", "#e0a21a", "#8c5a2b", "#b39ddb", "#56616c"]
    for name, col in zip(SIMSETS, MUC):
        m = S[name]["mu"]; c.plot(m["mu"], [x[2] for x in m["PRF"]], marker="o", ms=2.8, lw=1.0, color=col, mec="white", mew=0.3, label=name.replace("-", "/"))
        SD += [["c", name.replace("-", "/"), "", "", mu_, *prf, int(mu_ == 0.5)] for mu_, prf in zip(m["mu"], m["PRF"])]
    c.axvline(0.5, color="#bbbbbb", lw=0.6, zorder=0); c.set_xlabel("Penalty μ"); c.set_ylabel("Median F1 score of the clade pass")
    c.set_xticks([0.3, 0.5, 1.0, 1.5, 2.0]); c.set_xticklabels(["0.3", "0.5", "1", "1.5", "2"]); c.grid(axis="y", **GRID)
    c.legend(fontsize=FS - 1.2, handlelength=1.2, labelspacing=0.25, frameon=False, title="Depth n/t", title_fontsize=FS - 1.2, loc="center right"); tag(c, "c", 108, 5.3)
    fig.savefig(f"{OUT}/fig_joint.pdf", facecolor="white")
    plt.close(fig)
    src("fig_joint", ["panel", "depth_normal_tumor", "threshold_raise", "lowest_probability_with_support", "mu", "precision", "recall", "f1", "rule_used"], SD)


def sweepfig():
    from matplotlib.patches import Rectangle
    V = {(int(r["theta0"]), int(r["theta_t"])): r for r in csv.DictReader(open(f"{R}/sweep/rescore.txt"), delimiter="\t")}
    N0, NT = [5, 10, 15, 20, 25], [20, 30, 40, 50, 60]
    PAN = (("Full depth", ["with"]), ("Depth 30× to 5×, mean", ["30x", "20x", "15x", "10x", "5x"]), ("Purity 100% to 20%, mean", ["100p", "80p", "60p", "40p", "20p"]))
    fig, ax = page(W, 52)
    for n, (lab, ks) in enumerate(PAN):
        g = np.array([[st.mean(float(V[(a, b)][k]) for k in ks) for b in NT] for a in N0])
        SD = (SD if n else []) + [["abc"[n], lab, N0[i], NT[j], g[i, j], int((N0[i], NT[j]) == (20, 20))] for i in range(len(N0)) for j in range(len(NT))]
        x0 = 14 + n * 51; a = ax(x0, 10, 34, 32)
        lo, hi = g.min() - 0.002, g.max() + 0.002
        im = a.imshow(g, cmap=HEAT, vmin=lo, vmax=hi, aspect="auto", origin="lower")
        for i in range(len(N0)):
            for j in range(len(NT)):
                a.text(j, i, f"{g[i, j]:.3f}", ha="center", va="center", fontsize=FS - 1.4, color=ink(HEAT((g[i, j] - lo) / (hi - lo)), HEAT), fontweight="bold" if (N0[i], NT[j]) == (20, 20) else "normal")
        a.add_patch(Rectangle((NT.index(20) - 0.5, N0.index(20) - 0.5), 1, 1, fill=False, ec=JANUS, lw=1.4))
        a.set_xticks(range(len(NT))); a.set_xticklabels(NT); a.set_yticks(range(len(N0))); a.set_yticklabels(N0)
        a.set_xlabel("$\\tau_{\\mathrm{tumor}}$")
        if n == 0: a.set_ylabel("$\\tau_{\\mathrm{normal}}$")
        for sp in a.spines.values(): sp.set_visible(False)
        a.tick_params(length=0)
        title(a, lab, 1.03, 0.0)
        cb = fig.colorbar(im, cax=ax(x0 + 35.5, 10, 1.8, 32)); cb.ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4)); cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(ffmt))
        cb.ax.tick_params(length=1.5, pad=1, labelsize=FS - 1.4); cb.outline.set_linewidth(0.4); cb.ax.set_xlabel("F1", fontsize=FS - 1.4, labelpad=2)
        tag(a, "abc"[n], 2 if n == 0 else x0 - 7, 5.3)
        print("sweep", lab, "range %.4f-%.4f, at 20/20 %.4f" % (g.min(), g.max(), g[N0.index(20), NT.index(20)]))
    fig.savefig(f"{OUT}/fig_sweep.pdf", facecolor="white")
    plt.close(fig)
    src("fig_sweep", ["panel", "settings", "tau_normal", "tau_tumor", "f1", "thresholds_used"], SD)


def aucfig():
    from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, average_precision_score
    F = json.load(open(f"{H}/verdict/folds.json"))
    fig, ax = page(W, 64)
    a, b = ax(16, 10, 62, 44), ax(96, 10, 62, 44)
    for p in P5:
        pr, y, d = np.array(F[p]["prob"]), np.array(F[p]["label"]), np.array(F[p]["dupl"], bool)
        pr, y = pr[~d], y[~d]
        x1, y1, t1 = roc_curve(y, pr); SD = (SD if p != P5[0] else []) + [["a", PAIRNAME[p], *z, roc_auc_score(y, pr)] for z in zip(x1, y1, t1)]; a.plot(x1, y1, color=PAIRCOL[p], lw=1.0, label=f"{PAIRNAME[p]} ({roc_auc_score(y, pr):.3f})")
        pp, rr, t2 = precision_recall_curve(y, pr); SD += [["b", PAIRNAME[p], *z, average_precision_score(y, pr)] for z in zip(rr, pp, list(t2) + [""])]; b.plot(rr, pp, color=PAIRCOL[p], lw=1.0, label=f"{PAIRNAME[p]} ({average_precision_score(y, pr):.3f})")
    a.plot([0, 1], [0, 1], color="#bbbbbb", lw=0.6, ls=(0, (3, 2)))
    a.set_xlabel("False positive rate"); a.set_ylabel("True positive rate"); b.set_xlabel("Recall"); b.set_ylabel("Precision")
    for x, lab, loc, xl, yl in ((a, "ROC-AUC", "lower right", (0, 0.2), (0.75, 1.003)), (b, "Average precision", "lower left", (0.6, 1.003), (0.75, 1.003))):
        x.set_xlim(*xl); x.set_ylim(*yl); x.grid(**GRID); x.legend(title=lab, loc=loc, fontsize=FS - 1, title_fontsize=FS - 1, handlelength=1.4)
        fticks(x, "x", np.linspace(xl[0], round(xl[1], 1), 5)); fticks(x, "y", np.linspace(0.75, 1, 6))
    title(a, "ROC curve", 1.03); title(b, "Precision\u2013recall curve", 1.03); tag(a, "a", 2, 5.3); tag(b, "b", 82, 5.3)
    fig.savefig(f"{OUT}/fig_auc.pdf", facecolor="white")
    plt.close(fig)
    src("fig_auc", ["panel", "pair", "x", "y", "threshold", "area"], SD)


def splitsfig():
    SP = supp_tsv("splits_3_2.tsv"); A = supp_tsv("heldout_auc.tsv")
    fig, ax = page(W, 57)
    a, b = ax(16, 10.5, 48, 40), ax(84, 10.5, 74, 40)
    for i, p in enumerate(P5):
        v = [float(r["F1"]) for r in SP if r["pair"] == p and r["setting"] == "full"]
        f4 = next(float(r["F1"]) for r in A if r["pair"] == p and r["setting"] == "full")
        a.scatter([i] * len(v), v, s=10, color="#9a9a9a", zorder=2, label="three training pairs" if i == 0 else None)
        a.scatter([i], [f4], s=26, marker="D", facecolor="none", edgecolor=JANUS, lw=1.0, zorder=3, label="four training pairs" if i == 0 else None)
    a.set_xticks(range(5)); a.set_xticklabels([PAIRNAME[p] for p in P5], rotation=30, ha="right"); a.set_xlim(-0.6, 4.6)
    a.set_ylabel("F1 score at full depth"); a.grid(axis="y", **GRID); a.legend(loc="lower left", fontsize=FS - 1, handletextpad=0.3); a.yaxis.set_major_formatter(plt.FuncFormatter(ffmt))
    ks = SETTING + PURITY; sp = collections.defaultdict(list)
    for r in SP: sp[(r["train"], r["setting"])].append(float(r["F1"]))
    m3 = {k: [st.mean(v) for (tr, s_), v in sp.items() if s_ == k] for k in ks}
    four = {k: st.mean(float(r["F1"]) for r in A if r["setting"] == k) for k in ks}
    for seg, first in ((range(0, 6), True), (range(6, 11), False)):
        kk = [ks[i] for i in seg]; seg = [i + (i > 5) for i in seg]
        b.fill_between(list(seg), [min(m3[k]) for k in kk], [max(m3[k]) for k in kk], color="#e0e0e0", lw=0, label="three training pairs, range" if first else None)
        b.plot(list(seg), [st.mean(m3[k]) for k in kk], "-o", color="#7f7f7f", ms=2.4, lw=0.9, label="three training pairs, mean" if first else None)
        b.plot(list(seg), [four[k] for k in kk], "-D", color=JANUS, ms=2.4, lw=1.0, label="four training pairs" if first else None)
    b.axvline(6, color="#bbbbbb", lw=0.6)
    b.set_xticks([0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11]); b.set_xticklabels(["full", "30\u00d7", "20\u00d7", "15\u00d7", "10\u00d7", "5\u00d7", "100", "80", "60", "40", "20"])
    b.set_xlabel("Depth (left) and tumor purity (%) at 30\u00d7 (right)"); b.set_ylabel("Mean F1 score, held-out pairs"); b.grid(axis="y", **GRID); b.legend(loc="lower right", fontsize=FS - 1)
    b.yaxis.set_major_formatter(plt.FuncFormatter(ffmt)); tag(a, "a", 2, 2.5); tag(b, "b", 72, 2.5)
    fig.savefig(f"{OUT}/fig_splits.pdf", facecolor="white")
    plt.close(fig)
    LB = dict(zip(ks, ["full", "30x", "20x", "15x", "10x", "5x", "100%", "80%", "60%", "40%", "20%"])); PRF = ("P", "R", "F1"); SD = []
    SD += [["a", "full", "three training pairs", r["train"], PAIRNAME[p], *(float(r[x]) for x in PRF), "", ""] for p in P5 for r in SP if r["pair"] == p and r["setting"] == "full"]
    SD += [["a", "full", "four training pairs", r["train"], PAIRNAME[p], *(float(r[x]) for x in PRF), "", ""] for p in P5 for r in A if r["pair"] == p and r["setting"] == "full"]
    for k in ks:
        tr = sorted({r["train"] for r in SP if r["setting"] == k})
        SD += [["b", LB[k], "three training pairs", t, "mean of held-out pairs", *(st.mean(float(r[x]) for r in SP if r["train"] == t and r["setting"] == k) for x in PRF), "", ""] for t in tr]
        SD += [["b", LB[k], "three training pairs, mean", "", "", *(st.mean(st.mean(float(r[x]) for r in SP if r["train"] == t and r["setting"] == k) for t in tr) for x in PRF), min(m3[k]), max(m3[k])]]
        SD += [["b", LB[k], "four training pairs", "", "mean of held-out pairs", *(st.mean(float(r[x]) for r in A if r["setting"] == k) for x in PRF), "", ""]]
    src("fig_splits", ["panel", "setting", "series", "training_pairs", "held_out_pair", "precision", "recall", "f1", "f1_range_low", "f1_range_high"], SD)


def paramfig():
    GK = ["full", "20x", "5x", "20p", "all"]; GS = [s_ for s_ in SETTING + PURITY if s_ != "100p"]
    A = supp_tsv("heldout_auc.tsv")
    am = lambda k: st.mean(float(r["F1"]) for r in A if r["setting"] == k)
    used = {k: st.mean(am(s_) for s_ in GS) if k == "all" else am(k) for k in GK}
    fig, ax = page(W, 52)
    rg = np.random.default_rng(0)
    for n, (j, lab) in enumerate((("grid", "scikit-learn, 243 sets"),)):
        G = collections.defaultdict(dict)
        for r in supp_tsv(f"{j}.tsv"): G[r["variant"]][(r["pair"], r["setting"])] = float(r["F1"])
        GP = collections.defaultdict(dict)
        for r in supp_tsv(f"{j}.tsv"): GP[r["variant"]][(r["pair"], r["setting"])] = (float(r["P"]), float(r["R"]))
        gm = lambda d, k, x: st.mean(d[(p, s_)][x] for p in P5 for s_ in GS if k in ("all", s_))
        um = lambda k, x: st.mean(st.mean(float(r[x]) for r in A if r["setting"] == s_) for s_ in GS if k in ("all", s_))
        SD = [[v, k, gm(GP[v], k, 0), gm(GP[v], k, 1), st.mean(G[v][(p, s_)] for p in P5 for s_ in GS if k in ("all", s_)), um(k, "P"), um(k, "R"), used[k], st.mean(G[v][(p, s_)] for p in P5 for s_ in GS if k in ("all", s_)) - used[k]] for v in G for k in GK]
        a = ax(16, 6, 142, 44); lo_, hi_ = 0, 0
        for i, k in enumerate(GK):
            v = np.array([st.mean(d[(p, s_)] for p in P5 for s_ in GS if k in ("all", s_)) for d in G.values()]) - used[k]; lo_, hi_ = min(lo_, v.min()), max(hi_, v.max())
            a.scatter(i + rg.uniform(-0.28, 0.28, len(v)), v, s=2.5, color="#8a8a8a", lw=0, zorder=2)
        a.set_ylim(lo_ - 0.006, hi_ + 0.006); a.yaxis.set_major_formatter(plt.FuncFormatter(ffmt))
        a.axhline(0, color=JANUS, lw=1.0, zorder=3, label="parameters used in this paper")
        a.set_xticks(range(len(GK))); a.set_xticklabels(["Full", "20\u00d7", "5\u00d7", "20% purity", "All"]); a.grid(axis="y", **GRID)
        a.set_ylabel("F1 score minus that of\nthe parameters used") if n == 0 else None
        a.legend(loc="upper right", fontsize=FS - 1)
    fig.savefig(f"{OUT}/fig_params.pdf", facecolor="white")
    plt.close(fig)
    src("fig_params", ["variant", "settings", "precision", "recall", "f1", "used_precision", "used_recall", "used_f1", "f1_minus_used"], SD)


def specfig():
    fig, ax = page(W, 59.5)
    pan = (("tumor_only", "Tumor-only cores"), ("shared_tumor", "Cores shared with the normal sample"))
    for n, (key, lab) in enumerate(pan):
        a = ax(16 + n * 80, 10, 60, 40)
        for p in P5:
            h = {l.split()[0]: np.array(l.split()[1:], float) for l in open(f"{R}/full/{p}/t.cnt.hist")}[key]
            x = np.arange(2, 121); a.plot(x, np.where(h[2:121] > 0, h[2:121], np.nan), color=PAIRCOL[p], lw=0.9, label=PAIRNAME[p]); SD = (SD if n or p != P5[0] else []) + [["ab"[n], lab, PAIRNAME[p], i, h[i]] for i in x]
        a.set_yscale("log"); a.set_xlim(0, 120); a.grid(axis="y", **GRID)
        a.set_xlabel("Number of tumor reads with the core"); a.set_ylabel("Cores")
        a.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10))); a.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        title(a, lab, 1.03); tag(a, "ab"[n], 2 + n * 80, 5.3)
        if n == 0: a.legend(loc="upper right", fontsize=FS - 1, handlelength=1.4)
    fig.savefig(f"{OUT}/fig_spectra.pdf", facecolor="white")
    plt.close(fig)
    src("fig_spectra", ["panel", "cores", "pair", "tumor_reads", "count"], SD)



def datasets():
    P = P5 + ["hg008", "colo829"]
    R = {p: json.load(open(f"{H}/runs/full/{p}/run.json")) for p in P}
    X = {p: R[p]["params"]["bases_per_x"] for p in P}
    dep = {p: [R[p]["reads"][s]["bases"] / X[p] for s in "nt"] for p in P}
    ln = {p: [R[p]["reads"][s]["bases"] / R[p]["reads"][s]["reads"] / 1e3 for s in "nt"] for p in P}
    fig, ax = page(160, 64.5)
    for k, (D, yl, fmt, top) in enumerate(((dep, "Sequencing depth (×)", "%.0f", 90), (ln, "Mean read length (kb)", "%.1f", 25))):
        a = ax(13 + 81 * k, 18, 64, 40)
        for j, al in enumerate((1.0, 0.45)):
            xs = [i + (j - 0.5) * 0.4 for i in range(len(P))]
            a.bar(xs, [D[p][j] for p in P], width=0.38, color=[PAIRCOL[p] for p in P], alpha=al, lw=0)
            for x, p in zip(xs, P):
                a.text(x, D[p][j] + top * 0.015, fmt % D[p][j], ha="center", va="bottom", rotation=90, fontsize=FS - 1.2, bbox=dict(fc="white", ec="none", pad=0.2))
        a.set_xticks(range(len(P))); a.set_xticklabels([PAIRNAME[p] for p in P], rotation=35, ha="right", rotation_mode="anchor")
        a.set_xlim(-0.6, len(P) - 0.4); a.set_ylim(0, top); a.set_ylabel(yl); a.grid(axis="y", **GRID); a.set_axisbelow(True)
        tag(a, "ab"[k], 2 + 81 * k, 2.5)
    fig.legend(handles=[Patch(color="#4d4d4d", lw=0, label="Normal"), Patch(color="#4d4d4d", alpha=0.45, lw=0, label="Tumor")], loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, handlelength=1.2, columnspacing=1.6, borderaxespad=0.3)
    save(fig, "datasets")
    src("fig_datasets", ["panel", "pair", "sample", "value", "unit"], [[pn, PAIRNAME[p], sm, D[p][j], u] for pn, D, u in (("a", dep, "x"), ("b", ln, "kb")) for p in P for j, sm in enumerate(("normal", "tumor"))])
    for p in P:
        print(p, X[p], "depth n/t %.1f/%.1f" % tuple(dep[p]), "len n/t %.2f/%.2f kb" % tuple(ln[p]))


def consensus():
    NS = [10, 1, 4, 3]
    TN = ["HiFi", "ONT", "Illumina"]
    hist, combo = {}, {}
    for p in P5:
        r = csv.reader(open(f"{H}/eval_sets/{p}/all_with/None_support.tsv"), delimiter="\t")
        h = next(r); e = h.index("ensemble"); a = h.index("Minda_IDs") + 1
        assert len(h[a:a + 18]) == 18 and h[a + 18].startswith("condition"), h[a:]
        hist[p] = collections.Counter(); combo[p] = collections.Counter(); n = 0
        for f in r:
            if f[e] != "True":
                continue
            v = [x == "True" for x in f[a:a + 18]]
            assert all(x in ("True", "False") for x in f[a:a + 18])
            hi = any(v[:11]); on = any(v[11:15]); il = any(v[15:18])
            n += 1; hist[p][sum(v)] += 1; combo[p][(hi, on, il)] += 1
        print(p, "consensus", n, "min sets", min(hist[p]), "max", max(hist[p]))
        print("  hist", sorted(hist[p].items())); print("  combo", dict(combo[p]))
    fig, ax = page(160, 63.5)
    a = ax(13, 13, 64, 44)
    K = list(range(6, 19))
    for p in P5:
        a.plot(K, [hist[p][k] for k in K], marker="o", ms=2.6, lw=0.9, color=PAIRCOL[p], label=PAIRNAME[p])
    a.set_xticks(K); a.set_ylim(0, 260); a.set_xlim(5.5, 18.5); a.set_xlabel("Number of supporting call sets"); a.set_ylabel("Consensus SVs"); a.grid(axis="y", **GRID)
    a.legend(loc="upper left", handlelength=1.4)
    tag(a, "a", 2, 2.5)
    a = ax(93, 13, 65, 44)
    CB = [((True, True, True), "HiFi, ONT\nand Illumina"), ((True, True, False), "HiFi and\nONT only"), ((True, False, True), "HiFi and\nIllumina only"), ((False, True, True), "ONT and\nIllumina only")]
    top = max(max(c.values()) for c in combo.values()) * 1.2
    for j, p in enumerate(P5):
        xs = [i + (j - 2) * 0.17 for i in range(4)]
        a.bar(xs, [combo[p][k] for k, _ in CB], width=0.16, color=PAIRCOL[p], lw=0)
        for x, (k, _) in zip(xs, CB):
            a.text(x, combo[p][k] + top * 0.012, str(combo[p][k]), ha="center", va="bottom", fontsize=FS - 1.4, rotation=90, bbox=dict(fc="white", ec="none", pad=0.2))
    a.set_xticks(range(4)); a.set_xticklabels([l for _, l in CB]); a.set_ylabel("Consensus SVs"); a.set_xlabel("Technologies with a supporting call set")
    a.set_xlim(-0.5, 3.5); a.set_ylim(0, top); a.grid(axis="y", **GRID); a.set_axisbelow(True); a.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:,.0f}"))
    tag(a, "b", 82, 2.5)
    save(fig, "consensus")
    src("fig_consensus", ["panel", "pair", "item", "consensus_svs"], [["a", PAIRNAME[p], k, hist[p][k]] for p in P5 for k in K] + [["b", PAIRNAME[p], l.replace("\n", " "), combo[p][k]] for p in P5 for k, l in CB])
    tot = collections.Counter(); th = collections.Counter()
    for p in P5:
        tot.update(combo[p]); th.update(hist[p])
    print("pooled combo", {l.replace("\n", " "): tot[k] for k, l in CB}, "total", sum(tot.values()))
    print("pooled hist", sorted(th.items()), "all 18:", th[18], "exactly 6:", th[6])


def simgen():
    R = f"{H}/runs/sim/truth"
    TY = ["DEL", "INS", "DUP", "INV", "dDUP", "rTRA"]
    COL = {**SVTYPE, "dDUP": SVTYPE["DUP"], "rTRA": SVTYPE["BND"]}
    AL = {"dDUP": 0.45}
    LS = {"dDUP": (0, (2.2, 1.2))}
    T = ["T%02d" % i for i in range(1, 21)]
    per, sv, car = {}, {}, collections.Counter()
    for t in T:
        s = {}
        for l in open(f"{R}/{t}.tsv"):
            f = l.split("\t"); s[f[0]] = (f[1], int(f[2]))
        s = {m: v for m, v in s.items() if v[0] in ("dDUP", "rTRA") or v[1] >= 50}
        per[t] = collections.Counter(v[0] for v in s.values()); sv.update(s)
        for m in s: car[m] += 1
    fig, ax = page(160, 65)
    a = ax(14, 16, 44, 42)
    b = np.zeros(20)
    for ty in TY:
        v = np.array([per[t][ty] for t in T]); a.bar(range(20), v, bottom=b, color=COL[ty], alpha=AL.get(ty, 1), width=0.8, lw=0); b += v
    a.set_xticks(range(20)); a.set_xticklabels([t[1:] for t in T], rotation=90); a.set_xlim(-0.7, 19.7); a.set_ylim(0, 1500)
    a.set_xlabel("Tumor"); a.set_ylabel("Simulated SVs"); a.grid(axis="y", **GRID); a.set_axisbelow(True); a.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:,.0f}"))
    tag(a, "a", 2, 2.5)
    a = ax(71, 16, 37, 42)
    bins = np.logspace(np.log10(50), np.log10(3e6), 40)
    for ty in TY:
        a.hist([v[1] for v in sv.values() if v[0] == ty], bins=bins, histtype="step", color=COL[ty], ls=LS.get(ty, "-"), lw=0.9)
    a.set_xscale("log"); a.set_xlabel("SV size (bp)"); a.set_ylabel("Distinct SVs"); a.set_xlim(30, 6e6); a.set_xticks([1e2, 1e4, 1e6]); a.set_ylim(0, 155); a.grid(axis="y", **GRID)
    tag(a, "b", 61, 2.5)
    a = ax(121, 16, 37, 42)
    h = collections.Counter(car.values())
    a.plot([i for i in range(1, 21) if h[i]], [h[i] for i in range(1, 21) if h[i]], "o", color=NEUTRAL, ms=3.2, mew=0, zorder=3); a.grid(axis="y", **GRID)
    a.set_yscale("log"); a.set_ylim(8, 30000); a.set_xlim(0.2, 20.8); a.set_xticks([1, 5, 10, 15, 20]); a.set_xlabel("Tumors carrying the SV"); a.set_ylabel("Distinct SVs")
    for i in range(1, 21):
        if h[i]: a.text(i, h[i] * 1.35, f"{h[i]:,}", ha="center", fontsize=FS - 1.4, rotation=90, va="bottom")
    tag(a, "c", 110, 2.5)
    fig.legend(handles=[(Patch(facecolor=COL[t], alpha=AL.get(t, 1), lw=0), Line2D([], [], color=COL[t], ls=LS.get(t, "-"), lw=0.9)) for t in TY], labels=TY, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)}, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=6, handlelength=3.0, columnspacing=2.0, handletextpad=0.5, borderaxespad=0.3)
    save(fig, "simgen")
    SD = [["a", ty, t, "", "", per[t][ty]] for t in T for ty in TY]
    for ty in TY:
        hh = np.histogram([v[1] for v in sv.values() if v[0] == ty], bins=bins)[0]; SD += [["b", ty, "", bins[i], bins[i + 1], int(hh[i])] for i in range(len(hh))]
    SD += [["c", "", i, "", "", h[i]] for i in range(1, 21)]
    src("fig_simgen", ["panel", "svtype", "item", "size_low_bp", "size_high_bp", "svs"], SD)
    n = [sum(per[t].values()) for t in T]
    print("distinct", len(sv), dict(collections.Counter(v[0] for v in sv.values())), "per tumor", min(n), max(n), "carriers", sorted(h.items()))
    VAL["simgen"] = dict(n=len(sv), types=dict(collections.Counter(v[0] for v in sv.values())), per=n)


def align():
    F = f"{RC}/full"
    D = {}
    for p in P5:
        L = lambda s: open(f"{F}/{p}/logs/log.{s}").read()
        reads, pieces, bp = map(int, re.search(r"on (\d+) reads; (\d+) pieces, (\d+) bp", L("select")).groups())
        placed = int(re.search(r"(\d+) placed by cores", L("place")).group(1))
        win = int(re.search(r"for (\d+) target", L("local")).group(1))
        loc = sum(int(x) for x in re.findall(r"mapped (\d+) sequences", L("local")))
        gaps = sum(int(x) for x in re.findall(r"mapped (\d+) sequences", L("gaps")))
        r = json.load(open(f"{R}/full/{p}/run.json"))["reads"]["t"]
        D[p] = dict(reads=reads, pieces=pieces, bp=bp, placed=placed, windows=win, local=loc, gaps=gaps, tbases=r["bases"], treads=r["reads"])
        assert loc == placed
    fig, ax = page(160, 62)
    X = np.arange(5); N = [PAIRNAME[p] for p in P5]
    def fin(a, lab, tg):
        a.set_xticks(X); a.set_xticklabels(N, rotation=30, ha="right", rotation_mode="anchor"); a.set_xlim(-0.6, 4.6); a.set_ylabel(lab); a.grid(axis="y", **GRID); a.set_axisbelow(True); tag(a, tg, a.get_position().x0 * 160 - 14, 2)
    a = ax(15, 15, 38, 40)
    a.bar(X - 0.2, [D[p]["pieces"] / 1e3 for p in P5], width=0.38, color=STAGE["Selection"], label="Selected")
    a.bar(X + 0.2, [D[p]["placed"] / 1e3 for p in P5], width=0.38, color=STAGE["Placement"], label="Placed by cores")
    for i, p in enumerate(P5): a.text(i + 0.2, D[p]["placed"] / 1e3 + 1.5, f"{100 * D[p]['placed'] / D[p]['pieces']:.0f}%", ha="center", va="bottom", rotation=90, fontsize=FS - 1.2, bbox=dict(fc="white", ec="none", pad=0.2))
    a.set_ylim(0, 120); a.legend(loc="upper right", handlelength=1.2); fin(a, "Snippets (thousands)", "a")
    a = ax(67, 15, 38, 40)
    v = [100 * D[p]["bp"] / D[p]["tbases"] for p in P5]
    a.bar(X, v, width=0.6, color=STAGE["Selection"])
    for i, p in enumerate(P5): a.text(i, v[i] + 0.015, f"{D[p]['bp'] / 1e6:,.0f} Mbp", ha="center", va="bottom", rotation=90, fontsize=FS - 1.2, bbox=dict(fc="white", ec="none", pad=0.2))
    a.set_ylim(0, 0.8); fin(a, "Share of tumor read bases (%)", "b")
    a = ax(121, 15, 38, 40)
    v = [D[p]["gaps"] for p in P5]
    a.bar(X, v, width=0.6, color=STAGE["minimap2"])
    for i, p in enumerate(P5): a.text(i, v[i] + 150, f"{D[p]['gaps']:,}", ha="center", va="bottom", rotation=90, fontsize=FS - 1.2, bbox=dict(fc="white", ec="none", pad=0.2))
    a.set_ylim(0, 10500); a.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:,.0f}")); fin(a, "Parts realigned", "c")
    save(fig, "align")
    for p in P5: print(p, D[p], "placed %.3f share %.4f%% mean piece %.0f bp" % (D[p]["placed"] / D[p]["pieces"], 100 * D[p]["bp"] / D[p]["tbases"], D[p]["bp"] / D[p]["pieces"]))


def signatures():
    N = collections.Counter()
    for p in P5:
        for l in open(f"{H}/verdict/{p}.vcf"):
            if l[0] == "#": continue
            f = l.split("\t"); info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            N[(p, f[6], info["SVTYPE"])] += 1
    TY = [t for t in ["DEL", "INS", "DUP", "INV", "BND"] if any(k[2] == t for k in N)] + sorted({k[2] for k in N} - {"DEL", "INS", "DUP", "INV", "BND"})
    assert set(TY) <= set(SVTYPE), TY
    fig, ax = page(160, 62)
    a = ax(14, 10.5, 64, 42)
    b = np.zeros(5)
    for t in TY:
        v = np.array([N[(p, "PASS", t)] for p in P5]); a.bar(range(5), v, bottom=b, color=SVTYPE[t], width=0.6, label=t, lw=0); b += v
    for i in range(5): a.text(i, b[i] + 20, f"{int(b[i]):,}", ha="center", va="bottom", fontsize=FS - 1.0, bbox=dict(fc="white", ec="none", pad=0.2))
    a.set_xticks(range(5)); a.set_xticklabels([PAIRNAME[p] for p in P5], rotation=30, ha="right", rotation_mode="anchor"); a.set_xlim(-0.55, 4.55); a.set_ylim(0, 1500); a.set_ylabel("PASS records"); a.grid(axis="y", **GRID); a.set_axisbelow(True); a.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:,.0f}"))
    a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=len(TY), handlelength=1.0, columnspacing=1.2, handletextpad=0.4)
    tag(a, "a", 2, 5.5)
    a = ax(94, 10.5, 64, 42)
    for j, fl in enumerate(("PASS", "LOW", "DUPL")):
        v = [sum(N[(p, fl, t)] for p in P5) for t in TY]
        a.bar(np.arange(len(TY)) + (j - 1) * 0.27, v, width=0.25, color=FILTER[fl], label=fl, lw=0)
        for i, x in enumerate(v): a.text(i + (j - 1) * 0.27, x + 80, f"{x:,}", ha="center", va="bottom", rotation=90, fontsize=FS - 1.2, bbox=dict(fc="white", ec="none", pad=0.2))
    a.set_xticks(range(len(TY))); a.set_xticklabels(TY); a.set_xlim(-0.55, len(TY) - 0.45); a.set_ylim(0, 6500); a.set_ylabel("Records, five pairs pooled"); a.grid(axis="y", **GRID); a.set_axisbelow(True)
    a.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f"{y:,.0f}"))
    a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handlelength=1.0, columnspacing=1.6, handletextpad=0.4)
    tag(a, "b", 82, 5.5)
    save(fig, "signatures")
    src("fig_signatures", ["panel", "pair", "filter", "svtype", "records"], [["a", PAIRNAME[p], "PASS", t, N[(p, "PASS", t)]] for p in P5 for t in TY] + [["a", PAIRNAME[p], "PASS", "Total", sum(N[(p, "PASS", t)] for t in TY)] for p in P5]
        + [["b", "five pairs", fl, t, sum(N[(p, fl, t)] for p in P5)] for fl in ("PASS", "LOW", "DUPL") for t in TY])
    for p in P5: print(p, {fl: sum(v for k, v in N.items() if k[0] == p and k[1] == fl) for fl in ("PASS", "LOW", "DUPL")}, {t: N[(p, "PASS", t)] for t in TY})
    print({fl: {t: sum(N[(p, fl, t)] for p in P5) for t in TY} for fl in ("PASS", "LOW", "DUPL")}, set(k[1] for k in N))


def simeval():
    g = {"__name__": "x", "SCRIPTS": paper.SCRIPTS, "run_score": paper.run_score}
    exec(paper.SCRIPTS["simjoint"].split("out = {}")[0], g)
    T, P, TR, G, ns, load, prob, features, cutfor = (g[k] for k in ("T", "P", "TR", "G", "ns", "load", "prob", "features", "cutfor"))
    V = f"{TAB}/simeval"
    S = "20_20"
    TY = ["DEL", "INS", "DUP", "INV", "dDUP", "rTRA"]
    SZ = [("50–100 bp", 50, 100), ("0.1–1 kb", 100, 1000), ("1–10 kb", 1000, 10000), ("10–100 kb", 10000, 100000), (">100 kb", 100000, 10**12)]
    tk, hk, ts, hs = [collections.Counter() for _ in range(4)]
    for t in T:
        d = f"{P}/{S}/{t}"; L, K, X, J = load(f"{d}/calls.vcf"); mt = round(2 ** X[:, 1].min() - 1)
        pr = prob(features(X, f"{d}/t.cnt.hist", mt)); du = G["mark"](K, pr); cut = cutfor(mt)
        v = f"{V}/{S}/{t}"; os.makedirs(v, exist_ok=True)
        with open(f"{v}/verdict.vcf", "w") as o:
            o.write("".join(l for l in open(f"{d}/calls.vcf") if l[0] == "#"))
            for (g, info), q, u in zip(L, pr, du):
                g = list(g); g[6] = "DUPL" if u else "PASS" if q >= cut else "LOW"; g[7] += ";VP=%.4f" % q; o.write("\t".join(g) + "\n")
        sv, TJ = ns["truth"](f"{TR}/{t}.tsv")
        ln = {l.split("\t")[0]: int(l.split("\t")[2]) for l in open(f"{TR}/{t}.tsv")}
        f = ns["score"](f"{TR}/{t}.tsv", f"{v}/verdict.vcf")["found_ids"]
        for m, ty in sv.items():
            tk[ty] += 1; hk[ty] += m in f
            s = next(n for n, lo, hi in SZ if lo <= ln[m] < hi); ts[s] += 1; hs[s] += m in f
    fig, ax = page(160, 57.5)
    for x0, w, keys, T_, H_, lab, tg in ((14, 64, TY, tk, hk, "SV type", "a"), (94, 64, [s[0] for s in SZ], ts, hs, "SV size", "b")):
        a = ax(x0, 12.5, w, 38)
        r = [H_[q] / T_[q] for q in keys]
        a.bar(range(len(keys)), r, color=JANUS, width=0.65)
        for j, q in enumerate(keys): a.text(j, r[j] + 0.015, f"{r[j]:.3f}", ha="center", va="bottom", fontsize=FS - 1.0, bbox=dict(fc="white", ec="none", pad=0.2))
        a.set_xticks(range(len(keys))); a.set_xticklabels([f"{q}\n({T_[q]:,})" for q in keys]); a.set_xlim(-0.55, len(keys) - 0.45); a.set_ylim(0, 0.95); fticks(a, "y", [0, 0.25, 0.5, 0.75])
        a.set_ylabel("Recall"); a.set_xlabel(lab + " (simulated SVs)"); a.grid(axis="y", **GRID); a.set_axisbelow(True); tag(a, tg, x0 - 12, 2.5)
        print(lab, {q: round(H_[q] / T_[q], 3) for q in keys})
    fig.savefig(f"{OUT}/fig_simeval.pdf", facecolor="white", dpi=200)
    src("fig_simeval", ["panel", "class", "simulated_svs", "found", "recall"], [["a", q, tk[q], hk[q], hk[q] / tk[q]] for q in TY] + [["b", q, ts[q], hs[q], hs[q] / ts[q]] for q, _, _ in SZ])
    VAL["simeval"] = {q: hk[q] / tk[q] for q in TY}




def ms(x):
    x = round(x); return f"{x // 60}:{x % 60:02d}"


def hm(x):
    x = round(x / 60); return f"{x // 60}:{x % 60:02d}"


def hmin(x):
    x = round(x / 60); return f"{x // 60} h {x % 60} min"


def rng(v, d=1, sep="--"):
    v = list(v); a, b = f"{min(v):.{d}f}", f"{max(v):.{d}f}"
    return a if a == b else f"{a}{sep}{b}"


def irng(v, sep="--"):
    v = list(v); return f"{min(v):,}" if min(v) == max(v) else f"{min(v):,}{sep}{max(v):,}"


def tex(n):
    return "\\tool{}" if n == "janus" else n


def bold(vals, fmt="{:.3f}"):
    s = [fmt.format(v) for v in vals]; m = max(s, key=float)
    return [f"\\textbf{{{x}}}" if x == m else x for x in s]


def tabular(cols, head, rows):
    return "\\begin{tabular}{" + cols + "}\n\\toprule\n" + head + "\n\\midrule\n" + "\n".join(r if r.startswith(("\\midrule", "\\cmidrule")) or r.endswith("\\\\") else r + " \\\\" for r in rows) + "\n\\bottomrule\n\\end{tabular}\n"


def run(d, p):
    return json.load(open(f"{d}/full/{p}/run.json"))


def logs(p, s):
    return open(f"{RC}/full/{p}/logs/log.{s}").read()


def grab(p, s, pat):
    return int(re.search(pat, logs(p, s)).group(1))


def timing():
    c = cost(); tot = {p: sum(c["time"][p].values()) for p in P5}
    aln = lambda p: sum(json.load(open(f"{R}/align/{p}/{x}/run.json"))["total"]["wall_s"] for x in "nt")
    gbp = lambda p: sum(run(R, p)["reads"][x]["bases"] for x in "nt") / 1e9
    tp = {p: aln(p) / gbp(p) for p in ("hcc1954", "colo829")}; thr = st.mean(tp.values())
    mapt = {p: aln(p) if p == "hcc1954" else gbp(p) * thr for p in P5}
    phase = sum(float(l.split("\t")[1]) for l in open(f"{R}/phase/hcc1954/times.tsv") if not l.startswith("step"))
    pcpu = sum(float(l.split("\t")[2]) for l in open(f"{R}/phase/hcc1954/times.tsv") if not l.startswith("step"))
    acpu = sum(json.load(open(f"{R}/align/hcc1954/{x}/run.json"))["total"]["cpu_s"] for x in "nt")
    apk = max(json.load(open(f"{R}/align/hcc1954/{x}/run.json"))["total"]["peak_rss_kib"] for x in "nt") / 2**20
    scpu = sum(x["cpu_s"] for x in run(R, "hcc1954")["steps"] if x["step"].startswith("pack")) + run(RC, "hcc1954")["total"]["cpu_s"]
    pipe = {}
    for n, t, ph in [("Severus", "severus_phased", 1), ("SAVANA", "savana_phased", 1), ("nanomonsv", "nanomonsv", 0), ("dysgu", "dysgu", 0), ("SVision-pro", "svision-pro", 0),
                     ("cuteSV", "cutesv", 0), ("DeBreak", "debreak", 0), ("Sniffles2", "sniffles2", 0), ("kled", "kled", 0), ("SVDSS", "svdss", 0)]:
        j = json.load(open(f"{R}/competitors/{t}/hcc1954/run.json"))["total"]
        pipe[n] = dict(wall=mapt["hcc1954"] + phase * ph + j["wall_s"], cpu=acpu + pcpu * ph + j["cpu_s"])
    return dict(cost=c, tot=tot, mapt=mapt, tp=tp, gbp={p: gbp(p) for p in P5}, phase=phase, acpu=acpu, apk=apk, scpu=scpu, pipe=pipe)


def hmh(s):
    t = s / 3600; return f"{int(t)}:{int(round((t - int(t)) * 60)):02d}"


def supp_tsv(f):
    return list(csv.DictReader(open(f"{H}/supp/{f}"), delimiter="\t"))


SETTING = ["full", "30x", "20x", "15x", "10x", "5x"]; PURITY = ["100p", "80p", "60p", "40p", "20p"]


def tables():
    T, Z = {}, timing()
    rows = []
    for p in P5:
        r = D["all_with"][p]["janus"]; j = run(R, p); X = j["params"]["bases_per_x"]
        rows.append(f"{PAIRNAME[p]} & {j['reads']['n']['bases'] / X:.1f}/{j['reads']['t']['bases'] / X:.1f} & {r['truth']:,} & {r['calls']:,} & {r['tp']:,} & {r['fp']:,} & {r['fn']:,} & {r['P']:.3f} & {r['R']:.3f} & {r['F1']:.3f}")
    T["tab:castle"] = tabular("lrrrrrrrrr", "Pair & Coverage & Truth & Calls & TP & FP & FN & Precision & Recall & F1 \\\\\n     & n/t      &       &       &    &    &    &           &        &    \\\\", rows)
    ex = {}
    for l in open(f"{R}/external_eval/scores.tsv"):
        f = l.split("\t")
        if f[1] in EXT: ex[(f[0], EXT[f[1]])] = (int(f[2]), float(f[5]), float(f[6]), float(f[7]))
    names = sorted({n for _, n in ex}, key=lambda n: -ex[("colo829", n)][3])
    b = {q: bold([ex[(q, n)][3] for n in names]) for q in ("colo829", "hg008")}
    rows = [f"{tex(n)} & " + " & ".join(f"{ex[(q, n)][0]} & {ex[(q, n)][1]:.3f} & {ex[(q, n)][2]:.3f} & {b[q][i]}" for q in ("colo829", "hg008")) for i, n in enumerate(names)]
    xe = SIG["external"]; bn = {v: k for k, v in TOOLKEY.items()}
    rows += ["\\midrule", "\\tool{} $-$ best other, F1 & " + " & ".join(f"\\multicolumn{{4}}{{r}}{{{xe[q]['diff']:+.3f} [{xe[q]['lo']:+.3f}, {xe[q]['hi']:+.3f}] vs.\\ {bn[xe[q]['best']]}}}" for q in ("colo829", "hg008"))]
    T["tab:external"] = tabular("lrrrrrrrr", " & \\multicolumn{4}{c}{COLO829} & \\multicolumn{4}{c}{HG008} \\\\\nCall set & Calls & Precision & Recall & F1 & Calls & Precision & Recall & F1 \\\\", rows)
    hd = " & ".join(PAIRNAME[p] for p in P5)
    pk = lambda t: json.load(open(f"{R}/competitors/{t}/hcc1954/run.json"))["total"]["peak_rss_kib"] / 2**20
    TK = {"Severus": "severus_phased", "SAVANA": "savana_phased", "nanomonsv": "nanomonsv", "dysgu": "dysgu", "SVision-pro": "svision-pro", "cuteSV": "cutesv", "DeBreak": "debreak", "Sniffles2": "sniffles2", "kled": "kled", "SVDSS": "svdss"}
    jt = Z["tot"]["hcc1954"]; jc = Z["scpu"]; jm = max(Z["cost"]["mem"]["hcc1954"].values())
    rows = [f"\\tool{{}} & -- & -- & -- & {hm(jt)} & 1 & {jc / 3600:.1f} & {jm:.1f}"]
    for n in sorted(OURS, key=lambda n: Z["pipe"][n]["wall"]):
        w = Z["pipe"][n]["wall"]; ph = Z["phase"] if n in ("Severus", "SAVANA") else 0; cl = w - Z["mapt"]["hcc1954"] - ph
        rows.append(f"{tex(n)} & {hm(Z['mapt']['hcc1954'])} & {hm(ph) if ph else '--'} & {hm(cl)} & {hm(w)} & {w / jt:.0f} & {Z['pipe'][n]['cpu'] / 3600:.0f} & {max(pk(TK[n]), Z['apk']):.1f}")
    T["stab:pipelines"] = tabular("lrrrrrrr", "Call set & Mapping & Phasing & Calling & Total & Factor & CPU (h) & Peak RSS (GB) \\\\", rows)
    c = Z["cost"]
    rows = [f"{g} & " + " & ".join(ms(c["time"][p][g]) for p in P5) + f" & {rng(c['mem'][p][g] for p in P5)}" for g, _ in GROUP]
    rows += ["\\midrule", "Total & " + " & ".join(ms(Z["tot"][p]) for p in P5) + f" & {rng(max(c['mem'][p].values()) for p in P5)}"]
    T["tab:time"] = tabular("lrrrrrr", f"Stage & {hd} & Peak RSS (GB) \\\\", rows)
    rows = []
    for p in P5:
        e = "" if p == "hcc1954" else "$^\\dagger$"
        sn_ = int(re.search(r"pieces, (\d+) bp", open(f"{RC}/full/{p}/logs/log.select").read()).group(1)) / 1e9
        rows.append(f"{PAIRNAME[p]} & {Z['gbp'][p]:.0f} & {sn_:.2f} & {Z['gbp'][p] / sn_:,.0f} & {hm(Z['mapt'][p])}{e} & {Z['mapt'][p] / Z['tot'][p]:.1f}{e}")
    T["tab:map"] = tabular("lrrrrr", "Pair & Reads (Gbp) & Snippets (Gbp) & Ratio & Mapping (h:mm) & Ratio \\\\", rows)
    S = simdata(); TR, MP = simtrees(); rows = []
    f1 = lambda d: [2 * a * b / (a + b) for a, b in zip(d["P"], d["R"])]
    for s in SIMSETS:
        t = S[s]["time"]; q = [st.median(v) for k_ in ("pair", "first", "joint") for v in (S[s][k_]["P"], S[s][k_]["R"], f1(S[s][k_]))]
        rows.append(f"{s.replace('-', '/')} & " + " & ".join(f"{x:.3f}" for x in q) + f" & {MP[s]['shared']}/{MP[s]['true']} & {MP[s]['rf']} & {sum(t['pair_wall']) / 60:.1f} & {t['joint_wall'] / 60:.1f} & {max(t['pair_peak_kb']) / 2**20:.1f} & {t['joint_peak_kb'] / 2**20:.1f}")
    T["tab:sim"] = tabular("l" + "r" * 15, " & \\multicolumn{3}{c}{Pairwise runs} & \\multicolumn{3}{c}{Joint run, first pass} & \\multicolumn{3}{c}{Joint run, clade pass} & \\multicolumn{2}{c}{Tree} & \\multicolumn{2}{c}{Wall clock (min)} & \\multicolumn{2}{c}{Peak RSS (GB)} \\\\\n"
                           "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-12}\\cmidrule(lr){13-14}\\cmidrule(lr){15-16}\nDepth n/t & Precision & Recall & F1 & Precision & Recall & F1 & Precision & Recall & F1 & Splits & RF & pairwise & joint & pairwise & joint \\\\", rows)
    rows = []
    for s in SIMSETS:
        a_, b_ = S[s]["first"], S[s]["joint"]; n1, t1, n2, t2 = sum(a_["carried"]), sum(a_["carried_true"]), sum(b_["carried"]), sum(b_["carried_true"])
        rows.append(f"{s.replace('-', '/')} & {n1:,} & {t1:,} & {t1 / max(n1, 1):.3f} & {n2:,} & {t2:,} & {t2 / max(n2, 1):.3f} & " + " & ".join(f"{x[2]:.3f}" for x in S[s]["mu"]["PRF"]))
    PO = json.load(open(f"{R}/sim/pooled/scores.json"))["20-5"]; CL = list(PO["joint"]["classes"])
    rows = [f"{lab} & {o['calls']:,} & {o['tp']:,} & {o['found']:,} & {o['P']:.3f} & {o['R']:.3f} & {o['F1']:.3f} & " + " & ".join(f"{o['classes'][c][1]:.3f}" for c in CL)
            for key, lab in (("pooled", "20 tumors pooled into one sample"), ("pairwise", "20 pairwise runs"), ("joint", "Joint run")) for o in [PO[key]]]
    T["stab:pooled"] = tabular("lrrrrrrrrrr", " & & & & & & & \\multicolumn{4}{c}{Recall by the number of tumors carrying the SV} \\\\\n\\cmidrule(lr){8-11}\nRun & Calls & True & SVs found & Precision & Recall & F1 & " + " & ".join(f"{c} ({PO['joint']['classes'][c][0]:,})" for c in CL) + " \\\\", rows)
    T["stab:clade"] = tabular("l" + "r" * 11, " & \\multicolumn{3}{c}{Carriers, first pass} & \\multicolumn{3}{c}{Carriers, clade pass} & \\multicolumn{5}{c}{Median F1 of the clade pass by $\\mu$} \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-12}\nDepth n/t & Reported & True & Precision & Reported & True & Precision & 0.3 & 0.5 & 1.0 & 1.5 & 2.0 \\\\", rows)
    T["stab:simtumor"] = tabular("l" + "r" * 15, " & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{{s_.replace('-', '/')}}}" for s_ in SIMSETS) + " \\\\\n" + "".join(f"\\cmidrule(lr){{{2 + 3 * k}-{4 + 3 * k}}}" for k in range(5)) + "\nTumor" + " & pairwise & first pass & clade pass" * 5 + " \\\\",
                                 [f"T{i + 1:02d} & " + " & ".join(f"{f1(S[s_][m])[i]:.3f}" for s_ in SIMSETS for m in ("pair", "first", "joint")) for i in range(20)])
    TECHN = {"HiFi": "PacBio HiFi", "ONT": "ONT", "Illumina": "Illumina"}
    al = sorted(D["all_with"]["hcc1954"], key=lambda n: -mean("all_with", n, "F1"))
    rows = [f"{tex(base(n))} & {TECHN[tech(n)]} & " + " & ".join(f"{D['all_with'][p][n]['F1']:.3f}" for p in P5) + " & " + " & ".join(f"{mean('all_with', n, k):.3f}" for k in ("P", "R", "F1")) for n in al]
    T["stab:allwith"] = tabular("llrrrrrrrr", f"Call set & Technology & {hd} & Precision & Recall & F1 \\\\", rows)
    rows = []
    for kind, keys, labs in (("type", TYPES, ["Breakend", "Deletion", "Duplication", "Insertion", "Inversion"]), ("size", [k for k, _ in SIZES], ["no length", "50--100\\,bp", "0.1--1\\,kb", "1--10\\,kb", "10--100\\,kb", "$>$100\\,kb"])):
        if rows: rows.append("\\midrule")
        for i, (k, lab) in enumerate(zip(keys, labs)):
            v = [pooled("all_with", n, kind, k) for n in CAL5]
            rows.append(f"{(kind.capitalize() if i == 0 else '')} & {lab} & {v[0][1]:,} & " + " & ".join(bold([x[0] for x in v])))
    T["stab:types"] = tabular("llrrrrrr", " & Class & Consensus SVs & \\tool{} & Severus & SAVANA & nanomonsv & Sniffles2 \\\\", rows)
    L = loo(); lm = lambda n, k: st.mean(L[p][n][k] for p in P5)
    rows = [f"{tex(n)} & " + " & ".join(f"{mean('all_without', n, k):.3f}" for k in ("P", "R", "F1")) + " & " + " & ".join(f"{lm(n, k):.3f}" for k in ("P", "R", "F1")) for n in sorted(L["hcc1954"], key=lambda n: -lm(n, "F1"))]
    T["stab:loo"] = tabular("lrrrrrr", " & \\multicolumn{3}{c}{Consensus of the 17 call sets other than \\tool{}} & \\multicolumn{3}{c}{Consensus of the other 17 call sets} \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\nCall set & Precision & Recall & F1 & Precision & Recall & F1 \\\\", rows)
    C_ = cov(); sg = SIG["settings"]; rows = []
    for k, lab in [("full", "Full")] + [(f"{x}x", f"{x}$\\times$") for x in (30, 20, 15, 10, 5)] + [("mid", None)] + [(f"{q}p", f"{q}\\%") for q in (100, 80, 60, 40, 20)]:
        if k == "mid": rows += ["\\midrule", "Purity (30$\\times$) & & & & & & & \\\\", "\\midrule"]; continue
        g = sg[k]; bn = {v: k_ for k_, v in TOOLKEY.items()}[g["best"]]
        rows.append(f"{lab} & " + " & ".join(f"{series(C_, 'janus', [k], P5, j)[0]:.3f}" for j in range(3)) + f" & {g['snape_lo']:.3f}--{g['snape_hi']:.3f} & " + " & ".join(f"{series(C_, n, [k], P5, 2)[0]:.3f}" for n in COV) + f" & {g['diff']:+.3f} [{g['lo']:+.3f}, {g['hi']:+.3f}] vs.\\ {bn}")
    T["stab:depth"] = tabular("lrrrrrrrl", "Depth & \\multicolumn{4}{c}{\\tool{}} & Severus & SAVANA & nanomonsv & \\tool{} $-$ best other \\\\\n\\cmidrule(lr){2-5}\n & Precision & Recall & F1 & 95\\% interval & F1 & F1 & F1 & F1 [95\\% interval] \\\\", rows)
    tot, bam, stat = VAL["files"]
    T["stab:files"] = tabular("lrrrrrrr", "Pair & Packed reads & Core positions & Counts & Graph & Snippets, alignments, other & Total & BAM files \\\\",
                              [f"{PAIRNAME[p]} & " + " & ".join(f"{x:.1f}" for x in tot[p]) + f" & {sum(bam[p]):.0f}" for p in P5])
    ST = [("Nodes removed as tips or bubbles", lambda p: grab(p, "prune", r"(\d+) dropped nodes")), ("Sites", lambda p: grab(p, "prune", r"dropped nodes, (\d+) sites")),
          ("Additional rows from distance groups", lambda p: grab(p, "select", r"(\d+) more rows")), ("Reads from which snippets were cut", lambda p: grab(p, "select", r"on (\d+) reads")),
          ("Snippets", lambda p: grab(p, "select", r"reads; (\d+) pieces")), ("Sequence in snippets (Mbp)", lambda p: round(grab(p, "select", r"pieces, (\d+) bp") / 1e6)),
          ("Snippets placed on the reference", lambda p: grab(p, "place", r"(\d+) placed")), ("Sites with an event", lambda p: grab(p, "call", r"\[call\] (\d+) of \d+ sites carry")),
          ("Records", lambda p: grab(p, "call", r"at sites, (\d+) records")), ("PASS calls after classification", lambda p: sum(1 for l in open(f"{H}/verdict/{p}.vcf") if "\tPASS\t" in l))]
    T["stab:steps"] = tabular("lrrrrr", f" & {hd} \\\\", [f"{n} & " + " & ".join(f"{fn(p):,}" for p in P5) for n, fn in ST])
    A = [r for r in supp_tsv("heldout_auc.tsv")]
    am = lambda k, x: st.mean(float(r[x]) for r in A if r["setting"] == k)
    rows = [f"{lab} & " + " & ".join(f"{am(k, x):.3f}" for k in SETTING) for lab, x in (("ROC-AUC", "roc_auc"), ("Average precision", "pr_auc"))]
    rows += ["\\midrule", "Purity (30$\\times$) & 100\\% & 80\\% & 60\\% & 40\\% & 20\\% & ", "\\midrule"]
    rows += [f"{lab} & " + " & ".join(f"{am(k, x):.3f}" for k in PURITY) + " & " for lab, x in (("ROC-AUC", "roc_auc"), ("Average precision", "pr_auc"))]
    sh = "Setting & full & 30$\\times$ & 20$\\times$ & 15$\\times$ & 10$\\times$ & 5$\\times$ \\\\"
    T["stab:auc"] = tabular("lrrrrrr", sh, rows)
    SP = supp_tsv("splits_3_2.tsv"); sp = collections.defaultdict(list)
    for r in SP: sp[(r["train"], r["setting"])].append(float(r["F1"]))
    spm = lambda k: [st.mean(v) for (tr, s), v in sp.items() if s == k]
    def three(keys, pad):
        return [f"Four training pairs & " + " & ".join(f"{am(k, 'F1'):.3f}" for k in keys) + pad, f"Three training pairs, mean & " + " & ".join(f"{st.mean(spm(k)):.3f}" for k in keys) + pad,
                f"Three training pairs, range & " + " & ".join(rng(spm(k), 3) for k in keys) + pad]
    T["stab:splits"] = tabular("lrrrrrr", sh, three(SETTING, "") + ["\\midrule", "Purity (30$\\times$) & 100\\% & 80\\% & 60\\% & 40\\% & 20\\% & ", "\\midrule"] + three(PURITY, " & "))
    AB = supp_tsv("ablation.tsv"); VN = {"max_iter=100": "100 iterations", "max_iter=600": "600 iterations", "learning_rate=0.02": "Learning rate 0.02", "learning_rate=0.1": "Learning rate 0.1",
                                         "max_leaf_nodes=7": "At most 7 leaves", "max_leaf_nodes=31": "At most 31 leaves", "min_samples_leaf=10": "At least 10 records per leaf",
                                         "min_samples_leaf=160": "At least 160 records per leaf", "l2_regularization=0.0": "L2 regularization 0", "l2_regularization=5.0": "L2 regularization 5", "seeds=1": "One seed instead of five"}
    abm = lambda v, k: st.mean(float(r["F1"]) for r in AB if r["variant"] == v and r["setting"] == k)
    AK = ["full", "20x", "5x", "20p"]
    rows = ["Parameters used in this paper & " + " & ".join(f"{am(k, 'F1'):.3f}" for k in AK)] + [f"{lab} & " + " & ".join(f"{abm(v, k):.3f}" for k in AK) for v, lab in VN.items()]
    T["stab:ablation"] = tabular("lrrrr", "Variant & full & 20$\\times$ & 5$\\times$ & 20\\% purity \\\\", rows)
    GR = {j: collections.defaultdict(dict) for j in ("grid",)}; GS = [s for s in SETTING + PURITY if s != "100p"]; GK = AK + ["all"]
    for j in GR:
        for r in supp_tsv(f"{j}.tsv"): GR[j][r["variant"]][(r["pair"], r["setting"])] = (float(r["F1"]), float(r["inner_f1"]))
    gm = lambda d, k: st.mean(d[(p, s)][0] for p in P5 for s in GS if k in ("all", s))
    gpick = lambda j: {p: max(GR[j], key=lambda v: GR[j][v][(p, "full")][1]) for p in P5}
    gnest = lambda j, k: st.mean(GR[j][g[p]][(p, s)][0] for g in [gpick(j)] for p in P5 for s in GS if k in ("all", s))
    gbest = lambda j: max(GR[j], key=lambda v: gm(GR[j][v], "all"))
    used = {k: st.mean(am(s, "F1") for s in GS) if k == "all" else am(k, "F1") for k in GK}
    Z["grid"] = {"used": used, **{j: {"n": len(GR[j]), "best": {k: gm(GR[j][gbest(j)], k) for k in GK}, "nested": {k: gnest(j, k) for k in GK}} for j in GR}}
    rows = ["Parameters used in this paper & " + " & ".join(f"{used[k]:.3f}" for k in GK)]
    for j, lab in (("grid", "scikit-learn"),):
        rows += ["\\midrule", f"\\multicolumn{{6}}{{l}}{{{lab}, {len(GR[j])} sets}} \\\\",
                 "\\quad Best set (chosen on the held-out pairs) & " + " & ".join(f"{Z['grid'][j]['best'][k]:.3f}" for k in GK) + " \\\\",
                 "\\quad Sets chosen by nested cross-validation & " + " & ".join(f"{Z['grid'][j]['nested'][k]:.3f}" for k in GK) + " \\\\",
                 f"\\quad Range of the {len(GR[j])} sets & " + " & ".join(rng((gm(d, k) for d in GR[j].values()), 3) for k in GK) + " \\\\"]
    T["stab:grid"] = tabular("lrrrrr", "Parameter set & Full & 20$\\times$ & 5$\\times$ & 20\\% purity & All \\\\", rows)
    TD = os.environ.get("TABTEX", f"{os.path.dirname(OUT)}/tables"); os.makedirs(TD, exist_ok=True)
    for k, v in T.items():
        with open(f"{TD}/{k.split(':')[1]}.tex", "w") as o: o.write(v)
    with open(f"{TAB}/tables.tex", "w") as o:
        for k, v in T.items(): o.write(f"% {k}\n{v}\n")
    return T, Z


def numbers(Z):
    N = []
    put = lambda k, v: N.append((k, v))
    j = {p: run(R, p) for p in P5}; X = 3.1e9
    put("depth of the normal samples", rng((j[p]["reads"]["n"]["bases"] / X for p in P5), 0, "--") + "$\\times$")
    put("depth of the tumor samples", rng((j[p]["reads"]["t"]["bases"] / X for p in P5), 0, "--") + "$\\times$")
    f = {p: D["all_with"][p]["janus"] for p in P5}
    lo, hi = min(P5, key=lambda p: f[p]["F1"]), max(P5, key=lambda p: f[p]["F1"])
    put("F1 range of janus", f"from {f[lo]['F1']:.3f} for {PAIRNAME[lo]} to {f[hi]['F1']:.3f} for {PAIRNAME[hi]}")
    put("lowest precision / recall", f"{PAIRNAME[min(P5, key=lambda p: f[p]['P'])]} / {PAIRNAME[min(P5, key=lambda p: f[p]['R'])]} (consensus {f[min(P5, key=lambda p: f[p]['R'])]['truth']} SVs)")
    lab = {}
    for l in open(f"{H}/supp/labels11/scores.tsv"):
        x = l.rstrip("\n").split("\t"); lab[(x[0], x[1])] = x[2:]
    l17 = lambda k: st.mean(float(lab[(p, k)][2]) for p in P5)
    put("mean F1 with the 11-call-set labels, against 18 / against 17 call sets", f"{l17('with'):.3f} instead of {mean('all_with', 'janus', 'F1'):.3f} / {l17('without'):.3f} instead of {mean('all_without', 'janus', 'F1'):.3f}")
    put("labels11 HG008 and COLO829 F1, COLO829 calls", f"HG008 {float(lab[('hg008', 'truth')][2]):.3f}, COLO829 {float(lab[('colo829', 'truth')][2]):.3f}, {lab[('colo829', 'truth')][3]} calls")
    C_ = cov(); ks = ["30x", "20x", "15x", "10x", "5x"]
    put("labels11 mean F1 at 30x..5x", ", ".join(f"{l17(k):.3f}" for k in ks) + " compared with " + ", ".join(f"{x:.3f}" for x in series(C_, "janus", ks, P5, 2)))
    put("labels11 largest depth difference", f"{max(abs(round(l17(k), 3) - round(x, 3)) for k, x in zip(ks, series(C_, 'janus', ks, P5, 2))):.3f}")
    allc = list(D["all_with"]["hcc1954"])
    for n in ["janus", "Severus", "Severus ONT", "SAVANA", "nanomonsv", "Manta", "GRIPSS", "Sniffles2 ONT", "SvABA"]:
        put(f"mean F1 against 18 call sets, {n}", f"{mean('all_with', n, 'F1'):.3f}")
    put("pairs where janus ranks first (18 call sets)", ", ".join(PAIRNAME[p] for p in P5 if max(allc, key=lambda n: D["all_with"][p][n]["F1"]) == "janus"))
    put("highest precision of all call sets", max(allc, key=lambda n: mean("all_with", n, "P")))
    put("highest recall of the HiFi call sets", max(["janus"] + OURS, key=lambda n: mean("all_with", n, "R")))
    wo = sorted(OURS + ["janus"], key=lambda n: -mean("all_without", n, "F1"))[:3]
    put("top three against 17 call sets", ", ".join(f"{n} {mean('all_without', n, 'F1'):.3f}" for n in wo))
    L = loo(); lm = lambda n: st.mean(L[p][n]["F1"] for p in P5)
    put("leave-one-out top three", ", ".join(f"{n} {lm(n):.3f}" for n in sorted(L["hcc1954"], key=lambda n: -lm(n))[:3]))
    put("pairs where janus ranks first (leave one out)", ", ".join(PAIRNAME[p] for p in P5 if max(L[p], key=lambda n: L[p][n]["F1"]) == "janus"))
    for q in ("colo829", "hg008"):
        e = SIG["external"][q]
        put(f"{q}: difference to {e['best']} and 95% interval", f"{e['diff']:.3f} ({e['lo']:.3f} to {e['hi']:.3f})")
    put("consensus SVs, five pairs", f"{sum(D['all_with'][p]['janus']['truth'] for p in P5):,}")
    put("consensus without janus, sizes", ", ".join(f"{D['all_without'][p]['Severus']['truth']:,}" for p in P5))
    s30, s5 = series(C_, "janus", ["30x", "5x"], P5, 2), series(C_, "janus", ["30x", "5x"], P5, 1)
    v30, v5 = series(C_, "Severus", ["30x", "5x"], P5, 2), series(C_, "Severus", ["30x", "5x"], P5, 1)
    put("janus mean F1 and recall 30x to 5x", f"F1 {s30[0]:.3f} to {s30[1]:.3f}, recall {s5[0]:.3f} to {s5[1]:.3f}")
    put("Severus mean F1 30x to 5x, recall at 5x", f"{v30[0]:.3f} to {v30[1]:.3f}, recall {v5[1]:.3f}")
    put("20% purity, mean over five pairs: janus precision and F1, Severus F1", ", ".join(f"{st.mean(C_[(p, '20p', t)][j] for p in P5):.3f}" for t, j in (("snape", 0), ("snape", 2), ("severus", 2))))
    put("20% purity, Severus F1 on NCI-H2009, HCC1937, NCI-H1437", rng((C_[(p, "20p", "severus")][2] for p in ("h2009", "hcc1937", "h1437")), 3))
    put("purity, mean F1 100% to 20%, janus and Severus", "; ".join(", ".join(f"{st.mean(C_[(p, q, t)][2] for p in P5):.3f}" for q in ("100p", "80p", "60p", "40p", "20p")) for t in ("snape", "severus")))
    put("20% purity F1 of the other three pairs", rng((C_[(p, "20p", "snape")][2] for p in ("h2009", "hcc1937", "h1437")), 3))
    c, tot = Z["cost"], Z["tot"]
    def mins(x):
        s = round(x); return f"{s // 60} min {s % 60} s" if s % 60 else f"{s // 60} min"
    put("run time per pair", f"{mins(min(tot.values()))} to {mins(max(tot.values()))} (mean {st.mean(tot.values()) / 60:.1f} min)")
    put("packing and core extraction", rng(((c["time"][p]["Packing"] + c["time"][p]["Core extraction"]) / 60 for p in P5), 1) + " min")
    put("graph and site detection, longest", f"{max((c['time'][p]['Counting'] + c['time'][p]['Graph and pruning']) / 60 for p in P5):.1f} min")
    put("selection, placement and alignment", rng(((c["time"][p]["Selection"] + c["time"][p]["Placement"] + c["time"][p]["minimap2"]) / 60 for p in P5), 0) + " min")
    put("peak memory (minimap2)", rng(max(c["mem"][p].values()) for p in P5) + " GB")
    put("peak memory of the janus steps", f"{max(v for p in P5 for g, v in c['mem'][p].items() if g != 'minimap2'):.1f} GB")
    F = json.load(open(f"{H}/verdict/folds.json"))
    put("training one model", rng((F[p]["train_s"] for p in F), 0) + " s")
    put("classification step", f"{float(re.search(r'([\d.]+) s$', open(f'{H}/verdict/log.txt').read().strip().splitlines()[-1]).group(1)) / 60:.1f} minutes")
    tot_f, bam, _ = VAL["files"]
    put("intermediate files per pair", rng((tot_f[p][5] for p in P5), 0) + " GB")
    put("two BAM files per pair", rng((sum(bam[p]) for p in P5), 0) + f" GB, on average {100 * (st.mean(sum(bam[p]) / tot_f[p][5] for p in P5) - 1):.0f}% more")
    put("packed reads / core positions / rest", " / ".join(rng((tot_f[p][i] for p in P5), 0) for i in (0, 1)) + " / " + rng((sum(tot_f[p][2:5]) for p in P5), 0) + " GB")
    put("intermediate files as share of the BAM files", rng((100 * tot_f[p][5] / sum(bam[p]) for p in P5), 0) + "%")
    put("mapping HCC1954", f"{hmin(Z['mapt']['hcc1954'])}, {Z['apk']:.1f} GB")
    put("janus share of mapping time and memory (HCC1954)", f"{100 * tot['hcc1954'] / Z['mapt']['hcc1954']:.0f}% and {100 * max(c['mem']['hcc1954'].values()) / Z['apk']:.0f}% ({max(c['mem']['hcc1954'].values()):.1f} GB)")
    put("estimated mapping time, other four pairs", rng((Z["mapt"][p] / 3600 for p in P5 if p != "hcc1954"), 1) + " hours")
    put("mapping throughput HCC1954, COLO829", " and ".join(f"{Z['tp'][p]:.1f}" for p in ("hcc1954", "colo829")) + " s per Gbp")
    pw = {n: v["wall"] for n, v in Z["pipe"].items()}
    put("pipelines on HCC1954", f"from {hmin(min(pw.values()))} ({', '.join(n for n in pw if hmin(pw[n]) == hmin(min(pw.values())))}) to {hmin(max(pw.values()))} ({max(pw, key=pw.get)}); Severus {hmin(pw['Severus'])}, SAVANA {hmin(pw['SAVANA'])}, phasing {Z['phase'] / 60:.0f} min")
    gb = {p: grab(p, "select", r"pieces, (\d+) bp") / 1e9 for p in P5}
    put("aligned pieces", rng(gb.values(), 2) + " Gbp, " + rng((100 * gb[p] / Z["gbp"][p] for p in P5), 2) + "% of " + rng(Z["gbp"].values(), 0) + " Gbp")
    put("CPU time on HCC1954: janus, mapping, pipelines", f"{Z['scpu'] / 3600:.1f} h, {Z['acpu'] / 3600:.0f} h, " + rng((v["cpu"] / 3600 for v in Z["pipe"].values()), 0) + " h")
    put("nodes removed as tips or bubbles", rng((grab(p, "prune", r"(\d+) dropped nodes") / 1e6 for p in P5), 1) + " million")
    put("tumor reads, sites, pieces, records", "; ".join([rng((j[p]["reads"]["t"]["reads"] / 1e6 for p in P5), 1) + " million"] + [irng(round(grab(p, s, pat), -3) for p in P5) for s, pat in (("prune", r"dropped nodes, (\d+) sites"), ("select", r"reads; (\d+) pieces"))] + [irng(round(grab(p, "call", r"at sites, (\d+) records"), -2) for p in P5)]))
    sg = VAL["simgen"]
    put("simulated SVs, per tumor", f"{sg['n']:,}; " + irng(sg["per"]))
    put("simulated SVs by type", ", ".join(f"{v:,} {k}" for k, v in sorted(sg["types"].items(), key=lambda x: -x[1])))
    S = simdata(); TR, MP = simtrees()
    f1 = lambda d: np.array([2 * a * b / (a + b) for a, b in zip(d["P"], d["R"])])
    for s in SIMSETS:
        a, b = f1(S[s]["pair"]), f1(S[s]["joint"])
        put(f"{s}: joint run better for tumors, median F1 pairwise/joint", f"{int((b > a).sum())} of 20, {np.median(a):.3f}/{np.median(b):.3f}")
    t5 = S["20-5"]["time"]
    put("20/5 joint time, speedup, pairwise sum", f"{t5['joint_wall'] / 60:.1f} min, {sum(t5['pair_wall']) / t5['joint_wall']:.1f} times, {sum(t5['pair_wall']) / 60:.1f} min")
    put("joint memory above the largest pairwise run at", f"{sum(S[s]['time']['joint_peak_kb'] > max(S[s]['time']['pair_peak_kb']) for s in SIMSETS)} of 5 depths")
    M2 = mm2()
    for s, v in M2.items():
        put(f"minimap2 on the 21 read sets at {s}", f"{v['wall_min']:.1f} min ({hmin(v['wall_min'] * 60)}), {v['peak_gb']:.1f} GB")
    d = TR["sets"]["20-5"]
    put("20/5 joint SVs, on one clade", f"{d['events']:,}, {d['exact']:,} ({100 * d['exact'] / d['events']:.0f}%)")
    for s in SIMSETS:
        put(f"{s}: neighbor-joining splits, RF; compatible-set tree splits, RF", f"{MP[s]['shared']}/{MP[s]['true']}, {MP[s]['rf']}; {TR['sets'][s]['shared']}/{MP[s]['true']}, {TR['sets'][s]['rf']}")
    se = VAL["simeval"]
    put("20/20 pairwise recall by type", ", ".join(f"{k} {v:.2f}" for k, v in se.items()))
    put("joint rule sweep, F1 range at 20/20 and 20/5", "; ".join(rng(VAL["joint"][s], 3) for s in ("20-20", "20-5")))
    m = VAL["missed"]
    put("consensus SVs without a PASS call, with a record, without", f"{sum(m[p][1] for p in P5)}, {sum(m[p][2] for p in P5)}, {sum(m[p][3] for p in P5)}")
    put("Minda false negatives, five pairs", f"{sum(D['all_with'][p]['janus']['fn'] for p in P5)}")
    ev = VAL["evidence"]
    put("records matching / not matching the labels; share with a normal read |D|<20", f"{ev[1][0]:,} / {ev[0][0]:,}; {100 * ev[1][1]:.0f}% / {100 * ev[0][1]:.0f}%")
    au = [r for r in supp_tsv("heldout_auc.tsv") if r["setting"] == "full"]
    put("full-depth ROC-AUC range, average precision range", rng((float(r["roc_auc"]) for r in au), 3) + "; " + rng((float(r["pr_auc"]) for r in au), 3))
    SP = supp_tsv("splits_3_2.tsv")
    put("three-pair F1 range per pair at full depth", "; ".join(f"{PAIRNAME[p]} " + rng((float(r["F1"]) for r in SP if r["pair"] == p and r["setting"] == "full"), 3) for p in P5))
    ex = ACC["extra"]
    put("Sniffles2 mosaic mean F1, Delly F1 on HCC1954", f"{st.mean(ex['all_with'][p]['Sniffles2 (mosaic)']['F1'] for p in P5):.3f}, {ex['all_with']['hcc1954']['Delly']['F1']:.3f}")
    hifi = ["janus"] + OURS; mw = lambda n, k="F1": mean("all_with", n, k)
    put("call sets above janus in HCC1395", ", ".join(n for n in allc if D["all_with"]["hcc1395"][n]["F1"] > D["all_with"]["hcc1395"]["janus"]["F1"]))
    put("ONT below HiFi in F1 and precision", ", ".join(f"{n}: {mw(n + ' ONT') < mw(n)} {mw(n + ' ONT', 'P') < mw(n, 'P')}" for n in ("SAVANA", "nanomonsv")))
    put("HiFi call sets below Manta, GRIPSS, SvABA", "; ".join(f"{n}: {', '.join(x for x in hifi if mw(x) < mw(n))}" for n in ("Manta", "GRIPSS", "SvABA")))
    put("Illumina call sets, best and worst pair", "; ".join(f"{n}: {PAIRNAME[max(P5, key=lambda p: D['all_with'][p][n]['F1'])]}, {PAIRNAME[min(P5, key=lambda p: D['all_with'][p][n]['F1'])]}" for n in ILL))
    put("HiFi call sets with higher recall than nanomonsv", ", ".join(n for n in hifi if mw(n, "R") > mw("nanomonsv", "R")))
    for kind, keys in (("type", TYPES), ("size", [k for k, _ in SIZES])):
        put(f"callers above janus by {kind}", "; ".join(f"{k}: " + ", ".join(n for n in ["Severus", "SAVANA", "nanomonsv", "Sniffles2"] if pooled("all_with", n, kind, k)[0] > pooled("all_with", "janus", kind, k)[0]) for k in keys))
    put("aligned sequence, times less than both read sets", irng(round(Z["gbp"][p] / gb[p]) for p in P5))
    xt = {q: int(overall(f"{R}/external_eval/{q}/snape/minda")[0][8]) for q in ("colo829", "hg008")}
    ns = {}; exec(paper.SCRIPTS["score"].split("T, C = recs")[0], ns)
    put("COLO829 validated SVs, at least 50 bp; HG008 truth SVs of at least 50 bp", f"{len(ns['recs'](f'{ROOT}/hifi/truth/colo829/truth.vcf'))}, {xt['colo829']}; {xt['hg008']}")
    lc = {p: [int(overall(f"{H}/eval_loo/{p}/{t}/score")[0][8]) for t in EXT] for p in ("hcc1954", "h1437")}
    put("leave-one-out consensus sizes, HCC1954 and NCI-H1437", "; ".join(irng(v) for v in lc.values()))
    pn = {p: round(min(j[p]["reads"]["n"]["bases"] - 24 * X - X / 10, 30 * X) / X) for p in ("hcc1954", "hcc1395")}
    put("normal depth of the purity runs", ", ".join(f"{PAIRNAME[p]} {v}" for p, v in pn.items()))
    nr = len(glob.glob(f"{H}/verdict/*.vcf")) + len(glob.glob(f"{H}/verdict/coverage_matched/*.vcf")) + len(glob.glob(f"{H}/verdict/purity_matched/*.vcf"))
    put("classified real runs", str(nr))
    ff = []
    for p in P5:
        for l in open(f"{H}/verdict/{p}.vcf"):
            if l[0] != "#": ff.append(re.search(r"\bF=([^;\t]+)", l).group(1).split(","))
    ff = np.array(ff, float)
    put("feature values per record, constant ones", f"{ff.shape[1]}, {int((ff.max(0) == ff.min(0)).sum())}")
    mm = collections.defaultdict(set)
    for f in glob.glob(f"{R}/full/*/t.cnt") + glob.glob(f"{R}/coverage/*/t.cnt") + glob.glob(f"{R}/purity/*/t.cnt"):
        d = os.path.dirname(f); k = os.path.basename(d).rsplit("_", 1)[1] if "_" in os.path.basename(d) else "full"
        mm[(k, "m0")].add(open(f"{d}/n.cnt", "rb").read(27)[26]); mm[(k, "mt")].add(open(f, "rb").read(27)[26])
    put("minimum read support m0 by setting", ", ".join(f"{k} {sorted(v)}" for (k, x), v in sorted(mm.items()) if x == "m0"))
    put("minimum read support mt by setting", ", ".join(f"{k} {sorted(v)}" for (k, x), v in sorted(mm.items()) if x == "mt"))
    put("full-depth pairs with mt 2", ", ".join(os.path.basename(os.path.dirname(f)) for f in glob.glob(f"{R}/full/*/t.cnt") if open(f, "rb").read(27)[26] == 2))
    put("simulated runs, minimum read support", str(sorted({m for s in S.values() for m in s["pair"]["minimum"]} | {int(open(f).read()) for f in glob.glob(f"{R}/sim/settings/*/joint/*.minimum")})))
    rr = collections.defaultdict(list)
    for r in SP: rr[r["pair"]].append(float(r["F1"])) if r["setting"] == "full" else None
    A4 = {r["pair"]: float(r["F1"]) for r in supp_tsv("heldout_auc.tsv") if r["setting"] == "full"}
    put("three-pair range at full depth, widest; four-pair F1 inside every range", f"{max(max(v) - min(v) for v in rr.values()):.3f}; {all(min(rr[p]) <= A4[p] <= max(rr[p]) for p in P5)}")
    sp = collections.defaultdict(list)
    for r in SP: sp[(r["train"], r["setting"])].append(float(r["F1"]))
    A = supp_tsv("heldout_auc.tsv"); am = lambda k: st.mean(float(r["F1"]) for r in A if r["setting"] == k)
    put("three-pair mean minus four-pair mean, largest", f"{max(abs(round(st.mean(st.mean(v) for (t, s_), v in sp.items() if s_ == k), 3) - round(am(k), 3)) for k in SETTING + PURITY):.3f}")
    AB = supp_tsv("ablation.tsv"); abm = lambda v, k: st.mean(float(r["F1"]) for r in AB if r["variant"] == v and r["setting"] == k)
    put("ablation, largest change at full, 20x, 20% purity", f"{max(abs(round(abm(v, k), 3) - round(am(k), 3)) for v in {r['variant'] for r in AB} for k in ('full', '20x', '20p')):.3f}")
    g = Z["grid"]
    put("grid, mean F1 over all runs: parameters used; best set; set chosen inside the training pairs", f"{g['used']['all']:.3f}; {g['grid']['best']['all']:.3f}; {g['grid']['nested']['all']:.3f}")
    put("grid, mean F1 at full depth: parameters used; set chosen inside the training pairs", f"{g['used']['full']:.3f}; {g['grid']['nested']['full']:.3f}")
    ss = {}
    for i in range(1, 21):
        t = f"T{i:02d}"; sv = {}
        for l in open(f"{R}/sim/truth/{t}.tsv"):
            f = l.split("\t"); sv[f[0]] = (f[1], int(f[2]))
        s_ = {m for m, v in sv.items() if v[0] in ("dDUP", "rTRA") or v[1] >= 50}
        ss[t] = 100 * len({l.split()[0] for l in open(f"{R}/sim/truth/{t}.lost") if l.strip()} & s_) / len(s_)
    put("share of a tumor's SVs on a lost chromosome copy", rng(ss.values(), 0) + "%")
    fc = falsecalls("20-5")
    put("20/5 joint SVs, false calls, false calls in one tumor", f"{fc[0]:,}, {fc[1]:,}, {fc[2]:,}")
    put("simulated tree nodes", str(len(json.load(open(f"{R}/sim/truth/tree.json"))["parent"])))
    gd = f"{ROOT}/simulation/phylo/genomes"
    if os.path.exists(f"{gd}/germline/germline.truth.vcf.gz"):
        import gzip
        g = collections.Counter(); het = n = 0
        for l in gzip.open(f"{gd}/germline/germline.truth.vcf.gz", "rt"):
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t"); n += 1; het += f[9].split(":")[0] in ("0|1", "1|0")
            i = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            L = abs(int(i["SVLEN"])) if "SVLEN" in i else abs(len(f[3]) - len(f[4]))
            g["SNV" if "SVLEN" not in i and len(f[3]) == len(f[4]) == 1 else (i.get("SVTYPE", "INS" if len(f[4]) > len(f[3]) else "DEL") + " SV" if L >= 50 else "small indel")] += 1
        put("germline variants of the normal genome, heterozygous share", f"{g['DEL SV'] + g['INS SV']:,} SVs (" + ", ".join(f"{v:,} {k}" for k, v in sorted(g.items())) + f"); {het / n:.2f}")
        so = collections.Counter()
        for v in glob.glob(f"{gd}/tumor/*/mutations.vcf"):
            for l in open(v):
                if l[0] != "#": so[re.search(r"MUTCLASS=([^;\t]+)", l).group(1) if "MUTCLASS=" in l else "other"] += 1
        put("somatic variants of the tree by class", ", ".join(f"{v:,} {k}" for k, v in sorted(so.items())))
    with open(f"{TAB}/numbers.txt", "w") as o:
        for k, v in N: o.write(f"{k}\t{v}\n")
    return N


def falsecalls(name):
    g = {"__name__": "x", "SCRIPTS": paper.SCRIPTS, "run_score": paper.run_score}
    exec(paper.SCRIPTS["simjoint"].split("out = {}")[0], g)
    ns = g["ns"]; TI = [ns["index"](ns["truth"](f"{R}/sim/truth/{t}.tsv")[1]) for t in SIMT]
    d = pickle.load(open(f"{R}/sim/settings/{name.replace('-', '_')}/joint/carry.pkl", "rb"))
    pr, du, own, NB, J = [d[k] for k in ("pr", "du", "own", "NB", "J")]
    N = len(pr); cut = np.broadcast_to(np.asarray(d["cut"], float), (N,))
    conf = np.array([len({own[j] for j in NB.get(i, ()) if own[j] != own[i] and pr[j] >= cut[j]}) for i in range(N)])
    acc = ~du & ((pr >= cut + 0.2) | ((pr >= 0.1) & (conf >= 1)))
    root = list(range(N))

    def find(i):
        while root[i] != i:
            root[i] = root[root[i]]
            i = root[i]
        return i

    for i in range(N):
        for j in NB.get(i, ()):
            root[find(i)] = find(j)
    comp = collections.defaultdict(list)
    for i in np.where(acc)[0]:
        comp[find(i)].append(i)
    bad = [c for c, v in comp.items() if not any(J[i] and ns["match"](J[i], TI[own[i]]) for i in v)]
    return len(comp), len(bad), sum(1 for c in bad if len({own[i] for i in comp[c]}) == 1)


VAL = {}
FIGS = [fig2, fig3, fig4, datasets, supp_allsets, consensus, missedfig, supp_without, supp_depth, grid, filesfig, simgen, simeval, align, signatures, featfig, evidencefig, clsfig, jointfig, sweepfig, aucfig, splitsfig, paramfig, specfig]

if __name__ == "__main__":
    for f in FIGS:
        f()
    T, Z = tables()
    numbers(Z)
