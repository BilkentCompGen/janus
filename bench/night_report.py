import glob, json, os, subprocess, sys
S = "/home/yuio/snape"; R = f"{S}/bench/runs/competitors"; T = f"{S}/hifi/truth"
PAIRS = ["hcc1954", "h2009", "hcc1937", "h1437", "hcc1395"]

def score(truth, vcf):
    o = subprocess.run(["python3", f"{S}/bench/score.py", truth, vcf], capture_output=True, text=True).stdout.split()
    return dict(zip(o[0::2], o[1::2])) if len(o) >= 10 else {}

rows = []
for p in PAIRS:
    v = f"/tmp/snape_{p}_pass.vcf"
    with open(v, "w") as o:
        o.writelines(l for l in open(f"{S}/bench/verdict/{p}.vcf") if l[0] == "#" or "\tPASS\t" in l)
    j = json.load(open(f"{S}/bench/runs/full/{p}/run.json"))["total"]
    rows.append(("snape", p, 0, j["wall_s"], j["peak_rss_kib"] / 2**20, "clean", v))
for d in sorted(glob.glob(f"{R}/*/*/run.json")):
    j = json.load(open(d)); t = j["total"]; D = os.path.dirname(d)
    if j["params"].get("pair") not in PAIRS: continue
    rows.append((j["params"]["tool"], j["params"]["pair"], t["exit"], t["wall_s"], t["peak_rss_kib"] / 2**20, j["params"].get("run", "clean"), f"{D}/calls.vcf"))
print("tool\tpair\texit\twall_s\tpeak_gb\trun\tcalls\tP\tR\tF1\tP_ws\tR_ws\tF1_ws")
for tool, p, ex, wall, gb, kind, vcf in sorted(rows, key=lambda r: (PAIRS.index(r[1]), r[0] != "snape", r[0])):
    a = b = {}
    if ex == 0 and os.path.exists(vcf):
        a = score(f"{T}/castle_without_snape/{p}/None_minda_ensemble.vcf", vcf)
        b = score(f"{T}/castle_without_severus/{p}/None_minda_ensemble.vcf", vcf)
    print(f"{tool}\t{p}\t{ex}\t{wall:.0f}\t{gb:.1f}\t{kind}\t{a.get('calls', '-')}\t{a.get('precision', '-')}\t{a.get('recall', '-')}\t{a.get('f1', '-')}\t{b.get('precision', '-')}\t{b.get('recall', '-')}\t{b.get('f1', '-')}")
