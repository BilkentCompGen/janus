#!/usr/bin/env python3
"""Orchestrate the variant<->graph trace report.

Consumes the tracer outputs in <trace_dir> (variants.tsv, bubbles.tsv,
error_clusters.tsv, summary.json, annot.bin) and:
  seeds   -> pick ranked neighborhood seeds (all SVs + top bubbles/errors/SNPs) -> seeds.tsv
  neigh   -> ./tools/trace neigh ... (per-seed subgraph JSON + GFA + colours.csv)
  bandage -> render a curated gallery of subgraph PNGs (Bandage, headless)
  data    -> report_data.json consumed by report.html
  all     -> seeds, neigh, bandage, data

Run from the svp repo root.
"""
import os, sys, json, subprocess, collections

ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRAPH     = os.path.join(ROOT, "outputs/svp/refactor_e2e/graph.bin")
TRACE_BIN = os.path.join(ROOT, "tools/trace")
TRACE_DIR = os.path.join(ROOT, "outputs/svp/trace")
NEIGH_DIR = os.path.join(TRACE_DIR, "neigh")
RADIUS    = 30

# how many of each class get a precomputed subgraph (JSON, in-browser viewer)
N_SNP_VIEW, N_BUB_CLEAN, N_BUB_CONTAM, N_BUB_MATCHED, N_ERR_VIEW = 4000, 4000, 4000, 3000, 2000
# how many get a high-fidelity Bandage PNG (a curated headline gallery)
BANDAGE = dict(INV=99, DUP=99, DEL=150, INS=150, SNP=60, bub_clean=200, bub_contam=200, cycle=120, err=120)
# how many rows ship to the browser tables (variants+errors: all; bubbles: ranked subset)
N_BUB_TABLE = 80000

def read_tsv(path):
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(head, line.rstrip("\n").split("\t")))

# ---------------------------------------------------------------- seeds
def pick_seeds():
    variants = list(read_tsv(os.path.join(TRACE_DIR, "variants.tsv")))
    bubbles  = list(read_tsv(os.path.join(TRACE_DIR, "bubbles.tsv")))
    errors   = list(read_tsv(os.path.join(TRACE_DIR, "error_clusters.tsv")))
    seeds = {}                                   # label -> node_id (dedup by label)
    # variants: every found SV, plus the top-coverage found SNPs
    snps = []
    for v in variants:
        if v["found"] != "1":
            continue
        if v["type"] == "SNP":
            snps.append(v)
        else:
            seeds[f"var_{v['svid']}"] = v["best_node"]
    snps.sort(key=lambda v: int(v["best_tcov"]), reverse=True)
    for v in snps[:N_SNP_VIEW]:
        seeds[f"var_{v['svid']}"] = v["best_node"]
    # bubbles: top clean / contaminated by coverage, all matched, all cycles
    def cov(b): return int(b["ncov"]) + int(b["tcov"])
    clean   = sorted((b for b in bubbles if b["clean"] == "1" and b["kind"] == "bubble"), key=cov, reverse=True)
    contam  = sorted((b for b in bubbles if b["clean"] == "0" and b["kind"] == "bubble"), key=cov, reverse=True)
    matched = sorted((b for b in bubbles if b["matched_svid"] != "NA"), key=cov, reverse=True)
    cycles  = [b for b in bubbles if b["kind"] == "cycle"]
    for grp, n in ((clean, N_BUB_CLEAN), (contam, N_BUB_CONTAM), (matched, N_BUB_MATCHED), (cycles, len(cycles))):
        for b in grp[:n]:
            seeds[f"bub_{b['id']}"] = b["s1_node"]
    # errors: top clusters by size (errors.tsv already size-sorted)
    for e in errors[:N_ERR_VIEW]:
        seeds[f"err_{e['cluster']}"] = e["seed_node"]
    os.makedirs(TRACE_DIR, exist_ok=True)
    with open(os.path.join(TRACE_DIR, "seeds.tsv"), "w") as f:
        f.write("#label\tnode\n")
        for label, node in seeds.items():
            f.write(f"{label}\t{node}\n")
    print(f"[report] picked {len(seeds)} seeds -> seeds.tsv")

# ---------------------------------------------------------------- neigh
def run_neigh():
    subprocess.run([TRACE_BIN, "neigh", GRAPH, TRACE_DIR,
                    os.path.join(TRACE_DIR, "seeds.tsv"), str(RADIUS)], cwd=ROOT, check=True)

# ---------------------------------------------------------------- bandage gallery
def bandage_one(label):
    gfa = os.path.join(NEIGH_DIR, label + ".gfa")
    png = os.path.join(NEIGH_DIR, label + ".png")
    csv = os.path.join(NEIGH_DIR, label + ".colors.csv")
    if not os.path.exists(gfa):
        return False
    if os.path.exists(png) and os.path.getmtime(png) >= os.path.getmtime(gfa):
        return True                                   # fresh (PNG newer than its GFA) -> skip
    env = dict(os.environ, QT_QPA_PLATFORM="offscreen")
    cmd = ["Bandage", "image", gfa, png, "--width", "900", "--height", "650"]
    if os.path.exists(csv):
        cmd += ["--colour", "custom", "--colors", csv]
    try:
        subprocess.run(cmd, env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=60)
        return os.path.exists(png)
    except Exception:
        return False

def gallery_labels():
    variants = list(read_tsv(os.path.join(TRACE_DIR, "variants.tsv")))
    bubbles  = list(read_tsv(os.path.join(TRACE_DIR, "bubbles.tsv")))
    errors   = list(read_tsv(os.path.join(TRACE_DIR, "error_clusters.tsv")))
    labels = []
    byt = collections.defaultdict(list)
    for v in variants:
        if v["found"] == "1":
            byt[v["type"]].append(v)
    for t in ("INV", "DUP", "DEL", "INS", "SNP"):
        rows = sorted(byt[t], key=lambda v: int(v["best_tcov"]), reverse=True)[:BANDAGE[t]]
        labels += [f"var_{v['svid']}" for v in rows]
    def cov(b): return int(b["ncov"]) + int(b["tcov"])
    clean  = sorted((b for b in bubbles if b["clean"] == "1" and b["kind"] == "bubble"), key=cov, reverse=True)[:BANDAGE["bub_clean"]]
    contam = sorted((b for b in bubbles if b["clean"] == "0" and b["kind"] == "bubble"), key=cov, reverse=True)[:BANDAGE["bub_contam"]]
    cycles = [b for b in bubbles if b["kind"] == "cycle"][:BANDAGE["cycle"]]
    labels += [f"bub_{b['id']}" for b in clean + contam + cycles]
    labels += [f"err_{e['cluster']}" for e in errors[:BANDAGE["err"]]]
    return labels

def run_bandage():
    labels = gallery_labels()
    from concurrent.futures import ThreadPoolExecutor
    done = 0
    with ThreadPoolExecutor(max_workers=8) as ex:
        for ok in ex.map(bandage_one, labels):
            done += 1 if ok else 0
            if done % 50 == 0:
                print(f"[report] bandage {done}/{len(labels)}", flush=True)
    print(f"[report] bandage gallery: {done}/{len(labels)} PNGs -> {NEIGH_DIR}")

# ---------------------------------------------------------------- report_data.json
def has_view(label): return os.path.exists(os.path.join(NEIGH_DIR, label + ".json"))
def has_png(label):
    png, gfa = os.path.join(NEIGH_DIR, label + ".png"), os.path.join(NEIGH_DIR, label + ".gfa")
    return os.path.exists(png) and os.path.exists(gfa) and os.path.getmtime(png) >= os.path.getmtime(gfa)

def build_data():
    summary = json.load(open(os.path.join(TRACE_DIR, "summary.json")))
    variants = list(read_tsv(os.path.join(TRACE_DIR, "variants.tsv")))
    bubbles  = list(read_tsv(os.path.join(TRACE_DIR, "bubbles.tsv")))
    errors   = list(read_tsv(os.path.join(TRACE_DIR, "error_clusters.tsv")))

    # aggregates
    recall = collections.defaultdict(lambda: [0, 0])
    coloragree = [0, 0]
    for v in variants:
        recall[v["type"]][1] += 1
        if v["found"] == "1":
            recall[v["type"]][0] += 1
            coloragree[1] += 1
            if v["best_color"] == "2":
                coloragree[0] += 1
    bub_agg = dict(total=len(bubbles), clean=0, contam=0, cycles=0, matched=0,
                   buckets=collections.Counter())
    for b in bubbles:
        if b["kind"] == "cycle": bub_agg["cycles"] += 1
        bub_agg["clean" if b["clean"] == "1" else "contam"] += 1
        if b["matched_svid"] != "NA": bub_agg["matched"] += 1
        bub_agg["buckets"][b["bucket"]] += 1
    err_hist = collections.Counter(); err_loc = 0
    for e in errors:
        s = int(e["size"]); err_hist["1" if s==1 else "2" if s==2 else "3-5" if s<=5 else "6-20" if s<=20 else ">20"] += 1
        if e["near_chrom"] != "NA": err_loc += 1

    # ---- table columns (compact) ----
    def col(rows, keys): return {k: [r[k] for r in rows] for k in keys}
    for v in variants: v["v"] = 1 if has_view(f"var_{v['svid']}") else 0; v["p"] = 1 if has_png(f"var_{v['svid']}") else 0
    vt = col(variants, ["svid","chrom","pos1","type","svlen","n_nodes","best_color","best_tcov","best_ncov","found","v","p"])

    # bubbles: ALL of them (matched first, then contaminated, then clean by coverage), shipped
    # as a separate lazily-loaded file so the main page stays light.
    def bkey(b): return (0 if b["matched_svid"]!="NA" else 1, 0 if b["clean"]=="0" else 1, -(int(b["ncov"])+int(b["tcov"])))
    bubbles.sort(key=bkey)
    for b in bubbles: b["v"] = 1 if has_view(f"bub_{b['id']}") else 0; b["p"] = 1 if has_png(f"bub_{b['id']}") else 0
    bt = col(bubbles, ["id","kind","bucket","chrom","pos","arm_members","clean","ncov","tcov","matched_svid","v","p"])
    json.dump(dict(bub=bt), open(os.path.join(TRACE_DIR, "bubbles.json"), "w"))

    for e in errors: e["v"] = 1 if has_view(f"err_{e['cluster']}") else 0; e["p"] = 1 if has_png(f"err_{e['cluster']}") else 0
    et = col(errors, ["cluster","size","n_c1","n_c2","n_c3","seed_node","near_chrom","near_pos","v","p"])

    data = dict(
        summary=summary,
        recall_by_type={k: recall[k] for k in recall},
        color_agreement=coloragree,
        bubbles=dict(total=bub_agg["total"], clean=bub_agg["clean"], contam=bub_agg["contam"],
                     cycles=bub_agg["cycles"], matched=bub_agg["matched"], buckets=dict(bub_agg["buckets"])),
        errors=dict(nodes=summary["read_error"], clusters=len(errors), located=err_loc, hist=dict(err_hist)),
        variants_total=len(variants), bubbles_shown=len(bubbles), errors_total=len(errors),
        var=vt, err=et,                                   # bub shipped separately (bubbles.json)
    )
    out = os.path.join(TRACE_DIR, "report_data.json")
    json.dump(data, open(out, "w"))
    print(f"[report] wrote {out} ({os.path.getsize(out)/1e6:.1f} MB); "
          f"views={sum(vt['v'])+sum(bt['v'])+sum(et['v'])} pngs={sum(vt['p'])+sum(bt['p'])+sum(et['p'])}")

# ----------------------------------------------------------------
if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    if cmd in ("seeds", "all"):   pick_seeds()
    if cmd in ("neigh", "all"):   run_neigh()
    if cmd in ("bandage", "all"): run_bandage()
    if cmd in ("data", "all"):    build_data()
