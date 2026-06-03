#!/usr/bin/env python3
"""inspect_viz — Stage 3 renderer for the svp graph inspector.

Reads <dir>/inspect/<key>/{raw,modeA,modeB}.gfa + meta.json. Emits:

  topology.png         matplotlib BFS-depth + color-stripe layout (uses modeA by default)
  linear.png           matplotlib 2-track linear layout from modeB
  bandage.png          (if bandage-ng on PATH) Bandage's render of modeA
  odgi.png             (if odgi on PATH) ODGI's 1D viz of raw
  graph.dot            Graphviz DOT export of modeA
  index.html           per-key report

Usage:  python3 inspect_viz.py <key_dir>
        python3 inspect_viz.py <key_dir1> <key_dir2> ...   (batch mode + aggregate report)
"""

from __future__ import annotations

import collections
import json
import os
import shutil
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ---------- GFA1 parser (lightweight) ----------
def parse_gfa(path: str) -> Tuple[Dict[str, dict], List[Tuple[str, str]]]:
    nodes: Dict[str, dict] = {}
    edges: List[Tuple[str, str]] = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line[0] == "H":
                continue
            t = line.rstrip("\n").split("\t")
            if t[0] == "S":
                name = t[1]
                seq = t[2]
                node = {"name": name, "seq": seq, "tags": {}}
                for tag in t[3:]:
                    parts = tag.split(":", 2)
                    if len(parts) == 3:
                        node["tags"][parts[0]] = (parts[1], parts[2])
                nodes[name] = node
            elif t[0] == "L":
                edges.append((t[1], t[3]))
    return nodes, edges


def tag_str(node: dict, key: str, default: str = "") -> str:
    if key in node["tags"]:
        return node["tags"][key][1]
    return default


def tag_int(node: dict, key: str, default: int = 0) -> int:
    if key in node["tags"]:
        return int(node["tags"][key][1])
    return default


# ---------- Topology layout (BFS depth from anchors) ----------
def color_for_label(lbl: str) -> str:
    return {"NORMAL": "#5B9BD5", "TUMOR": "#E07B5C", "SHARED": "#8E8E8E"}.get(lbl, "#000000")


def topology_layout(nodes, edges):
    # Build adjacency
    adj = collections.defaultdict(set)
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)

    # Anchor nodes = CL=SHARED
    anchors = [n for n, d in nodes.items() if tag_str(d, "CL") == "SHARED"]
    if not anchors:
        # Fallback: rightmost / leftmost by name
        anchors = sorted(nodes.keys())[:2]

    # BFS depth from any anchor
    depth: Dict[str, int] = {a: 0 for a in anchors}
    q = collections.deque(anchors)
    while q:
        u = q.popleft()
        for v in adj[u]:
            if v not in depth:
                depth[v] = depth[u] + 1
                q.append(v)
    # Disconnected: any node not reached gets max-depth+1
    max_d = max(depth.values()) if depth else 0
    for n in nodes:
        if n not in depth:
            depth[n] = max_d + 1

    # Identify dead-ends: a non-shared node that BFS reaches from one anchor side
    # but doesn't form a path to the OTHER anchor.  Simple heuristic: a node that has
    # only one path back to the anchor set.  We approximate: degree-1 non-shared nodes
    # whose only neighbor is also non-shared.
    deadends = set()
    for n, d in nodes.items():
        if tag_str(d, "CL") == "SHARED":
            continue
        nbrs = adj[n]
        if len(nbrs) == 1:
            deadends.add(n)

    # x = BFS depth; y = color stripe (NORMAL=+1, TUMOR=-1, SHARED=0) + small per-node offset
    pos = {}
    by_depth_color: Dict[Tuple[int, str], List[str]] = collections.defaultdict(list)
    for n in sorted(nodes.keys()):
        cl = tag_str(nodes[n], "CL")
        by_depth_color[(depth[n], cl)].append(n)
    for (d, cl), ns in by_depth_color.items():
        ystripe = {"NORMAL": 1.0, "TUMOR": -1.0, "SHARED": 0.0}.get(cl, 0.0)
        for j, n in enumerate(sorted(ns)):
            # Spread within stripe
            offset = (j - (len(ns) - 1) / 2.0) * 0.18
            pos[n] = (d * 1.0, ystripe + offset)
    return pos, deadends


def plot_topology(nodes, edges, pos, deadends, out_path: str, meta: dict):
    fig, ax = plt.subplots(figsize=(13, 6), dpi=160)
    # Edges
    for a, b in edges:
        if a not in pos or b not in pos:
            continue
        xs, ys = [pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]]
        # Edge color: if both endpoints have same CL, use that color; else grey
        ca = tag_str(nodes[a], "CL")
        cb = tag_str(nodes[b], "CL")
        ec = color_for_label(ca) if ca == cb else "#bbbbbb"
        ax.plot(xs, ys, "-", color=ec, alpha=0.6, lw=1.2, zorder=1)

    # Nodes
    for n, d in nodes.items():
        if n not in pos:
            continue
        x, y = pos[n]
        cl = tag_str(d, "CL")
        ln = tag_int(d, "LN")
        size = max(50, min(800, 8 * ln))
        edgecolor = "#cc0000" if n in deadends else "#222222"
        ax.scatter([x], [y], s=size, c=color_for_label(cl),
                   edgecolors=edgecolor, linewidths=1.5 if n in deadends else 0.8,
                   zorder=2)
        # Label
        ax.text(x, y - 0.18, f"{n}\nLN={ln}", fontsize=6, ha="center", va="top",
                color="#444444", zorder=3)

    ax.set_xlabel("BFS depth from anchors")
    ax.set_yticks([-1, 0, 1])
    ax.set_yticklabels(["TUMOR", "SHARED", "NORMAL"])
    ax.grid(True, alpha=0.2, axis="y")

    # Annotations
    title = meta.get("key", "?")
    if "truth_sv" in meta:
        tsv = meta["truth_sv"]
        title += f"  |  TRUTH {tsv.get('svtype','?')} {tsv.get('svlen',0)} bp @ {tsv.get('chrom','?')}:{tsv.get('pos',0)}"
    ax.set_title(title, fontsize=10)

    # Legend
    handles = [
        mpatches.Patch(color=color_for_label("NORMAL"), label="NORMAL (color=1)"),
        mpatches.Patch(color=color_for_label("TUMOR"),  label="TUMOR (color=2)"),
        mpatches.Patch(color=color_for_label("SHARED"), label="SHARED (color=3 anchor)"),
        mpatches.Patch(color="#cc0000", label="dead-end (only 1 neighbor)"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, framealpha=0.9)

    # Per-bubble SV tag overlays
    if "modeB_calls" in meta and meta["modeB_calls"]:
        info = "Mode B classification:\n"
        for c in meta["modeB_calls"]:
            info += f"  bubble {c.get('bubble_id','?')}: SV={c.get('sv_tag','?')}  dN={c.get('diff_n','?')}  dT={c.get('diff_t','?')}\n"
        ax.text(0.01, 0.01, info, transform=ax.transAxes, fontsize=7,
                verticalalignment="bottom", family="monospace",
                bbox=dict(facecolor="#f7f7f7", edgecolor="#888888", boxstyle="round,pad=0.3"))

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(fig)


# ---------- Linear layout (Mode B only) ----------
def plot_linear(nodes, edges, out_path: str, meta: dict):
    # Identify L_anchor and R_anchor: the two CL=SHARED nodes with PS:Z, by leftmost/rightmost pos
    anchors = []
    for n, d in nodes.items():
        if tag_str(d, "CL") == "SHARED":
            ps = tag_str(d, "PS")
            if ps:
                # PS:Z:<chrom>:<start>-<end>:<strand>
                anchors.append((n, ps))
    anchors.sort(key=lambda x: int(x[1].split(":")[1].split("-")[0]) if ":" in x[1] else 0)

    # Identify N-segment and T-segment per bubble
    n_segs = [n for n, d in nodes.items() if tag_str(d, "CL") == "NORMAL" and d["seq"] != "*" and d["seq"]]
    t_segs = [n for n, d in nodes.items() if tag_str(d, "CL") == "TUMOR"  and d["seq"] != "*" and d["seq"]]

    fig, ax = plt.subplots(figsize=(13, 4.5), dpi=160)

    # Compute layout: anchors at fixed positions, diff segments between them
    y_normal = 1.0
    y_tumor  = -1.0
    track_height = 0.45
    pad = 0.5

    # If we have at least one L+R anchor, lay them out at x=0 and x=END
    if len(anchors) >= 2:
        L = anchors[0]
        R = anchors[-1]
        L_len = tag_int(nodes[L[0]], "LN")
        R_len = tag_int(nodes[R[0]], "LN")
        # Total width = L_len + max(N_diff, T_diff, 1) + R_len + small padding
        N_diff = sum(tag_int(nodes[s], "LN") for s in n_segs)
        T_diff = sum(tag_int(nodes[s], "LN") for s in t_segs)
        mid_len = max(N_diff, T_diff, 1)

        # Anchors share both tracks (normal and tumor — same content)
        for y, lbl in [(y_normal, "NORMAL"), (y_tumor, "TUMOR")]:
            ax.add_patch(mpatches.Rectangle((0, y - track_height/2), L_len, track_height,
                                            facecolor="#8E8E8E", edgecolor="#444444"))
            ax.text(L_len/2, y, f"{L[0]}\nLN={L_len}", ha="center", va="center", fontsize=7)
            ax.add_patch(mpatches.Rectangle((L_len + mid_len, y - track_height/2), R_len, track_height,
                                            facecolor="#8E8E8E", edgecolor="#444444"))
            ax.text(L_len + mid_len + R_len/2, y, f"{R[0]}\nLN={R_len}", ha="center", va="center", fontsize=7)

        # Normal diff
        if n_segs:
            x = L_len
            for s in n_segs:
                w = tag_int(nodes[s], "LN")
                ax.add_patch(mpatches.Rectangle((x, y_normal - track_height/2), w, track_height,
                                                facecolor=color_for_label("NORMAL"), edgecolor="#222"))
                ax.text(x + w/2, y_normal, f"{s}\n{w} bp", ha="center", va="center", fontsize=7, color="white")
                x += w
        else:
            ax.text(L_len + mid_len/2, y_normal, "(collapsed)", ha="center", va="center", fontsize=8, color="#666666")

        if t_segs:
            x = L_len
            for s in t_segs:
                w = tag_int(nodes[s], "LN")
                ax.add_patch(mpatches.Rectangle((x, y_tumor - track_height/2), w, track_height,
                                                facecolor=color_for_label("TUMOR"), edgecolor="#222"))
                ax.text(x + w/2, y_tumor, f"{s}\n{w} bp", ha="center", va="center", fontsize=7, color="white")
                x += w
        else:
            ax.text(L_len + mid_len/2, y_tumor, "(collapsed)", ha="center", va="center", fontsize=8, color="#666666")

        ax.set_xlim(-pad, L_len + mid_len + R_len + pad)
        ax.set_ylim(-2.5, 2.5)
        # Coord ruler
        ps_L = anchors[0][1]
        chrom_pos = ps_L.split(":")[0]
        start_pos = int(ps_L.split(":")[1].split("-")[0])
        ax.text(0, 2.0, f"{chrom_pos}:{start_pos}", fontsize=7, ha="left")

    ax.set_yticks([y_tumor, 0, y_normal])
    ax.set_yticklabels(["TUMOR", "", "NORMAL"])
    ax.set_xlabel("bp")
    ax.grid(True, alpha=0.2, axis="x")

    title = meta.get("key", "?") + " — Mode B linear"
    if "modeB_calls" in meta and meta["modeB_calls"]:
        c = meta["modeB_calls"][0]
        title += f"  |  SV={c.get('sv_tag','?')}  dN={c.get('diff_n','?')}  dT={c.get('diff_t','?')}"
    ax.set_title(title, fontsize=10)

    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(fig)


# ---------- DOT export ----------
def emit_dot(nodes, edges, out_path: str):
    lines = ["digraph G {", '  rankdir="LR";', '  node [shape=box, style=filled];']
    for n, d in nodes.items():
        cl = tag_str(d, "CL")
        fill = color_for_label(cl)
        ln = tag_int(d, "LN")
        lines.append(f'  "{n}" [label="{n}\\nLN={ln}", fillcolor="{fill}"];')
    for a, b in edges:
        lines.append(f'  "{a}" -> "{b}";')
    lines.append("}")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ---------- External tool subprocess wrappers ----------
def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def run_bandage(gfa: str, png: str):
    if not have("bandage-ng"):
        return False
    try:
        subprocess.check_call(["bandage-ng", "image", gfa, png],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return os.path.exists(png)
    except subprocess.CalledProcessError:
        return False


def run_odgi(gfa: str, png: str):
    if not have("odgi"):
        return False
    try:
        og = "/tmp/insp_viz.og"
        sg = "/tmp/insp_viz_sorted.og"
        subprocess.check_call(["odgi", "build", "-g", gfa, "-o", og],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.check_call(["odgi", "sort", "-i", og, "-o", sg],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.check_call(["odgi", "viz", "-i", sg, "-o", png, "-x", "1500", "-y", "300"],
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return os.path.exists(png)
    except subprocess.CalledProcessError:
        return False


# ---------- Per-key HTML ----------
def emit_index(key_dir: str, files: List[str], meta: dict):
    title = meta.get("key", os.path.basename(key_dir))
    lines = [
        "<!doctype html>",
        f"<html><head><title>{title}</title>",
        "<style>",
        "body { font-family: -apple-system, BlinkMacSystemFont, sans-serif; margin: 1em; }",
        "img { max-width: 100%; border: 1px solid #ccc; margin: 1em 0; }",
        "pre { background: #f7f7f7; padding: 1em; border: 1px solid #ddd; overflow-x: auto; }",
        "</style></head><body>",
        f"<h1>{title}</h1>",
    ]
    if "truth_sv" in meta:
        tsv = meta["truth_sv"]
        lines.append(f"<p><b>Truth:</b> {tsv.get('svtype','?')} {tsv.get('svlen',0)} bp at {tsv.get('chrom','?')}:{tsv.get('pos',0)}-{tsv.get('end',0)}</p>")
    for fname in files:
        if fname.endswith(".png"):
            lines.append(f"<h2>{fname}</h2><img src=\"{fname}\">")
    lines.append("<h2>meta.json</h2>")
    lines.append(f"<pre>{json.dumps(meta, indent=2)}</pre>")
    lines.append("<h2>Files</h2><ul>")
    for fname in sorted(os.listdir(key_dir)):
        lines.append(f'<li><a href="{fname}">{fname}</a></li>')
    lines.append("</ul></body></html>")
    with open(os.path.join(key_dir, "index.html"), "w") as f:
        f.write("\n".join(lines))


# ---------- Aggregate ----------
def emit_aggregate(parent_dir: str, keys: List[str]):
    out = os.path.join(parent_dir, "index.html")
    lines = [
        "<!doctype html>",
        "<html><head><title>svp inspect — aggregate</title>",
        "<style>",
        "body { font-family: -apple-system, BlinkMacSystemFont, sans-serif; margin: 1em; }",
        "table { border-collapse: collapse; }",
        "td, th { border: 1px solid #ccc; padding: 0.3em 0.6em; }",
        "img.thumb { max-height: 100px; }",
        "</style></head><body>",
        "<h1>svp inspect — all inspected subgraphs</h1>",
        "<table><tr><th>key</th><th>truth</th><th>topology</th><th>linear</th></tr>",
    ]
    for k in keys:
        meta_path = os.path.join(parent_dir, k, "meta.json")
        meta = {}
        if os.path.exists(meta_path):
            try:
                with open(meta_path) as f: meta = json.load(f)
            except Exception:
                pass
        truth = ""
        if "truth_sv" in meta:
            tsv = meta["truth_sv"]
            truth = f"{tsv.get('svtype','?')} {tsv.get('svlen',0)}"
        top_png = os.path.join(k, "topology.png")
        lin_png = os.path.join(k, "linear.png")
        lines.append(
            f"<tr><td><a href='{k}/index.html'>{k}</a></td>"
            f"<td>{truth}</td>"
            f"<td>{'<img class=thumb src=\"'+top_png+'\">' if os.path.exists(os.path.join(parent_dir, top_png)) else '—'}</td>"
            f"<td>{'<img class=thumb src=\"'+lin_png+'\">' if os.path.exists(os.path.join(parent_dir, lin_png)) else '—'}</td>"
            "</tr>"
        )
    lines.append("</table></body></html>")
    with open(out, "w") as f:
        f.write("\n".join(lines))
    print(f"[viz]  aggregate → {out}")


def render_key(key_dir: str):
    print(f"[viz] {key_dir}")
    meta_path = os.path.join(key_dir, "meta.json")
    meta: dict = {}
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)

    raw_gfa = os.path.join(key_dir, "raw.gfa")
    mA_gfa = os.path.join(key_dir, "modeA.gfa")
    mB_gfa = os.path.join(key_dir, "modeB.gfa")

    out_files: List[str] = []

    # Topology (uses Mode A by default)
    if os.path.exists(mA_gfa):
        n, e = parse_gfa(mA_gfa)
        pos, deadends = topology_layout(n, e)
        out = os.path.join(key_dir, "topology.png")
        plot_topology(n, e, pos, deadends, out, meta)
        out_files.append("topology.png")
        emit_dot(n, e, os.path.join(key_dir, "graph.dot"))

    # Linear (Mode B)
    if os.path.exists(mB_gfa):
        n, e = parse_gfa(mB_gfa)
        out = os.path.join(key_dir, "linear.png")
        plot_linear(n, e, out, meta)
        out_files.append("linear.png")

    # Optional Bandage-NG
    if os.path.exists(mA_gfa) and run_bandage(mA_gfa, os.path.join(key_dir, "bandage.png")):
        out_files.append("bandage.png")
    if os.path.exists(raw_gfa) and run_odgi(raw_gfa, os.path.join(key_dir, "odgi.png")):
        out_files.append("odgi.png")

    emit_index(key_dir, out_files, meta)


def main():
    if len(sys.argv) < 2:
        print("usage: inspect_viz.py <key_dir> [<key_dir2> ...]", file=sys.stderr)
        sys.exit(1)
    key_dirs = [os.path.abspath(p) for p in sys.argv[1:]]
    for kd in key_dirs:
        render_key(kd)
    # Aggregate report at parent of all keys (they should share a parent for this to be meaningful)
    parents = {os.path.dirname(p) for p in key_dirs}
    if len(parents) == 1:
        parent = parents.pop()
        # List all sibling key dirs in parent
        keys = sorted(d for d in os.listdir(parent) if os.path.isdir(os.path.join(parent, d)))
        emit_aggregate(parent, keys)


if __name__ == "__main__":
    main()
