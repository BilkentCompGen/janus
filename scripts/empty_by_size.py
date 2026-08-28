#!/usr/bin/env python3
# Stratify graph shapes (resolved graph) by divergent-middle size and shape type, showing the
# empty (no-variant) surface as (a) an empty-RATE heatmap and (b) an empty-MASS heatmap.
# Input: outputs/svp/trace_xxh3/{structures.tsv,shape_census.tsv}. Self-contained HTML, no CDN.
import sys, html

DIR = "outputs/svp/trace_xxh3"
STRUCT = f"{DIR}/structures.tsv"
CENSUS = f"{DIR}/shape_census.tsv"
OUT = f"{DIR}/empty_by_size.html"

BINS = [("0 (n/a)", lambda s: s == 0), ("1", lambda s: s == 1), ("2-10", lambda s: 2 <= s <= 10),
        ("11-25", lambda s: 11 <= s <= 25), ("26-50", lambda s: 26 <= s <= 50),
        ("51-100", lambda s: 51 <= s <= 100), ("101-200", lambda s: 101 <= s <= 200),
        ("201-300", lambda s: 201 <= s <= 300), ("301-500", lambda s: 301 <= s <= 500),
        ("501-1k", lambda s: 501 <= s <= 1000), (">1k", lambda s: s > 1000)]
SHAPES = ["simple", "asym", "seam", "super", "cycle"]

def binof(s):
    for name, f in BINS:
        if f(s): return name
    return ">1k"

# (bin, shape) -> [total, matched]
cell = {}
for b, _ in BINS:
    for sh in SHAPES:
        cell[(b, sh)] = [0, 0]

with open(STRUCT) as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        shp, l1, l2, matched = f[0], int(f[2]), int(f[3]), int(f[5])
        if shp not in SHAPES: continue
        b = binof(max(l1, l2))
        c = cell[(b, shp)]
        c[0] += 1
        c[1] += matched

# super-granularity shapes from census (no size dimension)
supers = []
with open(CENSUS) as fh:
    next(fh)
    for line in fh:
        sh, gran, cnt, m, e = line.rstrip("\n").split("\t")
        if gran == "super":
            supers.append((sh, int(cnt), int(m), int(e)))

# totals
tot_all = sum(cell[k][0] for k in cell) + sum(s[1] for s in supers)
mat_all = sum(cell[k][1] for k in cell) + sum(s[2] for s in supers)
emp_all = tot_all - mat_all

def rate_color(pct):  # green(0% empty)->amber(50)->red(100)
    if pct <= 50:
        t = pct / 50.0
        r, g, b = int(46 + t * (245 - 46)), int(160 + t * (200 - 160)), int(67 + t * (66 - 67))
    else:
        t = (pct - 50) / 50.0
        r, g, b = int(245 + t * (211 - 245)), int(200 + t * (47 - 200)), int(66 + t * (47 - 66))
    return f"rgb({r},{g},{b})"

def mass_color(cnt, mx):  # white->deep red on sqrt scale
    t = (cnt / mx) ** 0.5 if mx else 0
    r, g, b = 255, int(255 - t * (255 - 30)), int(255 - t * (255 - 30))
    return f"rgb({r},{g},{b})"

max_emp = max((cell[k][0] - cell[k][1]) for k in cell)

def grid_html(mode):
    rows = []
    hdr = "".join(f"<th>{sh}</th>" for sh in SHAPES) + "<th class='tot'>row Σ</th>"
    for b, _ in BINS:
        tds = []
        rt = rm = 0
        for sh in SHAPES:
            t, m = cell[(b, sh)]
            e = t - m
            rt += t; rm += m
            if t == 0:
                tds.append("<td class='na'>·</td>"); continue
            pct = 100.0 * e / t
            if mode == "rate":
                tds.append(f"<td style='background:{rate_color(pct)}' title='{e}/{t} empty'>"
                           f"<b>{pct:.0f}%</b><br><span class='sub'>{e}/{t}</span></td>")
            else:
                tds.append(f"<td style='background:{mass_color(e, max_emp)}' title='{pct:.0f}% empty'>"
                           f"<b>{e:,}</b><br><span class='sub'>{pct:.0f}%</span></td>")
        re = rt - rm
        rpct = 100.0 * re / rt if rt else 0
        if mode == "rate":
            tot = f"<td class='tot' style='background:{rate_color(rpct)}'><b>{rpct:.0f}%</b><br><span class='sub'>{re:,}/{rt:,}</span></td>"
        else:
            tot = f"<td class='tot' style='background:{mass_color(re, max_emp)}'><b>{re:,}</b><br><span class='sub'>{rpct:.0f}%</span></td>"
        rows.append(f"<tr><th class='rb'>{b}</th>{''.join(tds)}{tot}</tr>")
    return f"<table><tr><th>size (bp)</th>{hdr}</tr>{''.join(rows)}</table>"

sup_rows = "".join(
    f"<tr><th class='rb'>{html.escape(sh)}</th><td>{cnt:,}</td><td>{m:,}</td>"
    f"<td style='background:{rate_color(100*e/cnt)}'><b>{e:,}</b> ({100*e/cnt:.1f}%)</td></tr>"
    for sh, cnt, m, e in supers)

page = f"""<!doctype html><meta charset=utf-8>
<title>svp · empty shapes by size</title>
<style>
 body{{font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;margin:0;background:#0f1115;color:#e6e6e6}}
 .wrap{{max-width:1100px;margin:0 auto;padding:28px}}
 h1{{font-size:21px;margin:0 0 4px}} h2{{font-size:16px;margin:30px 0 6px;color:#cdd3da}}
 .sub2{{color:#8b95a1;font-size:13px;margin:0 0 18px}}
 .kpis{{display:flex;gap:14px;margin:16px 0 4px;flex-wrap:wrap}}
 .kpi{{background:#1a1e26;border:1px solid #2a2f3a;border-radius:10px;padding:12px 16px;flex:1;min-width:150px}}
 .kpi b{{display:block;font-size:24px;color:#fff}} .kpi span{{color:#8b95a1;font-size:12px}}
 table{{border-collapse:collapse;margin:6px 0 4px;width:100%}}
 th,td{{border:1px solid #11141a;padding:6px 4px;text-align:center;font-size:12.5px}}
 th{{background:#1a1e26;color:#cdd3da;font-weight:600}} td b{{color:#0b0d10}}
 .rb{{text-align:right;background:#1a1e26;color:#cdd3da;white-space:nowrap}}
 .sub{{font-size:10.5px;color:#1a1a1a;opacity:.75}} .na{{background:#161922;color:#3a4150}}
 .tot{{font-weight:700}} .cap{{color:#9aa4b0;font-size:12.5px;margin:6px 0 2px}}
 .small{{max-width:520px}}
</style>
<div class=wrap>
<h1>Empty graph shapes, stratified by size × shape</h1>
<p class=sub2>Resolved graph (X-cut applied) · a shape is <b>empty</b> when none of its somatic arm supers
 trace to a VCF variant · size = divergent-middle span max(len1,len2) in bp.</p>
<div class=kpis>
 <div class=kpi><b>{tot_all:,}</b><span>total shapes</span></div>
 <div class=kpi><b>{emp_all:,}</b><span>empty — {100*emp_all/tot_all:.1f}% trace to no variant</span></div>
 <div class=kpi><b>{mat_all:,}</b><span>back ≥1 variant</span></div>
</div>

<h2>A · Empty <i>rate</i> — discrimination signal (green = mostly real, red = mostly artifact)</h2>
<p class=cap>Flat ~50% down every bubble column → <b>size barely separates variant from error inside bubbles</b>.
 Seam stays 92–100% empty at every size → a one-sided shape is almost never a variant.</p>
{grid_html("rate")}

<h2>B · Empty <i>mass</i> — where the false-shape surface actually lives</h2>
<p class=cap>Cell = empty count. The surface is dominated by <b>seams at 100–500 bp</b> (unrejoined read-error
 tails) plus the <b>size-1 simple-bubble pile</b> (SNP vs error, 50/50).</p>
{grid_html("mass")}

<h2>C · Super-granularity shapes (no size) — the artifact pile</h2>
<table class=small><tr><th>shape</th><th>count</th><th>matched</th><th>empty</th></tr>{sup_rows}</table>
<p class=cap>Tips (dead-end errors), tangles (unresolved repeats) and orphans (isolated components) are
 99.4–100% empty — filtered for free by colour + coverage before any call.</p>
</div>"""

open(OUT, "w").write(page)
print(f"wrote {OUT}")
print(f"total={tot_all:,} empty={emp_all:,} ({100*emp_all/tot_all:.1f}%) matched={mat_all:,}")
