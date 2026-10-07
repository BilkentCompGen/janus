#!/usr/bin/env python3
import zipfile, re, xml.etree.ElementTree as E
T = "/home/yuio/snape/hifi/truth"
HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

def xlsx_rows(path):
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(path)
    ss = ["".join(t.text or "" for t in si.iter(ns + "t")) for si in E.fromstring(z.read("xl/sharedStrings.xml")).findall(ns + "si")]
    for r in E.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(ns + "row"):
        out = []
        for c in r.findall(ns + "c"):
            v = c.find(ns + "v")
            out.append(ss[int(v.text)] if c.get("t") == "s" and v is not None else (v.text if v is not None else ""))
        yield out

n = 0
with open(f"{T}/hcc1395/truth.vcf", "w") as o:
    o.write(HEAD)
    for i, r in enumerate(xlsx_rows(f"{T}/hcc1395/13059_2022_2816_MOESM4_ESM.xlsx")):
        if i == 0 or len(r) < 7 or not r[2]: continue
        sid, c1, p1, c2, p2, size, ty = r[0], r[1], int(r[2]), r[3], int(r[4]), int(float(r[5] or 0)), r[6]
        if ty == "TRA":
            o.write(f"{c1}\t{p1}\t{sid}\tN\tN[{c2}:{p2}[\t.\tPASS\tSVTYPE=BND\n"); n += 1
        elif abs(size) >= 50 or ty in ("INV", "DUP") and p2 - p1 >= 50:
            end = p1 if ty == "INS" else p2
            o.write(f"{c1}\t{p1}\t{sid}\tN\t<{ty}>\t.\tPASS\tSVTYPE={ty};END={end};SVLEN={abs(size) if ty == 'INS' else p2 - p1}\n"); n += 1
print("hcc1395", n, "events")

n = 0
with open(f"{T}/colo829/truth.vcf", "w") as o:
    for l in open(f"{T}/colo829/truthset_somaticSVs_COLO829_hg38lifted.vcf"):
        if l.startswith("##"): o.write(l); continue
        if l.startswith("#"): o.write(l); continue
        f = l.rstrip("\n").split("\t")
        f[0] = "chr" + f[0] if not f[0].startswith("chr") else f[0]
        f[4] = re.sub(r"([\[\]])([^:\[\]]+):", lambda m: m.group(1) + (m.group(2) if m.group(2).startswith("chr") else "chr" + m.group(2)) + ":", f[4])
        o.write("\t".join(f) + "\n"); n += 1
print("colo829", n, "records")
