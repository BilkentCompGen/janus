#!/usr/bin/env python3
# truth.py: writes the truth sets that are not VCFs with chr names already, as hifi/truth/<pair>/truth.vcf (GRCh38, chr
# names, PASS, SVs of 50 bp and more; BND as bracket ALTs), so score.py and verdict.py read every truth set the same way:
#   hcc1395  SEQC2 high-confidence consensus, Talsania et al. 2022, Additional file 4 (Table S4, an .xlsx: read with the
#            standard library, it is zipped XML)
#   colo829  Valle-Inclan et al. 2022 truth set, the GRCh38 (hg38-lifted) VCF, chromosome names given the chr prefix
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

# ---- HCC1395: Table S4 rows: SV_id, Chrom1, Pos1, Chrom2, Pos2, size, type (INS DEL DUP INV TRA), ...
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

# ---- COLO829: chr prefix on CHROM and on the mate chromosome inside bracket ALTs
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
