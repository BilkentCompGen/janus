#!/usr/bin/env python3
"""Gercek (VCF 4.2) bir varyant dosyasini insilicoSV'nin import formatina cevirir.

    vcf_to_insilicosv.py <in.vcf[.gz]> <out.vcf> [--contig-map chr20=chr20] [--sample S]
                         [--drop-overlaps] [--report out.tsv]

NEDEN GEREKLI
-------------
insilicoSV'nin import yolu (`insilicosv/variant_set.py:1028-1075`) sunlari sart kosuyor:

  * Her kaydin **SVTYPE veya OP_TYPE**'i olmali. Gercek VCF'lerdeki sade indel kayitlari
    (REF=TG ALT=T) SVTYPE tasimaz -> **reddedilir**. Bu yuzden acikca etiketliyoruz.
  * SNP otomatik algilamasi yalnizca `len(REF)==1 and 1<=len(ALT)<=2` icin calisir; yani
    1 bp'lik insersiyon (REF=T ALT=TG) yanlislikla **SNP** sayilir. Onu da acikca INS diyoruz.
  * INS kayitlarinda **INSSEQ** yoksa insilicoSV RASTGELE dizi uretir (kaynak: ayni dosya,
    "a random sequence is used"). Gercek dizi elimizdeyken bunu birakmak veri kaybidir.

KOORDINAT SOZLESMESI (dokumandan degil, koddan turetildi + ampirik dogrulanacak)
------------------------------------------------------------------------------
insilicoSV `rec_len = rec_end.pos - rec_start.pos` kullaniyor; pysam'da `start` 0-tabanli,
`stop` INFO/END'ten geliyor. Dokumandaki DEL ornegi (POS=6, END=60, SVLEN=55) ile tutarli
tek okuma: **POS = silinen ilk baz, END = silinen son baz (1-tabanli, kapali aralik).**

  VCF delesyonu  REF=TG ALT=T  @P   ->  insilicoSV: POS=P+1, END=P+1,   SVLEN=1
  VCF delesyonu  REF=TGGG ALT=T @P  ->  insilicoSV: POS=P+1, END=P+3,   SVLEN=3
  VCF insersiyon REF=T ALT=TGGG @P  ->  insilicoSV: POS=P+1, END=P+1,   SVLEN=3, INSSEQ=GGG

Insersiyon icin dokuman "POS, novel dizinin ONUNE eklendigi bazdir" diyor; P bazindan SONRA
eklemek = P+1 bazinin ONUNE eklemek. Bu sozlesme `scripts/verify_conversion.py` ile
**dizi duzeyinde** dogrulanmadan hicbir sey uretilmeyecek.

ORTUSME
-------
insilicoSV import edilen varyantlarin ayni haplotipte ortusmesine izin vermiyor
(`simulate.py:820`). `--drop-overlaps` ayni haplotipte cakisanlari eler ve neyin neden
elendigini `--report` dosyasina yazar (sessiz kayip yok).
"""
import argparse
import gzip
import sys
from collections import Counter


def opener(p):
    return gzip.open(p, "rt") if p.endswith((".gz", ".bgz")) else open(p)


def parse_info(s):
    d = {}
    if s == ".":
        return d
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        elif kv:
            d[kv] = True
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--sample", help="cok ornekli VCF'te hangi ornek (varsayilan: ilk)")
    ap.add_argument("--contig-rename", action="append", default=[],
                    help="eski=yeni, birden fazla kez verilebilir")
    ap.add_argument("--drop-overlaps", action="store_true")
    ap.add_argument("--min-gap", type=int, default=1,
                    help="ayni haplotipte iki varyant arasinda birakilacak en az bosluk (bp)")
    ap.add_argument("--report", help="elenen/donusturulen kayitlarin dokumu")
    ap.add_argument("--keep-symbolic", action="store_true", default=True)
    ap.add_argument("--seqless-ins", choices=["drop", "random"], default="drop",
                    help="dizisi verilmemis <INS> kayitlari: dusur (varsayilan) veya insilicoSV rastgele uretsin")
    ap.add_argument("--ref", help="referans FASTA (tek kromozom). Verilirse: REF eslesmesi "
                                  "dogrulanir ve N iceren baglamlar elenir.")
    ap.add_argument("--max-n-frac", type=float, default=0.05,
                    help="insilicoSV th_proportion_N ile ayni olmali (varsayilan 0.05)")
    a = ap.parse_args()

    # --- referans (istege bagli ama siddetle onerilen) ---
    # NEDEN: insilicoSV, bir varyantin [start-1, end+1] penceresindeki N orani
    # th_proportion_N'i (varsayilan 0.05) asarsa yerlestirmeyi REDDEDER ve tum kosuyu
    # hata koduyla dusurur (simulate.py:is_placement_valid). hg38'de buyuk bosluklarin
    # DISINDA da tek tek serpistirilmis N bazlari var (ornek: chr2:31,591,787). Bu yuzden
    # once biz eliyoruz. Ayrica REF alaninin referansla gercekten eslestigini dogruluyoruz —
    # eslesmiyorsa varyant yanlis koordinat uzayindan gelmis demektir.
    refseq = None
    if a.ref:
        with open(a.ref) as fh:
            buf = []
            for line in fh:
                if line[0] == ">":
                    if buf:
                        break
                    continue
                buf.append(line.strip())
        refseq = "".join(buf).upper()

    rename = dict(x.split("=", 1) for x in a.contig_rename)
    stats = Counter()
    dropped = []
    out_recs = []          # (chrom, pos, end, hapmask, line_fields)
    contigs = {}

    gt_col = None
    with opener(a.inp) as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##contig="):
                    # ##contig=<ID=chr20,length=64444167>
                    body = line.strip()[len("##contig=<"):-1]
                    kv = dict(p.split("=", 1) for p in body.split(",") if "=" in p)
                    cid = rename.get(kv.get("ID"), kv.get("ID"))
                    if cid and "length" in kv:
                        contigs[cid] = int(kv["length"])
                continue
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
                gt_col = 9
                if a.sample and a.sample in hdr:
                    gt_col = hdr.index(a.sample)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                stats["no_genotype"] += 1
                continue
            chrom = rename.get(f[0], f[0])
            pos = int(f[1])
            ref, alt_field = f[3].upper(), f[4].upper()
            info = parse_info(f[7])
            gt = f[gt_col].split(":")[0]
            sep = "|" if "|" in gt else "/"
            parts = gt.replace("/", "|").split("|")
            if len(parts) != 2:
                stats["bad_gt"] += 1
                continue
            # HGSVC3 (assembly tabanli) '.' = o haplotipte HIZALAMA YOK, yani cagrilamaz --
            # "referans" demek degil. Iki alel de '.' ise varyant hakkinda hicbir sey bilmiyoruz:
            # eliyoruz. Tek alel '.' ve digeri ALT ise, bildigimiz haplotipte varyant GERCEKTEN
            # var; bilinmeyen tarafi referans kabul edip heterozigot yaziyoruz (kayip yerine
            # muhafazakar atama). Tek alel '.' ve digeri 0 ise tasidigi bilgi yok, eliyoruz.
            if parts[0] == "." and parts[1] == ".":
                stats["gt_both_uncallable"] += 1
                continue
            if parts[0] == "." or parts[1] == ".":
                other = parts[1] if parts[0] == "." else parts[0]
                if other == "0":
                    stats["gt_uncallable_noalt"] += 1
                    continue
                parts = ["0", other] if parts[0] == "." else [other, "0"]
                stats["gt_half_callable_as_het"] += 1
            hapA, hapB = parts[0] != "0", parts[1] != "0"
            if not (hapA or hapB):
                stats["hom_ref"] += 1
                continue
            # cok alelli: ornegin tasidigi ilk alt aleli sec
            alts = alt_field.split(",")
            idx = int(parts[0]) if parts[0] != "0" else int(parts[1])
            if idx - 1 >= len(alts):
                stats["bad_allele_index"] += 1
                continue
            alt = alts[idx - 1]
            # ayni bireyde iki FARKLI alt alel (1|2) -> desteklenmiyor, ele
            if parts[0] != "0" and parts[1] != "0" and parts[0] != parts[1]:
                stats["multiallelic_het"] += 1
                dropped.append((chrom, pos, "multiallelic_het", gt, ref[:20], alt_field[:20]))
                continue
            new_gt = ("1" if hapA else "0") + "|" + ("1" if hapB else "0")

            newinfo = []
            if alt.startswith("<"):
                # sembolik SV: SVTYPE/END/SVLEN INFO'dan gelmeli
                svt = info.get("SVTYPE") or alt.strip("<>").split(":")[0]
                svt = {"CNV": "DUP", "LINE1": "INS", "ALU": "INS", "SVA": "INS"}.get(svt, svt)
                if svt not in ("DEL", "DUP", "INV", "INS"):
                    stats[f"symbolic_unsupported_{svt}"] += 1
                    dropped.append((chrom, pos, f"symbolic_unsupported:{svt}", gt, "", alt))
                    continue
                if "END" in info:
                    end = int(info["END"])
                elif "SVLEN" in info:
                    end = pos + abs(int(str(info["SVLEN"]).split(",")[0]))
                else:
                    stats["symbolic_no_end"] += 1
                    dropped.append((chrom, pos, "symbolic_no_end", gt, "", alt))
                    continue
                if svt == "INS":
                    svlen = abs(int(str(info.get("SVLEN", end - pos)).split(",")[0])) or 1
                    start_out, end_out = pos + 1, pos + 1
                    if "INSSEQ" in info:
                        newinfo = [f"SVTYPE=INS", f"END={end_out}", f"SVLEN={svlen}",
                                   f"INSSEQ={info['INSSEQ']}"]
                        stats["symbolic_INS_with_seq"] += 1
                    elif a.seqless_ins == "drop":
                        # DIZISI OLMAYAN INSERSIYON. Kisa-okuma callsetlerinde (1kGP) <INS>
                        # kayitlari yalnizca uzunluk tasir. insilicoSV bu durumda RASTGELE dizi
                        # uretir; rastgele bir 300-mer ise gercek bir Alu/L1 insersiyonuna hic
                        # benzemez (tekrarli degil, poly-A kuyrugu yok) ve HiFi okumalarinda
                        # sahte gorunur. Sessizce yanlis veri uretmektense dusuruyoruz; yerine
                        # dizi-cozunurlu (uzun-okuma) bir callset konacak.
                        stats["seqless_ins_dropped"] += 1
                        dropped.append((chrom, pos, f"seqless_INS_len{svlen}", gt, "", alt))
                        continue
                    else:  # "random": insilicoSV uretsin, ama DOGRU uzunlukta
                        end_out = pos + 1 + svlen        # rec_len = END-START = svlen
                        newinfo = [f"SVTYPE=INS", f"END={end_out}", f"SVLEN={svlen}"]
                        stats["seqless_ins_random"] += 1
                else:
                    start_out, end_out = pos + 1, end
                    if end_out <= start_out:
                        stats["symbolic_zero_len"] += 1
                        continue
                    newinfo = [f"SVTYPE={svt}", f"END={end_out}", f"SVLEN={end_out-start_out}"]
                    stats[f"symbolic_{svt}"] += 1
                span = (start_out, end_out)
                ref_out, alt_out = "N", f"<{svt}>"
            elif len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
                stats["snv"] += 1
                start_out, end_out = pos, pos
                ref_out, alt_out = ref, alt
                newinfo = []                       # SNP otomatik algilanir
                span = (pos, pos)
            elif set(ref) <= set("ACGTN") and set(alt) <= set("ACGTN"):
                if len(ref) == 1 and len(alt) > 1 and alt[0] == ref[0]:
                    ins = alt[1:]
                    # insilicoSV INSSEQ'in TAMAMEN TCGA olmasini sart kosuyor
                    # (variant_set.py:1092). HGSVC3'te, icinde assembly boslugu (N dizisi) tasiyan
                    # birkac buyuk insersiyon var (genom capinda 3 tane; chr14'te 1, chrX'te 2) --
                    # ve bu 3 kayit iki kromozomun tamamini dusurdu. Dizisi bilinmeyen bir
                    # insersiyonu simule etmenin anlami da yok: eliyoruz ve raporluyoruz.
                    if any(b not in "ACGT" for b in ins):
                        stats["ins_seq_has_N_dropped"] += 1
                        dropped.append((chrom, pos, f"INS_seq_nonACGT_len{len(ins)}", gt,
                                        ref[:12], alt[:12]))
                        continue
                    stats["ins"] += 1
                    start_out = end_out = pos + 1
                    newinfo = [f"SVTYPE=INS", f"END={end_out}", f"SVLEN={len(ins)}", f"INSSEQ={ins}"]
                    ref_out, alt_out = "N", "<INS>"
                    span = (pos + 1, pos + 1)
                elif len(alt) == 1 and len(ref) > 1 and ref[0] == alt[0]:
                    dl = len(ref) - 1
                    stats["del"] += 1
                    start_out, end_out = pos + 1, pos + dl
                    newinfo = [f"SVTYPE=DEL", f"END={end_out}", f"SVLEN={dl}"]
                    ref_out, alt_out = "N", "<DEL>"
                    span = (start_out, end_out)
                else:
                    stats["complex_indel_dropped"] += 1
                    dropped.append((chrom, pos, "complex_indel", gt, ref[:20], alt[:20]))
                    continue
            else:
                stats["other_dropped"] += 1
                dropped.append((chrom, pos, "other", gt, ref[:20], alt[:20]))
                continue

            # --- referansa karsi dogrulama (insilicoSV'nin reddedecegi seyleri onceden ele) ---
            if refseq is not None:
                s0 = max(0, span[0] - 2)          # [start-1, end+1] penceresi, 0-tabanli
                e0 = min(len(refseq), span[1] + 1)
                win = refseq[s0:e0]
                if win:
                    nfrac = (win.count("N") + win.count("n")) / len(win)
                    if nfrac > a.max_n_frac:
                        stats["dropped_N_context"] += 1
                        dropped.append((chrom, pos, f"N_context_{nfrac:.2f}", gt, ref[:20], alt[:20]))
                        continue
                # SNV/indel icin REF gercekten referansla eslesiyor mu?
                if not alt.startswith("<"):
                    exp = refseq[pos - 1: pos - 1 + len(ref)]
                    if exp != ref:
                        stats["dropped_ref_mismatch"] += 1
                        dropped.append((chrom, pos, f"ref_mismatch_exp_{exp[:12]}", gt,
                                        ref[:12], alt[:12]))
                        continue

            out_recs.append((chrom, start_out, span, (hapA, hapB),
                             [chrom, str(start_out), ".", ref_out, alt_out, "100", "PASS",
                              ";".join(newinfo) if newinfo else ".", "GT", new_gt]))

    # ---- ayni haplotipte ortusme temizligi ----
    kept = []
    if a.drop_overlaps:
        by_chrom = {}
        for r in out_recs:
            by_chrom.setdefault(r[0], []).append(r)
        for chrom, recs in by_chrom.items():
            recs.sort(key=lambda r: (r[2][0], r[2][1]))
            lastA = lastB = -10**12
            for r in recs:
                s, e = r[2]
                hA, hB = r[3]
                conflict = (hA and s - lastA < a.min_gap) or (hB and s - lastB < a.min_gap)
                if conflict:
                    stats["overlap_dropped"] += 1
                    dropped.append((chrom, s, "same_hap_overlap", f"{hA}{hB}", "", ""))
                    continue
                if hA:
                    lastA = max(lastA, e)
                if hB:
                    lastB = max(lastB, e)
                kept.append(r)
    else:
        kept = out_recs
    kept.sort(key=lambda r: (r[0], r[2][0]))

    with open(a.out, "w") as fo:
        fo.write("##fileformat=VCFv4.2\n")
        for c, L in contigs.items():
            fo.write(f"##contig=<ID={c},length={L}>\n")
        for tag, num, typ in (("SVTYPE", "1", "String"), ("END", "1", "Integer"),
                              ("SVLEN", "1", "Integer"), ("INSSEQ", "1", "String")):
            fo.write(f'##INFO=<ID={tag},Number={num},Type={typ},Description="{tag}">\n')
        fo.write('##ALT=<ID=DEL,Description="Deletion">\n##ALT=<ID=INS,Description="Insertion">\n')
        fo.write('##ALT=<ID=DUP,Description="Duplication">\n##ALT=<ID=INV,Description="Inversion">\n')
        fo.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fo.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for i, r in enumerate(kept):
            fields = list(r[4])
            fields[2] = f"v{i}"
            fo.write("\t".join(fields) + "\n")

    if a.report:
        with open(a.report, "w") as fr:
            fr.write("chrom\tpos\treason\tgt\tref\talt\n")
            for d in dropped:
                fr.write("\t".join(str(x) for x in d) + "\n")

    print(f"girdi kayitlari islendi -> yazilan: {len(kept)}", file=sys.stderr)
    for k, v in sorted(stats.items(), key=lambda kv: -kv[1]):
        print(f"  {k:32s} {v:>10,}", file=sys.stderr)


if __name__ == "__main__":
    main()
