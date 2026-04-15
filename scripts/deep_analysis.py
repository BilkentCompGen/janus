#!/usr/bin/env python3
"""Cross-reference SVP output calls with ground truth at position level."""
import sys

def parse_vcf(filename):
    variants = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split("\t")
            if len(fields) < 8: continue
            chrom = fields[0].replace("chr", "")
            pos = int(fields[1])
            info = {}
            for entry in fields[7].split(";"):
                if "=" in entry:
                    k, v = entry.split("=", 1)
                    info[k] = v
            svtype = info.get("SVTYPE", "SNP")
            svlen = abs(int(info["SVLEN"])) if "SVLEN" in info else 0
            variants.append({"chrom": chrom, "pos": pos, "svtype": svtype, "svlen": svlen, "line": line.strip()})
    return variants

def find_nearest_truth(call, truth_variants, tolerance=200):
    """Find nearest ground truth variant within tolerance."""
    best = None
    best_dist = tolerance + 1
    for tv in truth_variants:
        if call["chrom"] != tv["chrom"]: continue
        dist = abs(call["pos"] - tv["pos"])
        if dist < best_dist:
            best_dist = dist
            best = tv
    if best_dist <= tolerance:
        return best, best_dist
    return None, None

def prf(tp, fp, fn):
    """Compute precision, recall, F1."""
    p = tp / (tp + fp) if (tp + fp) > 0 else 0
    r = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0
    return p, r, f1

def analyze(run_name, output_vcf, truth_vcf):
    calls = parse_vcf(output_vcf)
    truth = parse_vcf(truth_vcf)

    print(f"\n{'='*80}")
    print(f"  DETAILED ANALYSIS: {run_name}")
    print(f"  Output: {output_vcf} ({len(calls)} calls)")
    print(f"  Truth:  {truth_vcf} ({len(truth)} variants)")
    print(f"{'='*80}")

    # 1. For each output call, find nearest truth variant
    matched_truth = set()
    correct_type = []
    wrong_type = []
    no_match = []

    for call in calls:
        tv, dist = find_nearest_truth(call, truth)
        if tv is None:
            no_match.append(call)
        else:
            truth_idx = truth.index(tv)
            if truth_idx in matched_truth:
                no_match.append(call)
                continue
            # Size check for SV calls (>=50bp): require pctsize >= 0.5
            if call["svlen"] >= 50 and tv["svlen"] >= 50:
                size_ratio = min(call["svlen"], tv["svlen"]) / max(call["svlen"], tv["svlen"])
                if size_ratio < 0.5:
                    no_match.append(call)
                    continue
            matched_truth.add(truth_idx)
            if call["svtype"] == tv["svtype"]:
                correct_type.append((call, tv, dist))
            else:
                wrong_type.append((call, tv, dist))

    # 2. Find missed truth variants
    missed = [(i, tv) for i, tv in enumerate(truth) if i not in matched_truth]

    # Report wrong-type matches
    if wrong_type:
        print(f"\n--- WRONG TYPE ({len(wrong_type)} calls detected real variants but wrong SVTYPE) ---")
        type_pairs = {}
        for call, tv, dist in wrong_type:
            key = f"called {call['svtype']}, truth {tv['svtype']}"
            if key not in type_pairs:
                type_pairs[key] = []
            type_pairs[key].append((call, tv, dist))
        for pair_name, items in sorted(type_pairs.items()):
            print(f"\n  {pair_name}: {len(items)} calls")
            for call, tv, dist in items[:5]:
                print(f"    pos={call['pos']} len={call['svlen']} | truth pos={tv['pos']} len={tv['svlen']} (dist={dist}bp)")
            if len(items) > 5:
                print(f"    ... and {len(items)-5} more")

    # Report ghost calls (no nearby truth)
    if no_match:
        print(f"\n--- GHOST FPs ({len(no_match)} calls with no nearby truth variant) ---")
        by_type = {}
        for call in no_match:
            by_type.setdefault(call["svtype"], []).append(call)
        for svtype, items in sorted(by_type.items()):
            print(f"  {svtype}: {len(items)}")
            for call in items[:3]:
                print(f"    pos={call['pos']} len={call['svlen']}")
            if len(items) > 3:
                print(f"    ... and {len(items)-3} more")

    # Report missed variants
    if missed:
        print(f"\n--- MISSED FNs ({len(missed)} truth variants not detected) ---")
        by_type = {}
        for i, tv in missed:
            by_type.setdefault(tv["svtype"], []).append(tv)
        for svtype, items in sorted(by_type.items()):
            sizes = [v["svlen"] for v in items]
            small = sum(1 for s in sizes if s < 50)
            medium = sum(1 for s in sizes if 50 <= s < 1000)
            large = sum(1 for s in sizes if 1000 <= s < 10000)
            xlarge = sum(1 for s in sizes if s >= 10000)
            print(f"  {svtype}: {len(items)} missed  (<50bp: {small}, 50-1kb: {medium}, 1-10kb: {large}, >10kb: {xlarge})")
            for v in items[:3]:
                print(f"    pos={v['pos']} len={v['svlen']}")
            if len(items) > 3:
                print(f"    ... and {len(items)-3} more")

    # Summary with P/R/F1
    tp = len(correct_type)
    fp = len(wrong_type) + len(no_match)
    fn = len(missed)
    p, r, f1 = prf(tp, fp, fn)

    print(f"\n--- SUMMARY (all variants) ---")
    print(f"  Correct type (TP): {tp}")
    print(f"  Wrong type:        {len(wrong_type)}")
    print(f"  Ghost FPs:         {len(no_match)}")
    print(f"  FP total:          {fp}")
    print(f"  Missed (FN):       {fn}")
    print(f"  Precision: {p:.4f}  Recall: {r:.4f}  F1: {f1:.4f}")

    # SV-only metrics (>=50bp)
    sv_calls = [c for c in calls if c["svlen"] >= 50]
    sv_truth = [t for t in truth if t["svlen"] >= 50]

    sv_matched_truth = set()
    sv_correct = 0
    sv_wrong = 0
    sv_ghost = 0

    for call in sv_calls:
        tv, dist = find_nearest_truth(call, sv_truth)
        if tv is None:
            sv_ghost += 1
        else:
            idx = sv_truth.index(tv)
            if idx in sv_matched_truth:
                sv_ghost += 1
                continue
            sv_matched_truth.add(idx)
            if call["svtype"] == tv["svtype"]:
                sv_correct += 1
            else:
                sv_wrong += 1

    sv_missed = len(sv_truth) - len(sv_matched_truth)
    sv_tp = sv_correct
    sv_fp = sv_wrong + sv_ghost
    sv_fn = sv_missed
    sv_p, sv_r, sv_f1 = prf(sv_tp, sv_fp, sv_fn)

    print(f"\n--- SUMMARY (SVs only, >=50bp) ---")
    print(f"  GT SVs: {len(sv_truth)}  Detected SVs: {len(sv_calls)}")
    print(f"  Correct type (TP): {sv_tp}")
    print(f"  Wrong type:        {sv_wrong}")
    print(f"  Ghost FPs:         {sv_ghost}")
    print(f"  FP total:          {sv_fp}")
    print(f"  Missed (FN):       {sv_fn}")
    print(f"  Precision: {sv_p:.4f}  Recall: {sv_r:.4f}  F1: {sv_f1:.4f}")

# Run analyses
analyze("DEL-only", "outputs/vcf/bench_del.vcf", "simulation/deletion_simulation.vcf")
analyze("INS-only", "outputs/vcf/bench_ins.vcf", "simulation/insertion_simulation.vcf")
analyze("SNP-only", "outputs/vcf/bench_snp.vcf", "simulation/snp_simulation.vcf")
analyze("DUP-only", "outputs/vcf/bench_dup.vcf", "simulation/duplication_simulation.vcf")
analyze("INV-only", "outputs/vcf/bench_inv.vcf", "simulation/inversion_simulation.vcf")
analyze("Mixed", "outputs/vcf/bench_mixed.vcf", "simulation/mixed.vcf")
