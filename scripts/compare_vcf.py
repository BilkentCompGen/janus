#!/usr/bin/env python3
"""Compare VCF files to evaluate variant detection accuracy.

Reports separate metrics for:
  - SNPs + short indels (< 50bp)
  - Structural variants (>= 50bp)
  - Overall
"""

import sys


def parse_vcf(filename):
    """Parse VCF file and return list of variant dicts."""
    variants = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            chrom = fields[0].replace("chr", "")
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            # Classify by size: max of ref/alt allele length difference
            ref_len = len(ref) if ref != "." else 0
            alt_len = len(alt) if alt != "." else 0
            size = max(ref_len, alt_len, abs(ref_len - alt_len))
            variants.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "size": size,
            })
    return variants


def classify(variant):
    """Classify variant as 'small', 'del', or 'ins'."""
    if variant["size"] < 50:
        return "small"
    ref_len = len(variant["ref"]) if variant["ref"] != "." else 0
    alt_len = len(variant["alt"]) if variant["alt"] != "." else 0
    return "del" if ref_len >= alt_len else "ins"


def compare_variants(simulated, detected, tolerance=50):
    """
    Compare simulated vs detected variants.
    Returns per-category (tp, fp, fn) counts and lists.
    """
    matched_sim = set()
    matched_det = set()

    # Match detected to simulated
    for di, dv in enumerate(detected):
        for si, sv in enumerate(simulated):
            if si in matched_sim:
                continue
            if dv["chrom"] != sv["chrom"]:
                continue
            pos_diff = abs(dv["pos"] - sv["pos"])
            if pos_diff > tolerance:
                continue
            # For entries with actual sequences, require allele match
            if dv["ref"] not in (".", "N") and sv["ref"] not in (".", "N"):
                if dv["ref"].upper() == sv["ref"].upper() and dv["alt"].upper() == sv["alt"].upper():
                    matched_sim.add(si)
                    matched_det.add(di)
                    break
            else:
                matched_sim.add(si)
                matched_det.add(di)
                break

    # Categorize results
    categories = {"small": {"tp": 0, "fp": 0, "fn": 0}, "del": {"tp": 0, "fp": 0, "fn": 0}, "ins": {"tp": 0, "fp": 0, "fn": 0}}

    for di, dv in enumerate(detected):
        cat = classify(dv)
        if di in matched_det:
            categories[cat]["tp"] += 1
        else:
            categories[cat]["fp"] += 1

    for si, sv in enumerate(simulated):
        if si not in matched_sim:
            cat = classify(sv)
            categories[cat]["fn"] += 1

    return categories


def metrics(tp, fp, fn):
    """Calculate precision, recall, and F1 score."""
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    return precision, recall, f1


def print_section(label, tp, fp, fn):
    """Print a results section."""
    precision, recall, f1 = metrics(tp, fp, fn)
    print(f"  {label}:")
    print(f"    TP: {tp:>5}   FP: {fp:>5}   FN: {fn:>5}")
    print(f"    Precision: {precision:.4f}  Recall: {recall:.4f}  F1: {f1:.4f}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_vcf.py <simulated.vcf> <detected.vcf>")
        print("       make compare <simulated.vcf> <detected.vcf>")
        sys.exit(1)

    sim_file = sys.argv[1]
    det_file = sys.argv[2]

    simulated = parse_vcf(sim_file)
    detected = parse_vcf(det_file)

    print(f"Comparing:")
    print(f"  Simulated: {sim_file} ({len(simulated)} variants)")
    print(f"  Detected:  {det_file} ({len(detected)} variants)")
    print()

    cats = compare_variants(simulated, detected, tolerance=50)

    s, d, i = cats["small"], cats["del"], cats["ins"]
    total_tp = s["tp"] + d["tp"] + i["tp"]
    total_fp = s["fp"] + d["fp"] + i["fp"]
    total_fn = s["fn"] + d["fn"] + i["fn"]

    print("Results:")
    print_section("SNP + short indel (<50bp)", s["tp"], s["fp"], s["fn"])
    print_section("Deletions (>=50bp)", d["tp"], d["fp"], d["fn"])
    print_section("Insertions (>=50bp)", i["tp"], i["fp"], i["fn"])
    print_section("Overall", total_tp, total_fp, total_fn)
    print()


if __name__ == "__main__":
    main()
