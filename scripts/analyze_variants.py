#!/usr/bin/env python3
"""
Analyze variant calling accuracy (DEL, INS, SNP) by comparing predicted VCF to ground truth.
"""

import argparse
import sys
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class Variant:
    chrom: str
    start: int
    end: int
    svtype: str
    ref: str = ""
    alt: str = ""
    svlen: int = 0
    
    @property
    def length(self):
        if self.svtype == 'SNP':
            return 1
        return self.end - self.start if self.svlen == 0 else abs(self.svlen)


def parse_vcf(path: str, chrom_filter: str = None) -> list[Variant]:
    """Parse VCF file and extract variants."""
    variants = []
    
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            if chrom_filter and chrom != chrom_filter:
                continue
            
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
            
            svtype = info_dict.get('SVTYPE', 'SNP')
            # If alt is standard DNA and length is 1, it might be a SNP
            if svtype == 'SNP' or (len(ref) == 1 and len(alt) == 1 and alt in 'ACGTN'):
                variants.append(Variant(chrom, pos, pos, 'SNP', ref, alt, 0))
                continue

            end = int(info_dict.get('END', pos))
            svlen = int(info_dict.get('SVLEN', 0))
            
            variants.append(Variant(chrom, pos, end, svtype, ref, alt, svlen))
    
    return variants


def overlap_ratio(a: Variant, b: Variant) -> float:
    """Calculate overlap ratio between two variants."""
    if a.chrom != b.chrom or a.svtype != b.svtype:
        return 0.0
    
    if a.svtype == 'SNP':
        return 1.0 if a.start == b.start and a.alt == b.alt else 0.0

    overlap_start = max(a.start, b.start)
    overlap_end = min(a.end, b.end)
    overlap_len = max(0, overlap_end - overlap_start)
    
    if overlap_len == 0:
        return 0.0
    
    len_a = a.length
    len_b = b.length
    
    # Reciprocal overlap for SVs
    ratio_a = overlap_len / len_a if len_a > 0 else 0
    ratio_b = overlap_len / len_b if len_b > 0 else 0
    
    return min(ratio_a, ratio_b)


def evaluate(truth: list[Variant], pred: list[Variant], sv_type: str, min_overlap: float = 0.5):
    """Evaluate predictions against ground truth for a specific type."""
    t_list = [v for v in truth if v.svtype == sv_type]
    p_list = [v for v in pred if v.svtype == sv_type]
    
    matched_truth = set()
    matched_pred = set()
    
    for i, p in enumerate(p_list):
        best_match = -1
        best_score = 0
        
        for j, t in enumerate(t_list):
            if j in matched_truth:
                continue
            
            score = overlap_ratio(p, t)
            if score > 0 and (sv_type == 'SNP' or score >= min_overlap):
                if score > best_score:
                    best_score = score
                    best_match = j
        
        if best_match >= 0:
            matched_truth.add(best_match)
            matched_pred.add(i)
    
    tp = len(matched_pred)
    fp = len(p_list) - tp
    fn = len(t_list) - tp
    
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    return tp, fp, fn, precision, recall, f1


def main():
    parser = argparse.ArgumentParser(description='Analyze variant calling accuracy (DEL, INS, SNP)')
    parser.add_argument('--truth', '-t', required=True, help='Ground truth VCF file')
    parser.add_argument('--pred', '-p', required=True, help='Predicted VCF file')
    parser.add_argument('--chrom', '-c', default=None, help='Filter by chromosome')
    args = parser.parse_args()
    
    truth = parse_vcf(args.truth, args.chrom)
    pred = parse_vcf(args.pred, args.chrom)
    
    print(f"--- Statistics ---")
    print(f"Truth: {len(truth)} variants")
    print(f"Pred:  {len(pred)} variants")
    print()
    
    for vtype in ['DEL', 'INS', 'SNP']:
        tp, fp, fn, prec, rec, f1 = evaluate(truth, pred, vtype)
        print(f"=== {vtype} ===")
        print(f"  TP: {tp}, FP: {fp}, FN: {fn}")
        print(f"  Precision: {prec:.4f}")
        print(f"  Recall:    {rec:.4f}")
        print(f"  F1 Score:  {f1:.4f}")
        print()


if __name__ == '__main__':
    main()

