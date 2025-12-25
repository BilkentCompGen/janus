#!/usr/bin/env python3
"""
Analyze deletion calling accuracy by comparing predicted VCF to ground truth.
"""

import argparse
import sys
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class Deletion:
    chrom: str
    start: int
    end: int
    svlen: int = 0
    
    @property
    def length(self):
        return self.end - self.start if self.svlen == 0 else abs(self.svlen)


def parse_vcf(path: str, chrom_filter: str = None) -> list[Deletion]:
    """Parse VCF file and extract deletions."""
    deletions = []
    
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
            info = fields[7]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
            
            svtype = info_dict.get('SVTYPE', '')
            if svtype != 'DEL':
                continue
            
            end = int(info_dict.get('END', pos))
            svlen = int(info_dict.get('SVLEN', 0))
            
            deletions.append(Deletion(chrom, pos, end, svlen))
    
    return deletions


def intervals_overlap(a_start, a_end, b_start, b_end) -> bool:
    """Check if two intervals overlap at all."""
    return a_start < b_end and b_start < a_end


def is_contained(outer_start, outer_end, inner_start, inner_end) -> bool:
    """Check if inner is contained within outer."""
    return outer_start <= inner_start and inner_end <= outer_end


def find_best_match(truth_del: Deletion, predictions: list[Deletion]) -> tuple[Deletion | None, str]:
    """
    Find best matching prediction for a truth deletion.
    Returns (matched_pred, match_type) where match_type is:
    - 'contained': pred contains truth (pred_start <= truth_start < truth_end <= pred_end)
    - 'partial': some overlap but not contained
    - 'none': no overlap at all
    """
    best_pred = None
    best_overlap = 0
    
    for pred in predictions:
        if pred.chrom != truth_del.chrom:
            continue
        
        if not intervals_overlap(pred.start, pred.end, truth_del.start, truth_del.end):
            continue
        
        overlap_start = max(pred.start, truth_del.start)
        overlap_end = min(pred.end, truth_del.end)
        overlap_len = overlap_end - overlap_start
        
        if overlap_len > best_overlap:
            best_overlap = overlap_len
            best_pred = pred
    
    if best_pred is None:
        return None, 'none'
    
    if is_contained(best_pred.start, best_pred.end, truth_del.start, truth_del.end):
        return best_pred, 'contained'
    else:
        return best_pred, 'partial'


def get_size_range(length: int) -> str:
    """Get size range category for a deletion."""
    if length <= 50:
        return '1-50'
    elif length <= 100:
        return '51-100'
    elif length <= 500:
        return '101-500'
    elif length <= 1000:
        return '501-1000'
    elif length <= 10000:
        return '1K-10K'
    else:
        return '10K+'


def get_yaml_range(length: int) -> str:
    """Get YAML-defined size range category."""
    if length <= 100:
        return 'Small (1-100bp)'
    elif length <= 1000:
        return 'Medium (100-1000bp)'
    elif length <= 10000:
        return 'Large (1K-10K)'
    else:
        return 'Very Large (10K-100K)'


def main():
    parser = argparse.ArgumentParser(description='Analyze deletion calling accuracy')
    parser.add_argument('--truth', '-t', required=True, help='Ground truth VCF file')
    parser.add_argument('--pred', '-p', required=True, help='Predicted VCF file')
    parser.add_argument('--chrom', '-c', default=None, help='Filter by chromosome')
    args = parser.parse_args()
    
    truth = parse_vcf(args.truth, args.chrom)
    pred = parse_vcf(args.pred, args.chrom)
    
    print(f"Ground truth: {len(truth)} deletions")
    print(f"Predicted: {len(pred)} deletions")
    print()
    
    # Analyze each truth deletion
    completely_missed = []  # No overlap at all
    partial_detected = []   # Some overlap but not contained
    well_detected = []      # Contained in prediction
    
    match_info = []  # (truth, pred, match_type)
    
    for t in truth:
        best_pred, match_type = find_best_match(t, pred)
        match_info.append((t, best_pred, match_type))
        
        if match_type == 'none':
            completely_missed.append(t)
        elif match_type == 'partial':
            partial_detected.append((t, best_pred))
        else:
            well_detected.append((t, best_pred))
    
    # ============================================================
    # 1. Completely Missed Deletions
    # ============================================================
    print("=" * 60)
    print("COMPLETELY MISSED DELETIONS (no overlap at all)")
    print("=" * 60)
    if completely_missed:
        for t in sorted(completely_missed, key=lambda x: x.start):
            print(f"  {t.chrom}:{t.start}-{t.end} (len={t.length})")
    else:
        print("  None! All deletions have some overlap with predictions.")
    print()
    
    # ============================================================
    # 2. Partial Detection (overlap but not contained)
    # ============================================================
    print("=" * 60)
    print("PARTIAL DETECTION (overlap but truth not contained in pred)")
    print("=" * 60)
    if partial_detected:
        for t, p in sorted(partial_detected, key=lambda x: x[0].start):
            overlap_start = max(t.start, p.start)
            overlap_end = min(t.end, p.end)
            overlap_len = overlap_end - overlap_start
            overlap_pct = 100 * overlap_len / t.length
            print(f"  Truth: {t.chrom}:{t.start}-{t.end} (len={t.length})")
            print(f"  Pred:  {p.chrom}:{p.start}-{p.end} (len={p.length})")
            print(f"  Overlap: {overlap_len}bp ({overlap_pct:.1f}%)")
            print()
    else:
        print("  None!")
    print()
    
    # ============================================================
    # 3. Size-based Analysis
    # ============================================================
    print("=" * 60)
    print("SIZE-BASED ANALYSIS (YAML categories)")
    print("=" * 60)
    
    yaml_ranges = ['Small (1-100bp)', 'Medium (100-1000bp)', 'Large (1K-10K)', 'Very Large (10K-100K)']
    
    for rng in yaml_ranges:
        range_truth = [t for t in truth if get_yaml_range(t.length) == rng]
        range_missed = [t for t in completely_missed if get_yaml_range(t.length) == rng]
        range_partial = [(t, p) for t, p in partial_detected if get_yaml_range(t.length) == rng]
        range_well = [(t, p) for t, p in well_detected if get_yaml_range(t.length) == rng]
        
        if not range_truth:
            continue
        
        print(f"\n{rng}:")
        print(f"  Total in truth: {len(range_truth)}")
        print(f"  Completely missed: {len(range_missed)}")
        print(f"  Partial detection: {len(range_partial)}")
        print(f"  Well detected (contained): {len(range_well)}")
        
        # For well-detected, calculate average expansion
        if range_well:
            expansions = []
            for t, p in range_well:
                expansion = p.length / t.length
                expansions.append(expansion)
            avg_expansion = sum(expansions) / len(expansions)
            print(f"  Avg pred/truth size ratio: {avg_expansion:.2f}x")
        
        # TP/FP/FN for this range (50% reciprocal overlap)
        tp = 0
        for t, p in range_well + range_partial:
            overlap_start = max(t.start, p.start)
            overlap_end = min(t.end, p.end)
            overlap_len = overlap_end - overlap_start
            ratio_t = overlap_len / t.length if t.length > 0 else 0
            ratio_p = overlap_len / p.length if p.length > 0 else 0
            if min(ratio_t, ratio_p) >= 0.5:
                tp += 1
        
        fn = len(range_truth) - tp
        recall = tp / len(range_truth) if range_truth else 0
        print(f"  TP (50% overlap): {tp}, FN: {fn}, Recall: {recall:.2%}")
    
    # ============================================================
    # 4. Detailed Size Brackets
    # ============================================================
    print("\n" + "=" * 60)
    print("DETAILED SIZE BRACKETS")
    print("=" * 60)
    
    brackets = ['1-50', '51-100', '101-500', '501-1000', '1K-10K', '10K+']
    
    for bracket in brackets:
        bracket_truth = [t for t in truth if get_size_range(t.length) == bracket]
        if not bracket_truth:
            continue
        
        bracket_matched = []
        bracket_missed = 0
        pred_sizes = []
        
        for t in bracket_truth:
            best_pred, match_type = find_best_match(t, pred)
            if match_type == 'none':
                bracket_missed += 1
            else:
                bracket_matched.append((t, best_pred))
                pred_sizes.append(best_pred.length)
        
        print(f"\n{bracket}bp (n={len(bracket_truth)}):")
        print(f"  Completely missed: {bracket_missed}")
        print(f"  Some overlap: {len(bracket_matched)}")
        
        if pred_sizes:
            avg_pred = sum(pred_sizes) / len(pred_sizes)
            avg_truth = sum(t.length for t, _ in bracket_matched) / len(bracket_matched)
            print(f"  Avg truth size: {avg_truth:.0f}bp")
            print(f"  Avg pred size: {avg_pred:.0f}bp")
            print(f"  Avg expansion: {avg_pred/avg_truth:.2f}x")
    
    # ============================================================
    # 5. Summary Statistics
    # ============================================================
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total ground truth deletions: {len(truth)}")
    print(f"Completely missed (0% overlap): {len(completely_missed)} ({100*len(completely_missed)/len(truth):.1f}%)")
    print(f"Partial detection: {len(partial_detected)} ({100*len(partial_detected)/len(truth):.1f}%)")
    print(f"Well detected (contained): {len(well_detected)} ({100*len(well_detected)/len(truth):.1f}%)")


if __name__ == '__main__':
    main()
