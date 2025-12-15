#!/usr/bin/env python3
"""
Evaluate deletion calls against ground truth.

This script compares predicted deletions (from merged TSV) with ground truth
deletions (from VCF file) and calculates precision, recall, and F1 score.
"""

import sys
import argparse
import re
from typing import List, Tuple, Set
from collections import defaultdict


def parse_vcf(vcf_file: str) -> List[Tuple[str, int, int]]:
    """
    Parse VCF file and extract deletion intervals.
    
    Returns:
        List of (chrom, start, end) tuples for deletions
    """
    deletions = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            pos = int(fields[1])
            info = fields[7]
            
            # Check if it's a deletion
            if '<DEL>' not in fields[4]:
                continue
            
            # Extract END from INFO field
            end_match = re.search(r'END=(\d+)', info)
            if end_match:
                end = int(end_match.group(1))
                deletions.append((chrom, pos, end))
            else:
                # If no END, use POS+1 (single base deletion)
                deletions.append((chrom, pos, pos + 1))
    
    return deletions


def parse_tsv(tsv_file: str) -> List[Tuple[str, int, int]]:
    """
    Parse TSV file and extract deletion intervals.
    
    Returns:
        List of (chrom, start, end) tuples for deletions
    """
    deletions = []
    
    with open(tsv_file, 'r') as f:
        header = f.readline()  # Skip header
        if not header.startswith('type'):
            raise ValueError("Expected TSV header starting with 'type'")
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            
            type_str, chrom, start_str, end_str = fields[0], fields[1], fields[2], fields[3]
            
            if type_str == "DEL":
                try:
                    start = int(start_str)
                    end = int(end_str)
                    deletions.append((chrom, start, end))
                except ValueError:
                    print(f"Warning: Invalid coordinates in line: {line}", file=sys.stderr)
                    continue
    
    return deletions


def calculate_overlap(interval1: Tuple[int, int], interval2: Tuple[int, int]) -> int:
    """Calculate overlap length between two intervals."""
    start1, end1 = interval1
    start2, end2 = interval2
    
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    return 0


def calculate_intersection_over_union(interval1: Tuple[int, int], interval2: Tuple[int, int]) -> float:
    """Calculate Intersection over Union (IoU) for two intervals."""
    start1, end1 = interval1
    start2, end2 = interval2
    
    intersection = calculate_overlap(interval1, interval2)
    
    union_start = min(start1, start2)
    union_end = max(end1, end2)
    union = union_end - union_start
    
    if union == 0:
        return 0.0
    
    return intersection / union


def match_deletions(
    predicted: List[Tuple[str, int, int]],
    ground_truth: List[Tuple[str, int, int]],
    min_overlap_ratio: float = 0.5,
    min_iou: float = 0.0
) -> Tuple[List[Tuple], List[Tuple], List[Tuple], List[Tuple]]:
    """
    Match predicted deletions with ground truth deletions.
    
    Args:
        predicted: List of (chrom, start, end) for predicted deletions
        ground_truth: List of (chrom, start, end) for ground truth deletions
        min_overlap_ratio: Minimum overlap ratio to consider a match (default: 0.5)
        min_iou: Minimum IoU to consider a match (default: 0.0, i.e., any overlap)
    
    Returns:
        (true_positives, false_positives, false_negatives, matched_pairs)
    """
    # Group by chromosome
    pred_by_chrom = defaultdict(list)
    gt_by_chrom = defaultdict(list)
    
    for chrom, start, end in predicted:
        pred_by_chrom[chrom].append((start, end))
    
    for chrom, start, end in ground_truth:
        gt_by_chrom[chrom].append((start, end))
    
    true_positives = []
    false_positives = []
    false_negatives = []
    matched_pairs = []
    
    # Track which ground truth deletions have been matched
    gt_matched = set()
    
    # For each predicted deletion, find best matching ground truth
    for pred_chrom, pred_intervals in pred_by_chrom.items():
        gt_intervals = gt_by_chrom.get(pred_chrom, [])
        
        for pred_start, pred_end in pred_intervals:
            pred_interval = (pred_start, pred_end)
            best_match = None
            best_iou = 0.0
            best_overlap_ratio = 0.0
            best_gt_idx = -1
            
            for idx, (gt_start, gt_end) in enumerate(gt_intervals):
                if idx in gt_matched:
                    continue
                
                gt_interval = (gt_start, gt_end)
                iou = calculate_intersection_over_union(pred_interval, gt_interval)
                overlap = calculate_overlap(pred_interval, gt_interval)
                
                # Calculate overlap ratio relative to both intervals
                pred_len = pred_end - pred_start
                gt_len = gt_end - gt_start
                overlap_ratio_pred = overlap / pred_len if pred_len > 0 else 0
                overlap_ratio_gt = overlap / gt_len if gt_len > 0 else 0
                overlap_ratio = min(overlap_ratio_pred, overlap_ratio_gt)
                
                if iou >= min_iou and overlap_ratio >= min_overlap_ratio:
                    if iou > best_iou:
                        best_iou = iou
                        best_overlap_ratio = overlap_ratio
                        best_match = (gt_start, gt_end)
                        best_gt_idx = idx
            
            if best_match is not None:
                true_positives.append((pred_chrom, pred_start, pred_end))
                matched_pairs.append((
                    (pred_chrom, pred_start, pred_end),
                    (pred_chrom, best_match[0], best_match[1]),
                    best_iou,
                    best_overlap_ratio
                ))
                gt_matched.add(best_gt_idx)
            else:
                false_positives.append((pred_chrom, pred_start, pred_end))
    
    # Find unmatched ground truth deletions
    for chrom, gt_intervals in gt_by_chrom.items():
        for idx, (gt_start, gt_end) in enumerate(gt_intervals):
            if idx not in gt_matched:
                false_negatives.append((chrom, gt_start, gt_end))
    
    return true_positives, false_positives, false_negatives, matched_pairs


def calculate_metrics(tp: int, fp: int, fn: int) -> Tuple[float, float, float]:
    """Calculate precision, recall, and F1 score."""
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return precision, recall, f1


def print_detailed_report(
    matched_pairs: List[Tuple],
    false_positives: List[Tuple],
    false_negatives: List[Tuple],
    output_file: str = None
):
    """Print detailed comparison report."""
    out = open(output_file, 'w') if output_file else sys.stdout
    
    out.write("=" * 80 + "\n")
    out.write("DETAILED DELETION COMPARISON REPORT\n")
    out.write("=" * 80 + "\n\n")
    
    # Matched deletions
    out.write(f"MATCHED DELETIONS ({len(matched_pairs)}):\n")
    out.write("-" * 80 + "\n")
    out.write(f"{'Predicted':<30} {'Ground Truth':<30} {'IoU':<10} {'Overlap Ratio':<15}\n")
    out.write("-" * 80 + "\n")
    
    for pred, gt, iou, overlap_ratio in sorted(matched_pairs, key=lambda x: (x[0][0], x[0][1])):
        pred_str = f"{pred[0]}:{pred[1]}-{pred[2]}"
        gt_str = f"{gt[0]}:{gt[1]}-{gt[2]}"
        out.write(f"{pred_str:<30} {gt_str:<30} {iou:<10.3f} {overlap_ratio:<15.3f}\n")
    
    out.write("\n")
    
    # False positives
    if false_positives:
        out.write(f"FALSE POSITIVES ({len(false_positives)}):\n")
        out.write("-" * 80 + "\n")
        for chrom, start, end in sorted(false_positives, key=lambda x: (x[0], x[1])):
            out.write(f"{chrom}:{start}-{end}\n")
        out.write("\n")
    
    # False negatives
    if false_negatives:
        out.write(f"FALSE NEGATIVES ({len(false_negatives)}):\n")
        out.write("-" * 80 + "\n")
        for chrom, start, end in sorted(false_negatives, key=lambda x: (x[0], x[1])):
            out.write(f"{chrom}:{start}-{end}\n")
        out.write("\n")
    
    if output_file:
        out.close()


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate deletion calls against ground truth",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --vcf simulation/sim.vcf --tsv outputs/tsv/d4_merged.tsv
  %(prog)s --vcf simulation/sim.vcf --tsv outputs/tsv/d4_merged.tsv --min-overlap 0.8 --report report.txt
        """
    )
    parser.add_argument('--vcf', required=True, help='Ground truth VCF file')
    parser.add_argument('--tsv', required=True, help='Predicted deletions TSV file (merged)')
    parser.add_argument('--min-overlap', type=float, default=0.5,
                        help='Minimum overlap ratio to consider a match (default: 0.5)')
    parser.add_argument('--min-iou', type=float, default=0.0,
                        help='Minimum IoU to consider a match (default: 0.0)')
    parser.add_argument('--report', help='Output file for detailed report (optional)')
    
    args = parser.parse_args()
    
    print("Reading ground truth from VCF...", file=sys.stderr)
    ground_truth = parse_vcf(args.vcf)
    print(f"Found {len(ground_truth)} ground truth deletions", file=sys.stderr)
    
    print("Reading predictions from TSV...", file=sys.stderr)
    predicted = parse_tsv(args.tsv)
    print(f"Found {len(predicted)} predicted deletions", file=sys.stderr)
    
    print("Matching deletions...", file=sys.stderr)
    tp, fp, fn, matched_pairs = match_deletions(
        predicted, ground_truth,
        min_overlap_ratio=args.min_overlap,
        min_iou=args.min_iou
    )
    
    precision, recall, f1 = calculate_metrics(len(tp), len(fp), len(fn))
    
    # Print summary
    print("\n" + "=" * 80)
    print("EVALUATION SUMMARY")
    print("=" * 80)
    print(f"Ground Truth Deletions: {len(ground_truth)}")
    print(f"Predicted Deletions:    {len(predicted)}")
    print(f"\nTrue Positives (TP):    {len(tp)}")
    print(f"False Positives (FP):   {len(fp)}")
    print(f"False Negatives (FN):   {len(fn)}")
    print(f"\nPrecision:              {precision:.4f} ({len(tp)}/{len(tp) + len(fp)})")
    print(f"Recall:                 {recall:.4f} ({len(tp)}/{len(tp) + len(fn)})")
    print(f"F1 Score:               {f1:.4f}")
    print("=" * 80)
    
    # Print detailed report if requested
    if args.report:
        print(f"\nWriting detailed report to {args.report}...", file=sys.stderr)
        print_detailed_report(matched_pairs, fp, fn, args.report)
        print("Done!", file=sys.stderr)
    elif len(matched_pairs) <= 50:  # Only print details if not too many
        print_detailed_report(matched_pairs, fp, fn)


if __name__ == '__main__':
    main()
