#!/usr/bin/env python3
"""Compare predicted SV calls with ground truth VCF"""

import sys
from collections import defaultdict

def parse_gt_vcf(vcf_path):
    """Parse ground truth VCF deletions"""
    deletions = []
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, id_val, ref, alt, qual, filter_val, info = fields[:8]
            
            if alt != '<DEL>':
                continue
            
            # Extract END from INFO
            end = None
            for info_field in info.split(';'):
                if info_field.startswith('END='):
                    end = int(info_field.split('=')[1])
                    break
            
            if end is None:
                continue
            
            start_pos = int(pos)
            if end <= start_pos:
                continue
            
            # Remove 'chr' prefix if present
            if chrom.startswith('chr'):
                chrom = chrom[3:]
            
            deletions.append({
                'chrom': chrom,
                'start': start_pos,
                'end': end
            })
    
    return deletions

def parse_predicted_tsv(tsv_path):
    """Parse predicted SV calls"""
    predictions = []
    with open(tsv_path, 'r') as f:
        header = f.readline().strip()
        if not header:
            return predictions
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            
            type_val, chrom, start_str, end_str = fields[0], fields[1], fields[2], fields[3]
            
            if type_val != 'DEL':
                continue
            
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue
            
            predictions.append({
                'chrom': chrom,
                'start': start,
                'end': end
            })
    
    return predictions

def overlap(start1, end1, start2, end2, min_overlap_ratio=0.5):
    """Check if two intervals overlap with minimum ratio"""
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    
    if overlap_end <= overlap_start:
        return False
    
    overlap_len = overlap_end - overlap_start
    len1 = end1 - start1
    len2 = end2 - start2
    
    # Check if overlap is at least min_overlap_ratio of both intervals
    ratio1 = overlap_len / len1 if len1 > 0 else 0
    ratio2 = overlap_len / len2 if len2 > 0 else 0
    
    return ratio1 >= min_overlap_ratio or ratio2 >= min_overlap_ratio

def compare(gt_deletions, predictions):
    """Compare GT deletions with predictions"""
    # Group by chromosome
    gt_by_chrom = defaultdict(list)
    pred_by_chrom = defaultdict(list)
    
    for gt in gt_deletions:
        gt_by_chrom[gt['chrom']].append(gt)
    
    for pred in predictions:
        pred_by_chrom[pred['chrom']].append(pred)
    
    # Match predictions to GT
    tp = 0  # True positives
    fp = 0  # False positives
    fn = 0  # False negatives
    
    matched_gt = set()
    matched_pred = set()
    
    # For each prediction, find matching GT
    for pred_idx, pred in enumerate(predictions):
        chrom = pred['chrom']
        if chrom not in gt_by_chrom:
            fp += 1
            continue
        
        found_match = False
        for gt_idx, gt in enumerate(gt_by_chrom[chrom]):
            if (chrom, gt_idx) in matched_gt:
                continue
            
            if overlap(pred['start'], pred['end'], gt['start'], gt['end']):
                tp += 1
                matched_gt.add((chrom, gt_idx))
                matched_pred.add((chrom, pred_idx))
                found_match = True
                break
        
        if not found_match:
            fp += 1
    
    # Count unmatched GT as false negatives
    for chrom, gt_list in gt_by_chrom.items():
        for gt_idx in range(len(gt_list)):
            if (chrom, gt_idx) not in matched_gt:
                fn += 1
    
    # Calculate metrics
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return {
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'precision': precision,
        'recall': recall,
        'f1': f1
    }

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <predicted_tsv> <gt_vcf>")
        sys.exit(1)
    
    pred_path = sys.argv[1]
    gt_path = sys.argv[2]
    
    print(f"Loading predictions from {pred_path}...")
    predictions = parse_predicted_tsv(pred_path)
    print(f"Found {len(predictions)} predicted deletions")
    
    print(f"Loading GT from {gt_path}...")
    gt_deletions = parse_gt_vcf(gt_path)
    print(f"Found {len(gt_deletions)} GT deletions")
    
    print("\nComparing...")
    results = compare(gt_deletions, predictions)
    
    print("\nResults:")
    print(f"  True Positives (TP):  {results['tp']}")
    print(f"  False Positives (FP): {results['fp']}")
    print(f"  False Negatives (FN): {results['fn']}")
    print(f"  Precision:            {results['precision']:.4f}")
    print(f"  Recall:               {results['recall']:.4f}")
    print(f"  F1-Score:             {results['f1']:.4f}")

