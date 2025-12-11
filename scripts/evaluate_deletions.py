#!/usr/bin/env python3
"""
Evaluate deletion detection performance by merging detected deletions
and comparing against ground truth CNVs.
"""

import argparse
import pandas as pd
from collections import defaultdict


def merge_deletions(deletions, merge_distance=1000):
    """
    Merge nearby deletions on the same chromosome.
    
    Args:
        deletions: List of deletion dictionaries with 'chrom', 'start', 'end', 'tumor_samples'
        merge_distance: Maximum distance between deletions to merge (bp)
    
    Returns:
        List of merged deletions
    """
    if not deletions:
        return []
    
    # Group by chromosome and tumor sample
    by_chrom_sample = defaultdict(list)
    for del_region in deletions:
        chrom = del_region['chrom']
        samples = del_region['tumor_samples'].split(',') if isinstance(del_region['tumor_samples'], str) else [del_region['tumor_samples']]
        for sample in samples:
            by_chrom_sample[(chrom, sample)].append({
                'start': del_region['start'],
                'end': del_region['end'],
                'sample': sample
            })
    
    merged = []
    
    for (chrom, sample), regions in by_chrom_sample.items():
        # Sort by start position
        regions = sorted(regions, key=lambda x: x['start'])
        
        current = regions[0].copy()
        
        for next_region in regions[1:]:
            # Check if should merge
            if next_region['start'] <= current['end'] + merge_distance:
                # Merge: extend current region
                current['end'] = max(current['end'], next_region['end'])
            else:
                # Save current and start new region
                merged.append({
                    'chrom': chrom,
                    'start': current['start'],
                    'end': current['end'],
                    'tumor_sample': sample
                })
                current = next_region.copy()
        
        # Don't forget the last one
        merged.append({
            'chrom': chrom,
            'start': current['start'],
            'end': current['end'],
            'tumor_sample': sample
        })
    
    return merged


def load_detected_deletions(tsv_file):
    """Load deletions from pipeline output TSV."""
    df = pd.read_csv(tsv_file, sep='\t')
    
    # Filter for deletions only
    df = df[df['type'] == 'DEL'].copy()
    
    deletions = []
    for _, row in df.iterrows():
        deletions.append({
            'chrom': str(row['chrom']),
            'start': int(row['start']),
            'end': int(row['end']),
            'tumor_samples': str(row['tumor_samples'])
        })
    
    return deletions


def load_ground_truth(tsv_file, parental=None):
    """Load ground truth deletions (copy == -1) from simulation.
    
    Args:
        tsv_file: Path to ground truth TSV
        parental: If specified (0 or 1), filter to only that parental haplotype
    """
    # Read file, skipping the # in the header
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
    
    # Fix the header line by removing the leading #
    if lines and lines[0].startswith('#'):
        lines[0] = lines[0][1:]  # Remove the # character
    
    # Read into pandas
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
    
    # Filter for deletions (copy == -1)
    df = df[df['copy'] == -1].copy()
    
    # Filter by parental if specified
    if parental is not None:
        df = df[df['parental'] == parental].copy()
    
    deletions = []
    for _, row in df.iterrows():
        deletions.append({
            'chrom': str(row['chr']),
            'start': int(row['start']),
            'end': int(row['end']),
            'tumor_sample': str(row['carrier']),
            'parental': int(row['parental'])
        })
    
    return deletions


def calculate_overlap(region1, region2):
    """Calculate overlap between two regions."""
    if region1['chrom'] != region2['chrom']:
        return 0
    
    start = max(region1['start'], region2['start'])
    end = min(region1['end'], region2['end'])
    
    if start >= end:
        return 0
    
    return end - start


def reciprocal_overlap(region1, region2):
    """Calculate reciprocal overlap ratio between two regions."""
    overlap = calculate_overlap(region1, region2)
    if overlap == 0:
        return 0.0
    
    len1 = region1['end'] - region1['start']
    len2 = region2['end'] - region2['start']
    
    return min(overlap / len1, overlap / len2)


def evaluate_deletions(detected, ground_truth, overlap_threshold=0.5):
    """
    Evaluate deletion detection performance.
    
    Args:
        detected: List of detected deletions
        ground_truth: List of ground truth deletions
        overlap_threshold: Minimum reciprocal overlap to consider a match
    
    Returns:
        Dictionary with evaluation metrics
    """
    # Track which ground truth deletions were detected
    true_positives = 0
    false_positives = 0
    matched_gt = set()
    
    # For each detected deletion, check if it matches any ground truth
    for det in detected:
        best_match = None
        best_overlap = 0
        
        for i, gt in enumerate(ground_truth):
            # Only compare same chromosome and same sample
            if det['chrom'] != gt['chrom']:
                continue
            if str(det['tumor_sample']) != str(gt['tumor_sample']):
                continue
            
            overlap = reciprocal_overlap(det, gt)
            if overlap > best_overlap:
                best_overlap = overlap
                best_match = i
        
        if best_overlap >= overlap_threshold:
            true_positives += 1
            matched_gt.add(best_match)
        else:
            false_positives += 1
    
    false_negatives = len(ground_truth) - len(matched_gt)
    
    # Calculate metrics
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    return {
        'true_positives': true_positives,
        'false_positives': false_positives,
        'false_negatives': false_negatives,
        'detected_total': len(detected),
        'ground_truth_total': len(ground_truth),
        'precision': precision,
        'recall': recall,
        'f1_score': f1
    }


def main():
    parser = argparse.ArgumentParser(description='Evaluate deletion detection performance')
    parser.add_argument('experiment', help='Experiment name (e.g., trial2)')
    parser.add_argument('--merge-distance', type=int, default=1000, help='Maximum distance to merge deletions (bp, default: 1000)')
    parser.add_argument('--overlap-threshold', type=float, default=0.5, help='Minimum reciprocal overlap for a match (default: 0.5)')
    parser.add_argument('--chrom', help='Filter to specific chromosome (e.g., 1)')
    parser.add_argument('--parental', type=int, choices=[0, 1], default=0, help='Filter ground truth to parental haplotype (default: 0 for parental_0 only)')
    parser.add_argument('--output', '-o', help='Output file for merged deletions (optional)')
    
    args = parser.parse_args()
    
    # Construct paths
    detected_tsv = f"output/tsv/{args.experiment}.tsv"
    ground_truth_tsv = "simulation/ground_truth_cnvs.tsv"
    
    print(f"Loading detected deletions from {detected_tsv}...")
    detected = load_detected_deletions(detected_tsv)
    print(f"  Loaded {len(detected)} raw deletions")
    
    print(f"\nMerging deletions (merge distance: {args.merge_distance} bp)...")
    merged = merge_deletions(detected, merge_distance=args.merge_distance)
    print(f"  After merging: {len(merged)} deletions")
    
    # Filter by chromosome if specified
    if args.chrom:
        merged = [d for d in merged if d['chrom'] == args.chrom]
        print(f"  After filtering to chr{args.chrom}: {len(merged)} deletions")
    
    if args.output:
        # Save merged deletions
        df_merged = pd.DataFrame(merged)
        df_merged.to_csv(args.output, sep='\t', index=False)
        print(f"  Saved merged deletions to {args.output}")
    
    print(f"\nLoading ground truth from {ground_truth_tsv} (parental={args.parental})...")
    ground_truth = load_ground_truth(ground_truth_tsv, parental=args.parental)
    print(f"  Loaded {len(ground_truth)} ground truth deletions (parental={args.parental})")
    
    # Filter by chromosome if specified
    if args.chrom:
        ground_truth = [d for d in ground_truth if d['chrom'] == args.chrom]
        print(f"  After filtering to chr{args.chrom}: {len(ground_truth)} ground truth deletions")
    
    # Debug: Show sample data and sample distribution
    if merged:
        print(f"\n[DEBUG] Sample detected deletion:")
        print(f"  {merged[0]}")
        detected_samples = set(d['tumor_sample'] for d in merged)
        print(f"[DEBUG] Detected samples: {sorted(detected_samples)}")
        print(f"[DEBUG] Detected count by sample: {[(s, sum(1 for d in merged if d['tumor_sample'] == s)) for s in sorted(detected_samples)]}")
    
    if ground_truth:
        print(f"\n[DEBUG] Sample ground truth deletion:")
        print(f"  {ground_truth[0]}")
        gt_samples = set(d['tumor_sample'] for d in ground_truth)
        print(f"[DEBUG] Ground truth samples: {sorted(gt_samples)}")
        print(f"[DEBUG] Ground truth count by sample: {[(s, sum(1 for d in ground_truth if d['tumor_sample'] == s)) for s in sorted(gt_samples)]}")
    
    print(f"\nEvaluating (overlap threshold: {args.overlap_threshold})...")
    metrics = evaluate_deletions(merged, ground_truth, overlap_threshold=args.overlap_threshold)
    
    print("\n" + "="*60)
    print(f"EVALUATION RESULTS - {args.experiment}")
    print("="*60)
    print(f"Ground truth deletions: {metrics['ground_truth_total']}")
    print(f"Detected deletions:     {metrics['detected_total']}")
    print()
    print(f"True Positives:  {metrics['true_positives']}")
    print(f"False Positives: {metrics['false_positives']}")
    print(f"False Negatives: {metrics['false_negatives']}")
    print()
    print(f"Precision: {metrics['precision']:.4f}")
    print(f"Recall:    {metrics['recall']:.4f}")
    print(f"F1 Score:  {metrics['f1_score']:.4f}")
    print("="*60)


if __name__ == '__main__':
    main()
