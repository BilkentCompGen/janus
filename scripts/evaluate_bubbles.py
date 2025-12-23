#!/usr/bin/env python3
import pandas as pd
import argparse, re
from pathlib import Path

def load_gt(vcf_path, min_len):
    records = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = parts[:8]
            if alt != '<DEL>':
                continue
            m = re.search(r'END=(\d+)', info)
            if not m:
                continue
            start = int(pos)
            end = int(m.group(1))
            if end <= start or (end - start) < min_len:
                continue
            records.append((chrom.replace('chr',''), start, end))
    return pd.DataFrame(records, columns=['chrom','start','end'])

def overlap(a1,a2,b1,b2):
    return max(0, min(a2,b2) - max(a1,b1))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--bubbles', required=True)
    ap.add_argument('--min-overlap', type=float, default=0.5,
                    help='IoU threshold to count as match')
    ap.add_argument('--min-gt-length', type=int, default=50)
    args = ap.parse_args()

    gt = load_gt(args.vcf, args.min_gt_length)
    bubbles = pd.read_csv(args.bubbles, sep='\t')
    bubbles['chrom'] = bubbles['chrom'].astype(str).str.replace('^chr','',regex=True)

    gt_matched = [False]*len(gt)
    bubble_matched = [False]*len(bubbles)
    matches = []

    for i, row in bubbles.iterrows():
        chrom = str(row['chrom'])
        ps, pe = int(row['start']), int(row['end'])
        cand = gt[gt['chrom']==chrom]
        best_iou, best_idx = 0.0, None
        best_ovl = 0
        for j, g in cand.iterrows():
            ovl = overlap(ps, pe, g.start, g.end)
            if ovl == 0:
                continue
            iou = ovl / ((pe-ps) + (g.end-g.start) - ovl)
            if iou > best_iou:
                best_iou = iou
                best_idx = j
                best_ovl = ovl
        if best_iou >= args.min_overlap and best_idx is not None:
            bubble_matched[i] = True
            gt_matched[best_idx] = True
            matches.append((i, best_idx, best_iou, best_ovl))

    tp = sum(bubble_matched)
    fp = len(bubbles) - tp
    fn = len(gt) - sum(gt_matched)
    prec = tp / (tp + fp) if tp+fp else 0
    rec  = tp / (tp + fn) if tp+fn else 0
    f1   = (2*prec*rec)/(prec+rec) if prec+rec else 0

    print('Ground Truth Deletions:', len(gt))
    print('Predicted (bubbles):   ', len(bubbles))
    print('TP:', tp, 'FP:', fp, 'FN:', fn)
    print(f'Precision: {prec:.4f}  Recall: {rec:.4f}  F1: {f1:.4f}')

    matches_sorted = sorted(matches, key=lambda x: x[2], reverse=True)[:10]
    print('\nTop matches (bubble_idx, gt_idx, IoU, overlap_bp):')
    for m in matches_sorted:
        print(m)

if __name__ == '__main__':
    main()
