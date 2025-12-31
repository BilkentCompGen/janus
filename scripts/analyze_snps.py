#!/usr/bin/env python3
"""
Analyze SNP calling accuracy by comparing predicted VCF to ground truth.
"""

import argparse
from dataclasses import dataclass

@dataclass
class Variant:
    chrom: str
    pos: int
    type: str
    ref: str
    alt: str

def parse_vcf(path: str, chrom_filter: str = None) -> list[Variant]:
    variants = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 5: continue
            
            chrom = fields[0]
            if chrom_filter and chrom != chrom_filter: continue
            
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            # Simple logic: if REF/ALT len is 1, it's a SNP
            if len(ref) == 1 and len(alt) == 1 and ref != alt:
                variants.append(Variant(chrom, pos, "SNP", ref, alt))
            elif "SVTYPE=SNP" in line: # Or if marked explicitly
                variants.append(Variant(chrom, pos, "SNP", ref, alt))
                
    return variants

def match_snps(truth, pred):
    tp = []
    fp = []
    fn = []
    
    # Store predictions in a map for fast lookup: (chrom, pos) -> Variant
    pred_map = {(v.chrom, v.pos): v for v in pred}
    
    for t in truth:
        key = (t.chrom, t.pos)
        if key in pred_map:
            p = pred_map[key]
            if p.ref == t.ref and p.alt == t.alt:
                tp.append(t)
                del pred_map[key] # Consumed
            else:
                # Position match but allele mismatch? Count as FN for strict or TP for loose?
                # Usually strict.
                fn.append(t)
                # The prediction remains as FP unless matched to something else (unlikely for SNPs)
        else:
            fn.append(t)
            
    fp = list(pred_map.values())
    return tp, fp, fn

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--truth', required=True)
    parser.add_argument('--pred', required=True)
    parser.add_argument('--chrom', default=None)
    args = parser.parse_args()
    
    truth = parse_vcf(args.truth, args.chrom)
    pred = parse_vcf(args.pred, args.chrom)
    
    tp, fp, fn = match_snps(truth, pred)
    
    print("="*40)
    print("SNP ANALYSIS RESULTS")
    print("="*40)
    print(f"Ground Truth SNPs: {len(truth)}")
    print(f"Predicted SNPs:    {len(pred)}")
    print("-" * 40)
    print(f"True Positives (Correct): {len(tp)}")
    print(f"False Positives (Extra):  {len(fp)}")
    print(f"False Negatives (Missed): {len(fn)}")
    
    recall = len(tp) / len(truth) if truth else 0
    precision = len(tp) / len(pred) if pred else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    print("-" * 40)
    print(f"Recall:    {recall:.2%}")
    print(f"Precision: {precision:.2%}")
    print(f"F1 Score:  {f1:.4f}")
    
    if fn:
        print("\nFirst 10 Missed SNPs:")
        for v in fn[:10]: print(f"  {v.chrom}:{v.pos} {v.ref}->{v.alt}")

    if fp:
        print("\nFirst 10 Extra SNPs:")
        for v in fp[:10]: print(f"  {v.chrom}:{v.pos} {v.ref}->{v.alt}")

if __name__ == '__main__':
    main()
