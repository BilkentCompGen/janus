#!/usr/bin/env python3
"""
Analyze insertion calling accuracy by comparing predicted VCF to ground truth.
Now includes analysis of False Positive SNPs to check for artifacts near insertions.
"""

import argparse
import sys
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class Variant:
    chrom: str
    pos: int
    type: str  # INS, DEL, SNP
    length: int  # Length of event (0 or 1 for SNP)
    id: str = ""

    @property
    def start(self):
        return self.pos
    
    @property
    def end(self):
        return self.pos


def parse_vcf(path: str, chrom_filter: str = None) -> list[Variant]:
    """Parse VCF file and extract variants (INS and SNP)."""
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
            info = fields[7]
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
                else:
                    info_dict[item] = True
            
            svtype = info_dict.get('SVTYPE', '')
            
            # Helper to detect SNP if SVTYPE missing
            ref = fields[3]
            alt = fields[4]
            
            if not svtype:
                if len(ref) == 1 and len(alt) == 1 and ref != alt:
                    svtype = 'SNP'
                elif len(alt) > len(ref):
                    svtype = 'INS'
                elif len(ref) > len(alt):
                    svtype = 'DEL'
            
            if svtype not in ['INS', 'SNP']: 
                continue # We focus on INS and associated artifacts (SNPs)
            
            length = 0
            if svtype == 'INS':
                svlen = int(info_dict.get('SVLEN', 0))
                length = abs(svlen) if svlen != 0 else abs(len(alt) - len(ref))
            elif svtype == 'SNP':
                length = 1

            variants.append(Variant(chrom, pos, svtype, length))
    
    return variants


def find_matches(truth: list[Variant], pred: list[Variant], tolerance: int = 50):
    """
    Match truth insertions to predictions. 
    Returns:
      matches: list of (truth, pred, dist)
      missed: list of truth
      extra: list of pred (False Positives)
    """
    # Separate by type
    truth_ins = [v for v in truth if v.type == 'INS']
    pred_ins = [v for v in pred if v.type == 'INS']
    pred_snps = [v for v in pred if v.type == 'SNP']
    
    # Match INS to INS
    matches = []
    matched_pred_ids = set()
    
    missed = []
    
    for t in truth_ins:
        best_p = None
        best_dist = float('inf')
        best_idx = -1
        
        for i, p in enumerate(pred_ins):
            if i in matched_pred_ids: continue # Greedy matching? Or allow multiple? 
            # Simple greedy: find closest unused
            if p.chrom != t.chrom: continue
            
            dist = abs(p.pos - t.pos)
            if dist <= tolerance and dist < best_dist:
                best_dist = dist
                best_p = p
                best_idx = i
        
        if best_p:
            matches.append((t, best_p, best_dist))
            matched_pred_ids.add(best_idx)
        else:
            missed.append(t)

    # Extra INS
    extra_ins = [p for i, p in enumerate(pred_ins) if i not in matched_pred_ids]
    
    # Extra SNPs (assuming truth has no SNPs as per current simulation context)
    extra_snps = pred_snps
    
    return matches, missed, extra_ins, extra_snps


def analyze_artifacts(extra_snps: list[Variant], truth_ins: list[Variant], tolerance: int = 100):
    """
    Check if extra SNPs are close to truth insertions (artifacts).
    """
    artifacts = []
    independent = []
    
    for snp in extra_snps:
        nearest_ins = None
        min_dist = float('inf')
        
        for ins in truth_ins:
            if ins.chrom != snp.chrom: continue
            dist = abs(ins.pos - snp.pos)
            if dist < min_dist:
                min_dist = dist
                nearest_ins = ins
        
        if min_dist <= tolerance:
            artifacts.append((snp, nearest_ins, min_dist))
        else:
            independent.append((snp, min_dist))
            
    return artifacts, independent


def main():
    parser = argparse.ArgumentParser(description='Analyze insertion calling and artifact SNPs')
    parser.add_argument('--truth', '-t', required=True, help='Ground truth VCF')
    parser.add_argument('--pred', '-p', required=True, help='Predicted VCF')
    parser.add_argument('--chrom', '-c', default=None, help='Filter chromosome')
    parser.add_argument('--tolerance', '-T', type=int, default=50, help='Match tolerance for INS')
    args = parser.parse_args()
    
    truth = parse_vcf(args.truth, args.chrom)
    pred = parse_vcf(args.pred, args.chrom)
    
    truth_ins_count = sum(1 for v in truth if v.type == 'INS')
    pred_ins_count = sum(1 for v in pred if v.type == 'INS')
    pred_snp_count = sum(1 for v in pred if v.type == 'SNP')
    
    print(f"Loaded {len(truth)} truth variants ({truth_ins_count} INS)")
    print(f"Loaded {len(pred)} predicted variants ({pred_ins_count} INS, {pred_snp_count} SNP)")
    print()
    
    matches, missed, extra_ins, extra_snps = find_matches(truth, pred, args.tolerance)
    
    print("="*60)
    print("INSERTION ANALYSIS")
    print("="*60)
    print(f"Recall: {len(matches)} / {truth_ins_count} ({100*len(matches)/truth_ins_count if truth_ins_count else 0:.1f}%)")
    print(f"Precision (INS): {len(matches)} / {pred_ins_count} ({100*len(matches)/pred_ins_count if pred_ins_count else 0:.1f}%)")
    if missed:
        print(f"Missed: {len(missed)}")
        for m in missed[:10]: print(f"  {m.chrom}:{m.pos} len={m.length}")
    else:
        print("Missed: 0")
        
    print()
    
    print("="*60)
    print("FALSE POSITIVE SNP ANALYSIS (Artifact Check)")
    print("="*60)
    print(f"Total Extra SNPs: {len(extra_snps)}")
    
    if extra_snps:
        artifacts, independent = analyze_artifacts(extra_snps, [v for v in truth if v.type == 'INS'], tolerance=100)
        
        print(f"Artifacts (within 100bp of a Truth INS): {len(artifacts)} ({100*len(artifacts)/len(extra_snps):.1f}%)")
        print(f"Independent (likely errors): {len(independent)}")
        
        if artifacts:
            print("\nExample Artifacts:")
            for snp, ins, dist in artifacts[:10]:
                print(f"  SNP {snp.chrom}:{snp.pos} -> Near INS {ins.chrom}:{ins.pos} (dist {dist}bp)")
            if len(artifacts) > 10: print(f"  ... {len(artifacts)-10} more")
            
        if independent:
            print("\nExample Independent SNPs:")
            for snp, dist in independent[:10]:
                print(f"  SNP {snp.chrom}:{snp.pos} -> Nearest INS dist {dist}bp")


    # Cross-reference analysis: Do artifacts occur near Perfect Matches?
    print("\n" + "="*60)
    print("ARTIFACT CORRELATION WITH MATCH EXACTNESS")
    print("="*60)

    # matches is list of (truth, pred, dist)
    exact_matches = [m for m in matches if m[2] == 0]
    inexact_matches = [m for m in matches if m[2] > 0]
    
    def count_associated_artifacts(match_list, snps, tolerance=50):
        count_with_artifacts = 0
        total_assoc = 0
        for t, p, d in match_list:
            # any snp near p
            nearby = [s for s in snps if s.chrom == p.chrom and abs(s.pos - p.pos) <= tolerance]
            if nearby:
                count_with_artifacts += 1
                total_assoc += len(nearby)
        return count_with_artifacts, total_assoc

    em_cnt, em_art = count_associated_artifacts(exact_matches, extra_snps)
    im_cnt, im_art = count_associated_artifacts(inexact_matches, extra_snps)
    
    print(f"Exact Matches (dist=0): {len(exact_matches)}")
    if exact_matches:
        print(f"  -> with nearby SNPs: {em_cnt} ({100*em_cnt/len(exact_matches):.1f}%)")
        print(f"  -> total artifacts: {em_art}")
        
    print(f"Inexact Matches (dist>0): {len(inexact_matches)}")
    if inexact_matches:
        print(f"  -> with nearby SNPs: {im_cnt} ({100*im_cnt/len(inexact_matches):.1f}%)")
        print(f"  -> total artifacts: {im_art}")
        
    print("\n" + "="*60)
    print("INEXACT MATCH SIZE BREAKDOWN")
    print("="*60)
    
    # Check sizes of inexact matches
    small_inexact = [m for m in inexact_matches if m[0].length <= 1000]
    large_inexact = [m for m in inexact_matches if m[0].length > 1000]
    
    small_exact = [m for m in exact_matches if m[0].length <= 1000]
    large_exact = [m for m in exact_matches if m[0].length > 1000]
    
    print(f"Small/Medium (<= 1000bp): {len(small_inexact)} Inexact vs {len(small_exact)} Exact")
    print(f"  -> Inexact Rate: {100*len(small_inexact)/(len(small_inexact)+len(small_exact)):.1f}%")
    
    print(f"Large (> 1000bp): {len(large_inexact)} Inexact vs {len(large_exact)} Exact")
    print(f"  -> Inexact Rate: {100*len(large_inexact)/(len(large_inexact)+len(large_exact)):.1f}%")
    
    # Very Large > 30k
    huge_inexact = [m for m in inexact_matches if m[0].length > 30000]
    huge_exact = [m for m in exact_matches if m[0].length > 30000]
    print(f"Very Large (> 30kbp): {len(huge_inexact)} Inexact vs {len(huge_exact)} Exact")
    if huge_inexact:
        print("  -> NOTE: These used the boundary scanning logic (if thresholds matched).")


if __name__ == '__main__':
    main()
