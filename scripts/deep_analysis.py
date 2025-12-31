#!/usr/bin/env python3
"""
Deep Analysis of Variant Calling Failures
Focuses on: clustering, merging, positional shifts, size mismatches
"""

import argparse
from dataclasses import dataclass, field
from collections import defaultdict
import re

@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    svtype: str
    svlen: int
    
    @property
    def end(self):
        if self.svtype == "DEL":
            return self.pos + abs(self.svlen)
        return self.pos  # INS/SNP are point events

def parse_vcf(path: str, chrom_filter: str = None) -> dict[str, list[Variant]]:
    """Parse VCF and categorize by SVTYPE"""
    variants = {"DEL": [], "INS": [], "SNP": []}
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 8: continue
            
            chrom = fields[0]
            if chrom_filter and chrom != chrom_filter: continue
            
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            info = fields[7] if len(fields) > 7 else ""
            
            # Determine type and length
            svtype = None
            svlen = 0
            
            # Check INFO field
            if "SVTYPE=DEL" in info:
                svtype = "DEL"
                m = re.search(r'SVLEN=(-?\d+)', info)
                svlen = abs(int(m.group(1))) if m else len(ref) - len(alt)
            elif "SVTYPE=INS" in info:
                svtype = "INS"
                m = re.search(r'SVLEN=(-?\d+)', info)
                svlen = abs(int(m.group(1))) if m else len(alt) - len(ref)
            elif "SVTYPE=SNP" in info or (len(ref) == 1 and len(alt) == 1 and ref != alt):
                svtype = "SNP"
                svlen = 1
            elif len(ref) > len(alt):
                svtype = "DEL"
                svlen = len(ref) - len(alt)
            elif len(alt) > len(ref):
                svtype = "INS"
                svlen = len(alt) - len(ref)
            else:
                continue
                
            variants[svtype].append(Variant(chrom, pos, ref, alt, svtype, svlen))
    
    # Sort by position
    for k in variants:
        variants[k].sort(key=lambda v: v.pos)
    return variants

def find_clusters(variants: list[Variant], max_dist: int = 1000) -> list[list[Variant]]:
    """Find clusters of nearby variants"""
    if not variants:
        return []
    
    clusters = []
    current = [variants[0]]
    
    for v in variants[1:]:
        if v.pos - current[-1].end <= max_dist:
            current.append(v)
        else:
            if len(current) > 1:
                clusters.append(current)
            current = [v]
    
    if len(current) > 1:
        clusters.append(current)
    
    return clusters

def match_variants(truth: list[Variant], pred: list[Variant], pos_tol: int = 50, size_tol: float = 0.1):
    """Match variants with tolerance, return detailed match info"""
    tp = []
    fn = []
    fp = []

    # Filter TRUTH to match the user's focus on SVs (>50bp)
    # This prevents artificial deflation of recall
    filtered_truth = [v for v in truth if abs(v.svlen) >= 50]
    
    # We also keep a raw count for reference but use filtered/target for main metrics if requested or implied
    # Since the user explicitly requested only >50bp, let's base the MAIN stats on that
    
    matched_pred_indices = set() # Keep track of indices in the 'pred' list
    
    for t in filtered_truth:
        best_match = None
        best_score = 0
        best_idx = -1
        
        for i, p in enumerate(pred):
            if i in matched_pred_indices:
                continue
            
            pos_diff = abs(t.pos - p.pos)
            if pos_diff > pos_tol + max(t.svlen, 100):
                continue
            
            size_ratio = min(t.svlen, p.svlen) / max(t.svlen, p.svlen) if max(t.svlen, p.svlen) > 0 else 1
            pos_score = max(0, 1 - pos_diff / (pos_tol + 1))
            score = (pos_score + size_ratio) / 2
            
            if score > best_score and size_ratio >= (1 - size_tol):
                best_score = score
                best_match = p
                best_idx = i
        
        if best_match:
            tp.append((t, best_match, {
                "pos_diff": t.pos - best_match.pos,
                "size_diff": t.svlen - best_match.svlen,
                "size_ratio": t.svlen / best_match.svlen if best_match.svlen else float('inf')
            }))
            matched_pred_indices.add(best_idx)
        else:
            fn.append(t)
    
    fp = [p for i, p in enumerate(pred) if i not in matched_pred_indices]
    return tp, fp, fn

def analyze_fn_clusters(fn: list[Variant], pred: list[Variant]):
    """Analyze if FNs are being merged into larger predictions"""
    results = []
    
    for f in fn:
        # Find any prediction that might contain this FN
        containing = []
        for p in pred:
            # Check if FN falls within prediction's span
            p_start = p.pos
            p_end = p.pos + p.svlen if p.svtype == "DEL" else p.pos
            f_start = f.pos
            f_end = f.pos + f.svlen if f.svtype == "DEL" else f.pos
            
            # Overlaps or is contained?
            if p_start <= f_end + 100 and p_end >= f_start - 100:
                containing.append(p)
        
        if containing:
            results.append({
                "fn": f,
                "overlapping_preds": containing,
                "likely_merged": len(containing) > 0
            })
        else:
            results.append({
                "fn": f,
                "overlapping_preds": [],
                "likely_merged": False
            })
    
    return results

def analyze_fp_near_fn(fp: list[Variant], fn: list[Variant]):
    """Find FPs that are near FNs (potential position shifts)"""
    results = []
    
    for f in fp:
        nearby_fn = []
        for n in fn:
            dist = abs(f.pos - n.pos)
            if dist < 500 and abs(f.svlen - n.svlen) / max(f.svlen, n.svlen, 1) < 0.5:
                nearby_fn.append((n, dist))
        
        if nearby_fn:
            results.append({
                "fp": f,
                "nearby_fns": nearby_fn,
                "likely_shifted": True
            })
    
    return results

def size_bracket(svlen: int) -> str:
    if svlen <= 10:
        return "tiny (1-10bp)"
    elif svlen <= 50:
        return "small (11-50bp)"
    elif svlen <= 500:
        return "medium (51-500bp)"
    elif svlen <= 5000:
        return "large (501-5kb)"
    elif svlen <= 50000:
        return "very_large (5k-50kb)"
    else:
        return "mega (50kb+)"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--truth', required=True)
    parser.add_argument('--pred', required=True)
    parser.add_argument('--chrom', default=None)
    parser.add_argument('--type', default='DEL', choices=['DEL', 'INS', 'SNP'])
    args = parser.parse_args()
    
    truth_all = parse_vcf(args.truth, args.chrom)
    pred_all = parse_vcf(args.pred, args.chrom)
    
    truth = truth_all[args.type]
    pred = pred_all[args.type]
    
    print("="*80)
    print(f"DEEP ANALYSIS: {args.type}")
    print("="*80)
    print(f"Truth: {len(truth)}, Predicted: {len(pred)}")
    print()
    
    # Basic matching
    tp, fp, fn = match_variants(truth, pred)
    
    print(f"TP: {len(tp)}, FP: {len(fp)}, FN: {len(fn)}")
    precision = len(tp) / len(pred) if pred else 0
    recall = len(tp) / len(truth) if truth else 0
    print(f"Precision: {precision:.2%}, Recall: {recall:.2%}")
    print()
    
    # 1. SIZE BRACKET ANALYSIS
    print("-"*80)
    print("1. SIZE BRACKET ANALYSIS")
    print("-"*80)
    
    fn_by_size = defaultdict(list)
    for f in fn:
        fn_by_size[size_bracket(f.svlen)].append(f)
    
    tp_by_size = defaultdict(list)
    for t, p, _ in tp:
        tp_by_size[size_bracket(t.svlen)].append((t, p))
    
    truth_by_size = defaultdict(list)
    for t in truth:
        truth_by_size[size_bracket(t.svlen)].append(t)
    
    for bracket in ["tiny (1-10bp)", "small (11-50bp)", "medium (51-500bp)", 
                    "large (501-5kb)", "very_large (5k-50kb)", "mega (50kb+)"]:
        total = len(truth_by_size[bracket])
        missed = len(fn_by_size[bracket])
        found = len(tp_by_size[bracket])
        if total > 0:
            print(f"  {bracket}: {found}/{total} found ({found/total:.0%}), {missed} missed")
    print()
    
    # 2. CLUSTER ANALYSIS
    print("-"*80)
    print("2. CLUSTER ANALYSIS (nearby variants in truth)")
    print("-"*80)
    
    truth_clusters = find_clusters(truth, max_dist=5000)
    print(f"Found {len(truth_clusters)} clusters of nearby truth variants")
    
    cluster_fn_count = 0
    cluster_partial_count = 0
    
    for i, cluster in enumerate(truth_clusters[:10]):  # Show first 10
        cluster_fn = [c for c in cluster if c in fn]
        cluster_tp = [c for c in cluster if c in [t for t, _, _ in tp]]
        
        if cluster_fn:
            cluster_fn_count += 1
        if cluster_fn and cluster_tp:
            cluster_partial_count += 1
        
        print(f"\n  Cluster {i+1}: {len(cluster)} variants at {cluster[0].pos}-{cluster[-1].end}")
        print(f"    Sizes: {[c.svlen for c in cluster]}")
        print(f"    Found: {len(cluster_tp)}, Missed: {len(cluster_fn)}")
        
        if cluster_fn:
            # Check if there's a prediction covering this cluster
            cluster_start = min(c.pos for c in cluster)
            cluster_end = max(c.end for c in cluster)
            covering_preds = [p for p in pred if p.pos <= cluster_end + 500 and p.end >= cluster_start - 500]
            if covering_preds:
                print(f"    POTENTIAL MERGE: {len(covering_preds)} prediction(s) overlap cluster")
                for p in covering_preds[:3]:
                    print(f"      -> Pred at {p.pos}, len={p.svlen}")
    
    print(f"\n  Summary: {cluster_fn_count}/{len(truth_clusters)} clusters have missed variants")
    print()
    
    # 3. MERGE ANALYSIS
    print("-"*80)
    print("3. MERGE ANALYSIS (FNs potentially absorbed by larger calls)")
    print("-"*80)
    
    merge_analysis = analyze_fn_clusters(fn, pred)
    merged_count = sum(1 for m in merge_analysis if m["likely_merged"])
    
    print(f"FNs with overlapping predictions: {merged_count}/{len(fn)}")
    print("\nTop examples of likely merges:")
    
    for m in merge_analysis[:15]:
        if m["overlapping_preds"]:
            f = m["fn"]
            preds = m["overlapping_preds"]
            print(f"  FN at {f.pos} (len={f.svlen}) overlaps with:")
            for p in preds[:2]:
                print(f"    -> Pred at {p.pos} (len={p.svlen})")
    print()
    
    # 4. POSITION SHIFT ANALYSIS
    print("-"*80)
    print("4. POSITION SHIFT ANALYSIS (FPs near FNs)")
    print("-"*80)
    
    shift_analysis = analyze_fp_near_fn(fp, fn)
    shift_count = len(shift_analysis)
    
    print(f"FPs with nearby FNs (potential shifts): {shift_count}/{len(fp)}")
    print("\nTop examples of likely position shifts:")
    
    for s in shift_analysis[:10]:
        f = s["fp"]
        for n, dist in s["nearby_fns"]:
            print(f"  FP at {f.pos} (len={f.svlen}) <-> FN at {n.pos} (len={n.svlen}), dist={dist}bp")
    print()
    
    # 5. SIZE MISMATCH ANALYSIS (in TPs)
    print("-"*80)
    print("5. SIZE MISMATCH ANALYSIS (in matched variants)")
    print("-"*80)
    
    size_mismatches = [(t, p, info) for t, p, info in tp if abs(info["size_diff"]) > 10]
    print(f"TPs with size mismatch >10bp: {len(size_mismatches)}/{len(tp)}")
    for t, p, info in sorted(size_mismatches, key=lambda x: -abs(x[2]["size_diff"]))[:10]:
        print(f"  Truth {t.pos} (len={t.svlen}) vs Pred {p.pos} (len={p.svlen}), diff={info['size_diff']}")
    print()
    
    # 6. ISOLATED FNs (not explained by clusters or merges)
    print("-" * 80)
    print("7. TARGET SV ANALYSIS (>50bp)")
    print("-" * 80)
    
    target_truth = [v for v in truth_all[args.type] if abs(v.svlen) >= 50]
    target_pred = [v for v in pred_all[args.type] if abs(v.svlen) >= 50]
    
    target_tp = 0
    target_missed = []
    
    # Simple matching for target stats
    for t in target_truth:
        found = False
        for p in target_pred:
            if p.chrom == t.chrom and abs(p.pos - t.pos) < 500 and abs(abs(p.svlen) - abs(t.svlen)) < abs(t.svlen) * 0.5:
                found = True
                break
        if found:
            target_tp += 1
        else:
            target_missed.append(t)
            
    target_prec = 0.0
    if len(target_pred) > 0:
        target_prec = target_tp / len(target_pred)
        
    target_rec = 0.0
    if len(target_truth) > 0:
        target_rec = target_tp / len(target_truth)
        
    print(f"Target SVs (>50bp): Truth={len(target_truth)}, Pred={len(target_pred)}")
    print(f"TP: {target_tp}, Missed: {len(target_missed)}")
    print(f"Target Precision: {target_prec*100:.2f}%, Target Recall: {target_rec*100:.2f}%")
    
    if len(target_missed) > 0:
        print("\nMissed Target SVs:")
        for t in target_missed[:10]:
             print(f"  -> {t.pos} (len={t.svlen})")
        if len(target_missed) > 10:
            print(f"  ... and {len(target_missed)-10} more")

    print("-" * 80)
    print("8. ISOLATED FNs (unexplained misses)")
    print("-" * 80)
    
    isolated = [m for m in merge_analysis if not m["likely_merged"]]
    print(f"Isolated FNs (no overlapping prediction): {len(isolated)}/{len(fn)}")
    
    # Group by size
    isolated_by_size = defaultdict(list)
    for m in isolated:
        isolated_by_size[size_bracket(m["fn"].svlen)].append(m["fn"])
    
    for bracket, variants in isolated_by_size.items():
        print(f"  {bracket}: {len(variants)}")
        for v in variants[:3]:
            print(f"    -> {v.pos} (len={v.svlen})")
    print()
    
    # 7. SUMMARY & RECOMMENDATIONS
    print("="*80)
    print("SUMMARY & DIAGNOSIS")
    print("="*80)
    
    print(f"""
Total FNs: {len(fn)}

Breakdown:
  - Likely merged into larger calls: {merged_count} ({merged_count/len(fn)*100:.1f}%)
  - In dense clusters: {cluster_fn_count} clusters affected
  - Isolated misses: {len(isolated)} ({len(isolated)/len(fn)*100:.1f}%)

Top Issues:
1. VARIANT CLUSTERING: {len(truth_clusters)} regions have multiple nearby variants.
   The aligner may be struggling to resolve them correctly.

2. MERGE EFFECT: {merged_count} FNs have overlapping predictions, suggesting
   they are being absorbed into larger, imprecise calls.

3. SIZE BRACKET VULNERABILITY: Check the size bracket table above.
   Some size ranges may be systematically underperforming.

4. POSITION SHIFTS: {shift_count} FP/FN pairs suggest correct detection but
   incorrect positioning.
""")

if __name__ == '__main__':
    main()
