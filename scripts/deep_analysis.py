#!/usr/bin/env python3
"""
Deep Analysis Script for SVP Variant Calling Accuracy

This script performs comprehensive analysis of prediction vs ground truth variants,
identifying patterns in errors, clustering of false positives, and potential causes.

Usage:
    python3 deep_analysis.py --truth <truth.vcf> --pred <pred.vcf> --output <output_dir>
"""

import argparse
import re
import os
import json
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
import statistics

# ==============================================================================
# Data Structures
# ==============================================================================

@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    svtype: str
    svlen: int
    end: int = 0
    source: str = ""  # "truth" or "pred"
    matched: bool = False
    match_type: str = ""  # "exact", "overlap", "missed", "fp", "merged"
    match_partner: Optional['Variant'] = None
    merged_partners: List['Variant'] = field(default_factory=list)
    
    def __post_init__(self):
        if self.end == 0:
            if self.svtype == "DEL":
                self.end = self.pos + abs(self.svlen)
            else:
                self.end = self.pos + 1
    
    def overlaps(self, other: 'Variant', tolerance: int = 0) -> bool:
        """Check if two variants overlap (for DEL) or are within tolerance (for INS/SNP)"""
        if self.chrom != other.chrom:
            return False
        if self.svtype != other.svtype:
            return False
        
        if self.svtype == "DEL":
            return not (self.end + tolerance < other.pos or other.end + tolerance < self.pos)
        else:  # INS or SNP
            return abs(self.pos - other.pos) <= tolerance
    
    def reciprocal_overlap(self, other: 'Variant') -> float:
        """Calculate reciprocal overlap for DELs"""
        if self.svtype != "DEL" or other.svtype != "DEL":
            return 0.0
        overlap_start = max(self.pos, other.pos)
        overlap_end = min(self.end, other.end)
        overlap_len = max(0, overlap_end - overlap_start)
        
        len_self = self.end - self.pos
        len_other = other.end - other.pos
        
        if len_self == 0 or len_other == 0:
            return 0.0
        
        return min(overlap_len / len_self, overlap_len / len_other)
    
    def size_similarity(self, other: 'Variant') -> float:
        """Calculate size similarity"""
        if self.svlen == 0 and other.svlen == 0:
            return 1.0
        len_self = abs(self.svlen) if self.svlen != 0 else 1
        len_other = abs(other.svlen) if other.svlen != 0 else 1
        return min(len_self, len_other) / max(len_self, len_other)

# ==============================================================================
# Parsing Functions
# ==============================================================================

def parse_vcf(filepath: str, source: str) -> List[Variant]:
    """Parse VCF file and return list of Variant objects"""
    variants = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]
            pos = int(pos)
            
            # Parse INFO field
            info_dict = {}
            for item in info.split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    info_dict[k] = v
            
            # Determine svtype and svlen
            svtype = info_dict.get('SVTYPE', '')
            svlen = 0
            
            if 'SVLEN' in info_dict:
                svlen = int(info_dict['SVLEN'])
            
            if alt == '<DEL>' or svtype == 'DEL':
                svtype = 'DEL'
                if svlen == 0 and 'END' in info_dict:
                    svlen = int(info_dict['END']) - pos
            elif alt == '<INS>' or svtype == 'INS':
                svtype = 'INS'
            elif len(ref) == 1 and len(alt) == 1 and ref != alt:
                svtype = 'SNP'
                svlen = 1
            elif len(ref) > len(alt):
                svtype = 'DEL'
                svlen = len(ref) - len(alt)
            elif len(alt) > len(ref):
                svtype = 'INS'
                svlen = len(alt) - len(ref)
            else:
                continue  # Skip complex or unknown
            
            if svtype:
                v = Variant(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    svtype=svtype,
                    svlen=abs(svlen) if svlen else 1,
                    source=source
                )
                variants.append(v)
    
    return variants

# ==============================================================================
# Analysis Functions
# ==============================================================================

def analyze_fp_clustering(fps: List[Variant], window: int = 100) -> Dict:
    """Analyze clustering of false positives"""
    if not fps:
        return {
            "clusters": [], 
            "max_cluster_size": 0, 
            "chatter_regions": [],
            "num_clusters": 0,
            "fp_in_chatter": 0,
            "total_fp": 0
        }
    
    fps_sorted = sorted(fps, key=lambda v: (v.chrom, v.pos))
    clusters = []
    current_cluster = [fps_sorted[0]]
    
    for v in fps_sorted[1:]:
        if v.chrom == current_cluster[-1].chrom and v.pos - current_cluster[-1].pos <= window:
            current_cluster.append(v)
        else:
            if len(current_cluster) >= 3:
                clusters.append(current_cluster)
            current_cluster = [v]
    
    if len(current_cluster) >= 3:
        clusters.append(current_cluster)
    
    # Find chatter regions (extremely dense FP clusters)
    chatter_regions = []
    for cluster in clusters:
        start = cluster[0].pos
        end = cluster[-1].pos
        span = end - start + 1
        density = len(cluster) / span if span > 0 else 0
        
        if density > 0.1 or len(cluster) > 50:  # High density or large clusters
            chatter_regions.append({
                "chrom": cluster[0].chrom,
                "start": start,
                "end": end,
                "count": len(cluster),
                "density": density,
                "span": span
            })
    
    return {
        "clusters": clusters,
        "num_clusters": len(clusters),
        "max_cluster_size": max(len(c) for c in clusters) if clusters else 0,
        "chatter_regions": chatter_regions,
        "fp_in_chatter": sum(r["count"] for r in chatter_regions),
        "total_fp": len(fps)
    }

def analyze_missed_variants(missed: List[Variant], all_truth: List[Variant]) -> Dict:
    """Analyze why variants might have been missed"""
    if not missed:
        return {"patterns": [], "nearby_variants": []}
    
    patterns = defaultdict(list)
    
    for m in missed:
        # Check for size patterns
        if abs(m.svlen) > 50000:
            patterns["very_large_sv"].append(m)
        elif abs(m.svlen) > 10000:
            patterns["large_sv"].append(m)
        elif abs(m.svlen) > 1000:
            patterns["medium_sv"].append(m)
        elif abs(m.svlen) > 50:
            patterns["small_sv"].append(m)
        else:
            patterns["tiny_sv"].append(m)
        
        # Check if there are nearby truth variants (potential collision/masking)
        nearby = []
        for t in all_truth:
            if t.chrom == m.chrom and t != m:
                dist = abs(t.pos - m.pos)
                if dist < 10000 and dist > 0:  # Within 10kb
                    nearby.append({
                        "pos": t.pos,
                        "svtype": t.svtype,
                        "svlen": t.svlen,
                        "distance": dist
                    })
        
        if nearby:
            patterns["has_nearby_variants"].append({
                "variant": m,
                "nearby": sorted(nearby, key=lambda x: x["distance"])[:5]
            })
    
    return {
        "by_size": {k: [(v.pos, v.svlen, v.svtype) for v in vlist] for k, vlist in patterns.items() if k != "has_nearby_variants"},
        "with_nearby": patterns.get("has_nearby_variants", []),
        "size_distribution": {
            "very_large_sv": len(patterns["very_large_sv"]),
            "large_sv": len(patterns["large_sv"]),
            "medium_sv": len(patterns["medium_sv"]),
            "small_sv": len(patterns["small_sv"]),
            "tiny_sv": len(patterns["tiny_sv"])
        }
    }

def analyze_length_mismatches(matched_pairs: List[Tuple[Variant, Variant]]) -> Dict:
    """Analyze length mismatches between matched variants"""
    mismatches = []
    
    for truth, pred in matched_pairs:
        if truth.svlen > 0 and pred.svlen > 0:
            ratio = pred.svlen / truth.svlen
            diff = abs(pred.svlen - truth.svlen)
            
            if ratio < 0.5 or ratio > 2.0 or diff > 1000:
                mismatches.append({
                    "truth_pos": truth.pos,
                    "truth_len": truth.svlen,
                    "pred_pos": pred.pos,
                    "pred_len": pred.svlen,
                    "ratio": ratio,
                    "diff": diff,
                    "svtype": truth.svtype
                })
    
    return {
        "severe_mismatches": sorted(mismatches, key=lambda x: -abs(x["diff"]))[:20],
        "count": len(mismatches),
        "avg_ratio": statistics.mean([m["ratio"] for m in mismatches]) if mismatches else 1.0
    }

def find_overlapping_truth_events(variants: List[Variant]) -> List[Dict]:
    """Find truth events that overlap each other (complex regions)"""
    overlapping = []
    
    for i, v1 in enumerate(variants):
        for j, v2 in enumerate(variants[i+1:], i+1):
            if v1.chrom == v2.chrom:
                # Check for overlap
                if v1.svtype == "DEL" and v2.svtype == "DEL":
                    if not (v1.end < v2.pos or v2.end < v1.pos):
                        overlapping.append({
                            "v1": (v1.pos, v1.svlen, v1.svtype),
                            "v2": (v2.pos, v2.svlen, v2.svtype),
                            "overlap_type": "DEL-DEL"
                        })
                elif abs(v1.pos - v2.pos) < 1000:
                    overlapping.append({
                        "v1": (v1.pos, v1.svlen, v1.svtype),
                        "v2": (v2.pos, v2.svlen, v2.svtype),
                        "overlap_type": f"{v1.svtype}-{v2.svtype}"
                    })
    
    return overlapping

def check_merge_matches(truth: List[Variant], pred: List[Variant], results: Dict) -> int:
    """
    Check if a single prediction covers multiple ground truth variants (Merge Match).
    Returns number of merge matches found.
    
    Validation criteria:
    - Multiple GT events must be clustered (gaps <= 500bp)
    - Total GT length must fill >= 75% of prediction length
    """
    merge_count = 0
    
    # Group truth by type
    truth_by_chrom = defaultdict(list)
    for t in truth:
        if not t.matched and t.svtype == "DEL":
            truth_by_chrom[t.chrom].append(t)
            
    for p in pred:
        if p.matched or p.svtype != "DEL":
            continue
            
        # Find all GT variants this prediction overlaps
        candidates = []
        if p.chrom in truth_by_chrom:
            for t in truth_by_chrom[p.chrom]:
                if not t.matched and t.overlaps(p, tolerance=100):
                    candidates.append(t)
        
        # If it covers multiple close variants
        if len(candidates) >= 2:
            # Check if candidates are close to each other
            candidates.sort(key=lambda x: x.pos)
            is_valid_merge = True
            
            # Verify candidates are clustered
            for i in range(len(candidates)-1):
                if candidates[i+1].pos - candidates[i].end > 500:
                    is_valid_merge = False
                    break
            
            if is_valid_merge:
                # Calculate fill ratio: how much of prediction is explained by GTs
                total_gt_len = sum(abs(t.svlen) for t in candidates)
                pred_len = abs(p.svlen)
                fill_ratio = total_gt_len / pred_len if pred_len > 0 else 0
                
                # Require at least 75% fill
                if fill_ratio >= 0.75:
                    p.matched = True
                    p.match_type = "merged"
                    p.merged_partners = candidates
                    
                    for t in candidates:
                        t.matched = True
                        t.match_type = "merged_into"
                        t.match_partner = p
                    
                    results["DEL"]["tp"].append(p)
                    merge_count += 1

    return merge_count

def match_variants(truth: List[Variant], pred: List[Variant], tolerance: int = 50) -> Dict:
    """Match truth and predicted variants"""
    results = {
        "DEL": {"tp": [], "fp": [], "fn": [], "matched_pairs": [], "merges": []},
        "INS": {"tp": [], "fp": [], "fn": [], "matched_pairs": [], "merges": []},
        "SNP": {"tp": [], "fp": [], "fn": [], "matched_pairs": [], "merges": []}
    }
    
    truth_by_type = defaultdict(list)
    pred_by_type = defaultdict(list)
    
    for v in truth:
        truth_by_type[v.svtype].append(v)
    for v in pred:
        pred_by_type[v.svtype].append(v)
    
    for svtype in ["DEL", "INS", "SNP"]:
        truth_list = truth_by_type[svtype]
        pred_list = pred_by_type[svtype]
        
        used_pred = set()
        
        for t in truth_list:
            best_match = None
            best_score = 0
            
            for i, p in enumerate(pred_list):
                if i in used_pred:
                    continue
                
                if svtype == "DEL":
                    ro = t.reciprocal_overlap(p)
                    if ro > best_score:
                        best_score = ro
                        best_match = (i, p)
                elif svtype == "SNP":
                    if t.chrom == p.chrom and t.pos == p.pos:
                        best_score = 1.0
                        best_match = (i, p)
                else:  # INS
                    # Proportional tolerance: max(10, 30% of smaller variant size)
                    min_len = min(abs(t.svlen), abs(p.svlen)) if t.svlen and p.svlen else 50
                    prop_tolerance = max(10, int(min_len * 0.3))
                    
                    if t.chrom == p.chrom and abs(t.pos - p.pos) <= prop_tolerance:
                        size_sim = t.size_similarity(p)
                        if size_sim > best_score:
                            best_score = size_sim
                            best_match = (i, p)
            
            if best_match and best_score >= 0.5:
                used_pred.add(best_match[0])
                t.matched = True
                best_match[1].matched = True
                results[svtype]["tp"].append(t)
                results[svtype]["matched_pairs"].append((t, best_match[1]))
            else:
                results[svtype]["fn"].append(t)
        
        for i, p in enumerate(pred_list):
            if i not in used_pred:
                results[svtype]["fp"].append(p)
    
    # Check for Merge Matches (Many-to-One)
    check_merge_matches(truth, pred, results)
    
    return results

# ==============================================================================
# Report Generation
# ==============================================================================

def generate_report(truth_path: str, pred_path: str, output_dir: str):
    """Generate comprehensive analysis report"""
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 60)
    print("SVP Deep Analysis Report")
    print("=" * 60)
    
    # Parse files
    print("\n[1/7] Parsing VCF files...")
    truth = parse_vcf(truth_path, "truth")
    pred = parse_vcf(pred_path, "pred")
    
    print(f"    Truth variants: {len(truth)}")
    print(f"    Pred variants: {len(pred)}")
    
    # Match variants
    print("\n[2/7] Matching variants...")
    results = match_variants(truth, pred)
    
    # Summary stats
    print("\n[3/7] Computing summary statistics...")
    summary = {}
    for svtype in ["DEL", "INS", "SNP"]:
        tp = len(results[svtype]["tp"])
        fp = len(results[svtype]["fp"])
        fn = len(results[svtype]["fn"])
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        summary[svtype] = {
            "TP": tp, "FP": fp, "FN": fn,
            "Precision": round(precision, 4),
            "Recall": round(recall, 4),
            "F1": round(f1, 4)
        }
        print(f"    {svtype}: TP={tp}, FP={fp}, FN={fn}, P={precision:.3f}, R={recall:.3f}, F1={f1:.3f}")
        
        # Report Merges
        merges = [p for p in results[svtype]["tp"] if p.match_type == "merged"]
        if merges:
            print(f"    + Detected {len(merges)} merged {svtype} events (One Pred covering Multiple GTs):")
            for m in merges[:5]:
                print(f"      Pred {m.chrom}:{m.pos} (len={m.svlen}) covers:")
                for t in m.merged_partners:
                    print(f"        -> GT {t.chrom}:{t.pos} (len={t.svlen})")
    
    # Analyze FP clustering (SNP)
    print("\n[4/7] Analyzing FP clustering (SNP)...")
    snp_fp_analysis = analyze_fp_clustering(results["SNP"]["fp"])
    print(f"    Total SNP FPs: {snp_fp_analysis['total_fp']}")
    print(f"    Num clusters (>=3 FPs each): {snp_fp_analysis['num_clusters']}")
    print(f"    Max cluster size: {snp_fp_analysis['max_cluster_size']}")
    print(f"    FPs in chatter regions: {snp_fp_analysis['fp_in_chatter']}")
    
    if snp_fp_analysis['chatter_regions']:
        print("\n    Top chatter regions:")
        for i, r in enumerate(sorted(snp_fp_analysis['chatter_regions'], key=lambda x: -x['count'])[:5]):
            print(f"      {i+1}. {r['chrom']}:{r['start']}-{r['end']} ({r['count']} FPs, span={r['span']}bp, density={r['density']:.4f})")
    
    # Analyze missed variants
    print("\n[5/7] Analyzing missed variants...")
    for svtype in ["DEL", "INS"]:
        missed = results[svtype]["fn"]
        if missed:
            print(f"\n    {svtype} - Completely Missed ({len(missed)}):")
            miss_analysis = analyze_missed_variants(missed, [v for v in truth if v.svtype == svtype])
            
            for size_cat, count in miss_analysis["size_distribution"].items():
                if count > 0:
                    print(f"      {size_cat}: {count}")
            
            if miss_analysis["with_nearby"]:
                print(f"\n    {svtype} variants with nearby truth variants:")
                for item in miss_analysis["with_nearby"][:5]:
                    v = item["variant"]
                    print(f"      {v.chrom}:{v.pos} (len={v.svlen}) - {len(item['nearby'])} nearby variants")
                    for n in item["nearby"][:3]:
                        print(f"        -> {n['svtype']} at pos {n['pos']} (len={n['svlen']}, dist={n['distance']})")
    
    # Analyze length mismatches
    print("\n[6/7] Analyzing length mismatches...")
    for svtype in ["DEL", "INS"]:
        mismatch_analysis = analyze_length_mismatches(results[svtype]["matched_pairs"])
        if mismatch_analysis["severe_mismatches"]:
            print(f"\n    {svtype} - Severe length mismatches ({mismatch_analysis['count']}):")
            for m in mismatch_analysis["severe_mismatches"][:5]:
                print(f"      pos={m['truth_pos']}: truth_len={m['truth_len']} vs pred_len={m['pred_len']} (ratio={m['ratio']:.2f})")
    
    # Find overlapping truth events
    print("\n[7/7] Analyzing overlapping truth events...")
    overlapping = find_overlapping_truth_events(truth)
    if overlapping:
        print(f"    Found {len(overlapping)} overlapping truth event pairs")
        overlap_types = defaultdict(int)
        for o in overlapping:
            overlap_types[o["overlap_type"]] += 1
        for ot, count in overlap_types.items():
            print(f"      {ot}: {count}")
    
    # Save detailed JSON report
    report = {
        "summary": summary,
        "snp_fp_clusters": {
            "total_fp": snp_fp_analysis['total_fp'],
            "num_clusters": snp_fp_analysis['num_clusters'],
            "max_cluster_size": snp_fp_analysis['max_cluster_size'],
            "chatter_regions": snp_fp_analysis['chatter_regions']
        },
        "missed_del": [(v.pos, v.svlen) for v in results["DEL"]["fn"]],
        "missed_ins": [(v.pos, v.svlen) for v in results["INS"]["fn"]],
        "overlapping_truth": overlapping[:50]
    }
    
    report_path = os.path.join(output_dir, "deep_analysis.json")
    with open(report_path, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\n\nDetailed report saved to: {report_path}")
    
    # Generate position-based analysis for missed variants
    print("\n" + "=" * 60)
    print("DETAILED MISSED VARIANT ANALYSIS")
    print("=" * 60)
    
    # Group missed by potential cause
    print("\n### Missed DEL Analysis:")
    missed_dels = sorted(results["DEL"]["fn"], key=lambda v: -abs(v.svlen))
    
    # Check for clusters of missed variants (suggests regional issue)
    missed_del_positions = [v.pos for v in missed_dels]
    if missed_del_positions:
        missed_del_positions.sort()
        
        clusters = []
        current = [missed_del_positions[0]]
        for pos in missed_del_positions[1:]:
            if pos - current[-1] < 100000:  # Within 100kb
                current.append(pos)
            else:
                if len(current) >= 2:
                    clusters.append(current)
                current = [pos]
        if len(current) >= 2:
            clusters.append(current)
        
        if clusters:
            print("\nMissed DEL clusters (suggests problematic regions):")
            for cluster in clusters:
                print(f"  Region ~{cluster[0]}-{cluster[-1]}: {len(cluster)} missed variants")
    
    # Look for position collision patterns
    print("\n### Position Collision Analysis:")
    truth_positions = {(v.chrom, v.pos): v for v in truth}
    pred_positions = {(v.chrom, v.pos): v for v in pred}
    
    collision_issues = []
    for t in truth:
        key = (t.chrom, t.pos)
        if key in pred_positions:
            p = pred_positions[key]
            if t.svtype != p.svtype:
                collision_issues.append({
                    "pos": t.pos,
                    "truth_type": t.svtype,
                    "pred_type": p.svtype,
                    "truth_len": t.svlen,
                    "pred_len": p.svlen
                })
    
    if collision_issues:
        print(f"\nType collisions at same position ({len(collision_issues)}):")
        for c in collision_issues[:10]:
            print(f"  pos={c['pos']}: truth={c['truth_type']}(len={c['truth_len']}) vs pred={c['pred_type']}(len={c['pred_len']})")
    
    print("\n" + "=" * 60)
    print("RECOMMENDATIONS")
    print("=" * 60)
    
    # Generate recommendations based on analysis
    recommendations = []
    
    if snp_fp_analysis['fp_in_chatter'] > snp_fp_analysis['total_fp'] * 0.5:
        recommendations.append(
            "HIGH PRIORITY: Most SNP FPs ({:.1f}%) come from chatter regions. "
            "Consider implementing chatter detection and fallback to size-based calling, "
            "or filtering SNPs from regions with high variant density.".format(
                100 * snp_fp_analysis['fp_in_chatter'] / max(1, snp_fp_analysis['total_fp'])
            )
        )
    
    if len(results["DEL"]["fn"]) > 0:
        # Check if many missed are in clusters
        if clusters and sum(len(c) for c in clusters) > len(missed_dels) * 0.3:
            recommendations.append(
                "MEDIUM PRIORITY: Multiple missed DELs cluster together, suggesting "
                "problematic genomic regions. Consider: (1) checking bubble detection in these "
                "regions, (2) adjusting LCP depth, or (3) implementing special handling for complex regions."
            )
    
    large_missed = [v for v in results["DEL"]["fn"] if v.svlen > 5000]
    if large_missed:
        recommendations.append(
            f"MEDIUM PRIORITY: {len(large_missed)} large DELs (>5kb) were completely missed. "
            "These should be detectable via size-based calling. Check if the bubbles are being detected "
            "and extracted correctly in these regions."
        )
    
    for rec in recommendations:
        print(f"\n• {rec}")
    
    print("\n" + "=" * 60)

# ==============================================================================
# Main
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description='Deep analysis of SVP variant calling accuracy')
    parser.add_argument('--truth', '-t', required=True, help='Truth VCF file')
    parser.add_argument('--pred', '-p', required=True, help='Predicted VCF file')
    parser.add_argument('--output', '-o', default='analysis_output', help='Output directory')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.truth):
        print(f"Error: Truth file not found: {args.truth}")
        return
    if not os.path.exists(args.pred):
        print(f"Error: Pred file not found: {args.pred}")
        return
    
    generate_report(args.truth, args.pred, args.output)

if __name__ == "__main__":
    main()
