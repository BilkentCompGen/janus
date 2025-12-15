#!/usr/bin/env python3
"""
Merge overlapping deletion calls in TSV files.

This script reads a TSV file with deletion calls and merges overlapping
or adjacent intervals on the same chromosome with the same tumor_samples.
"""

import sys
import argparse
from collections import defaultdict
from typing import List, Tuple, Set, Dict, Optional


def parse_tsv_line(line: str) -> Tuple[str, str, int, int, str, int]:
    """Parse a TSV line and return (type, chrom, start, end, tumor_samples, count)."""
    fields = line.strip().split('\t')
    if len(fields) != 6:
        raise ValueError(f"Expected 6 fields, got {len(fields)}: {line}")
    
    type_str, chrom, start_str, end_str, samples_str, count_str = fields
    return (type_str, chrom, int(start_str), int(end_str), samples_str, int(count_str))


def intervals_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Check if two intervals overlap."""
    return start1 < end2 and end1 > start2


def load_collision_positions(collision_file: Optional[str]) -> Tuple[Dict[str, List[Tuple[int, int, int]]], Dict[Tuple[str, int, int], int]]:
    """
    Load collision positions from TSV file, indexed by chromosome.
    
    Returns:
        (collision_by_chrom, position_to_hash_map)
        collision_by_chrom: dict from chrom to list of (start, end, hash_decimal) tuples, sorted by start
        position_to_hash_map: map from (chrom, start, end) to hash_decimal
    """
    if not collision_file:
        return {}, {}
    
    collision_by_chrom = defaultdict(list)  # chrom -> [(start, end, hash), ...]
    position_to_hash = {}
    
    print(f"Reading collision file {collision_file}...", file=sys.stderr)
    with open(collision_file, 'r') as f:
        header = f.readline().strip()
        if not header.startswith('hash'):
            raise ValueError("Expected collision TSV header starting with 'hash'")
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            
            hash_hex, hash_decimal_str, collision_count_str, chrom, start_str, end_str = fields[:6]
            
            try:
                hash_decimal = int(hash_decimal_str)
                collision_count = int(collision_count_str)
                start = int(start_str)
                end = int(end_str)
                
                # Only track positions from hashes with collisions (count > 1)
                if collision_count > 1:
                    pos = (chrom, start, end)
                    collision_by_chrom[chrom].append((start, end, hash_decimal))
                    position_to_hash[pos] = hash_decimal
            except ValueError:
                continue
    
    # Sort collision positions by start for each chromosome (for efficient lookup)
    for chrom in collision_by_chrom:
        collision_by_chrom[chrom].sort(key=lambda x: x[0])
    
    total_collisions = sum(len(positions) for positions in collision_by_chrom.values())
    print(f"Loaded {total_collisions} collision positions from {len(set(position_to_hash.values()))} unique hashes", file=sys.stderr)
    return dict(collision_by_chrom), position_to_hash


def filter_collision_deletions(
    deletions: List[Tuple[int, int, int, str]],  # (start, end, count, chrom)
    collision_by_chrom: Dict[str, List[Tuple[int, int, int]]],  # chrom -> [(start, end, hash), ...]
    position_to_hash: Dict[Tuple[str, int, int], int]
) -> Tuple[List[Tuple[int, int, int, str]], List[Tuple[int, int, int, str]]]:
    """
    Filter collision-based deletions based on neighbor presence.
    
    If a collision position has neighbors (overlapping deletions), keep ALL deletions
    from that collision hash. If no neighbors, discard ALL deletions from that hash.
    
    Returns:
        (kept_deletions, discarded_deletions)
    """
    if not collision_by_chrom:
        # No collisions to filter, return all as kept
        return deletions, []
    
    # Group deletions by chromosome and sort
    deletions_by_chrom = defaultdict(list)
    for start, end, count, chrom in deletions:
        deletions_by_chrom[chrom].append((start, end, count, chrom))
    
    # For each chromosome, sort deletions
    for chrom in deletions_by_chrom:
        deletions_by_chrom[chrom].sort(key=lambda x: x[0])
    
    # Build exact match lookup for fast checking
    exact_collision_positions = set(position_to_hash.keys())
    
    # Track which hashes should be kept (have neighbors) or discarded (no neighbors)
    hash_decisions = {}  # hash_decimal -> True (keep) or False (discard)
    
    # Check each deletion for neighbors (only check collisions from same chromosome)
    for chrom, chrom_deletions in deletions_by_chrom.items():
        chrom_collisions = collision_by_chrom.get(chrom, [])
        
        for i, (start, end, count, _) in enumerate(chrom_deletions):
            pos = (chrom, start, end)
            
            # Check exact match first (fast)
            collision_hash = position_to_hash.get(pos)
            
            # If no exact match, check overlap with sorted collision list (only nearby ones)
            if collision_hash is None and chrom_collisions:
                # Binary search for nearby collisions (only check a small window)
                # Since collisions are sorted, we can do a limited search
                for coll_start, coll_end, coll_hash in chrom_collisions:
                    # Early termination if we've passed the relevant region
                    if coll_end < start:
                        continue
                    if coll_start > end:
                        break
                    if intervals_overlap(start, end, coll_start, coll_end):
                        collision_hash = coll_hash
                        break
            
            if collision_hash is None:
                continue
            
            # Check for neighbors (other deletions that overlap with this one)
            has_neighbors = False
            for j, (other_start, other_end, _, _) in enumerate(chrom_deletions):
                if i != j and intervals_overlap(start, end, other_start, other_end):
                    has_neighbors = True
                    break
            
            # Decision: if has neighbors, keep this hash; if no neighbors, discard
            if collision_hash not in hash_decisions:
                hash_decisions[collision_hash] = has_neighbors
            else:
                # If any position of this hash has neighbors, keep all
                if has_neighbors:
                    hash_decisions[collision_hash] = True
    
    # Apply filtering
    kept = []
    discarded = []
    
    for start, end, count, chrom in deletions:
        pos = (chrom, start, end)
        
        # Check if this position overlaps a collision (exact match first, then overlap)
        collision_hash = position_to_hash.get(pos)
        
        # If no exact match, check overlap with sorted collision list
        if collision_hash is None:
            chrom_collisions = collision_by_chrom.get(chrom, [])
            for coll_start, coll_end, coll_hash in chrom_collisions:
                # Early termination
                if coll_end < start:
                    continue
                if coll_start > end:
                    break
                if intervals_overlap(start, end, coll_start, coll_end):
                    collision_hash = coll_hash
                    break
        
        if collision_hash is not None:
            if hash_decisions.get(collision_hash, False):
                kept.append((start, end, count, chrom))
            else:
                discarded.append((start, end, count, chrom))
        else:
            # Not a collision, always keep
            kept.append((start, end, count, chrom))
    
    return kept, discarded


def merge_intervals(intervals: List[Tuple[int, int, int]]) -> List[Tuple[int, int, int]]:
    """
    Merge overlapping or adjacent intervals.
    
    Args:
        intervals: List of (start, end, count) tuples, sorted by start
    
    Returns:
        List of merged (start, end, count) tuples
    """
    if not intervals:
        return []
    
    merged = []
    current_start, current_end, current_count = intervals[0]
    
    for start, end, count in intervals[1:]:
        # Check if intervals overlap or are adjacent (end + 1 >= start)
        if start <= current_end + 1:
            # Merge: extend end and sum counts
            current_end = max(current_end, end)
            current_count += count
        else:
            # No overlap: save current and start new
            merged.append((current_start, current_end, current_count))
            current_start, current_end, current_count = start, end, count
    
    # Don't forget the last interval
    merged.append((current_start, current_end, current_count))
    
    return merged


def process_tsv(input_file: str, output_file: str, collision_file: Optional[str] = None, discarded_file: Optional[str] = None, min_length: int = 0, min_count: int = 0):
    """Process TSV file and merge overlapping deletions with optional collision filtering."""
    
    # Load collision positions if provided
    collision_by_chrom, position_to_hash = load_collision_positions(collision_file)
    
    # Group intervals by (chrom, tumor_samples)
    groups = defaultdict(list)
    header = None
    
    print(f"Reading {input_file}...", file=sys.stderr)
    with open(input_file, 'r') as f:
        # Read header
        header = f.readline().strip()
        if not header.startswith('type'):
            raise ValueError("Expected header line starting with 'type'")
        
        # Read all intervals
        line_num = 1
        for line in f:
            line_num += 1
            line = line.strip()
            if not line:
                continue
            
            try:
                type_str, chrom, start, end, samples, count = parse_tsv_line(line)
                
                # Only process DEL type
                if type_str == "DEL":
                    key = (chrom, samples)
                    groups[key].append((start, end, count, chrom))
                else:
                    # Non-DEL types: keep as-is (we'll write them separately if needed)
                    print(f"Warning: Skipping non-DEL type '{type_str}' at line {line_num}", file=sys.stderr)
            except Exception as e:
                print(f"Error parsing line {line_num}: {e}", file=sys.stderr)
                continue
    
    print(f"Found {len(groups)} unique (chrom, tumor_samples) groups", file=sys.stderr)
    
    # Process each group: filter collisions, then merge
    merged_results = []
    discarded_results = []
    total_before = 0
    total_after_merge = 0
    total_after_filter = 0
    total_discarded = 0
    
    for (chrom, samples), intervals in groups.items():
        total_before += len(intervals)
        
        # Extract just (start, end, count) for sorting and merging
        intervals_for_merge = [(start, end, count) for start, end, count, _ in intervals]
        intervals_for_filter = [(start, end, count, chrom) for start, end, count, chrom in intervals]
        
        # Apply collision filtering if collisions are provided
        if collision_by_chrom:
            kept_intervals, discarded_intervals = filter_collision_deletions(
                intervals_for_filter, collision_by_chrom, position_to_hash
            )
            total_discarded += len(discarded_intervals)
            
            # Convert back to (start, end, count) for merging
            intervals_for_merge = [(start, end, count) for start, end, count, _ in kept_intervals]
            
            # Store discarded for output
            for start, end, count, _ in discarded_intervals:
                discarded_results.append(("DEL", chrom, start, end, samples, count))
        else:
            kept_intervals = intervals_for_filter
        
        # Sort by start position
        intervals_for_merge.sort(key=lambda x: x[0])
        
        # Merge overlapping intervals
        merged = merge_intervals(intervals_for_merge)
        total_after_merge += len(merged)
        total_after_filter += len(merged)
        
        # Add to results
        for start, end, count in merged:
            merged_results.append(("DEL", chrom, start, end, samples, count))
    
    # Sort results: by chromosome (numeric if possible), then by start
    def sort_key(row):
        type_str, chrom, start, end, samples, count = row
        try:
            chrom_num = int(chrom)
        except ValueError:
            chrom_num = float('inf')  # Non-numeric chromosomes go last
        return (chrom_num, start)
    
    merged_results.sort(key=sort_key)
    discarded_results.sort(key=sort_key)
    
    # Apply length and count filtering
    filtered_by_length_count = []
    for type_str, chrom, start, end, samples, count in merged_results:
        length = end - start
        if length >= min_length and count >= min_count:
            filtered_by_length_count.append((type_str, chrom, start, end, samples, count))
    
    total_filtered = len(filtered_by_length_count)
    
    # Print statistics
    print(f"Merged {total_before} intervals into {total_after_merge} intervals", file=sys.stderr)
    if collision_by_chrom:
        print(f"After collision filtering: {total_after_filter} kept, {total_discarded} discarded", file=sys.stderr)
        print(f"Collision filtering removed {100 * total_discarded / total_before:.1f}% of intervals", file=sys.stderr)
    if min_length > 0 or min_count > 0:
        length_filtered = len(merged_results) - total_filtered
        print(f"After length/count filtering (min_length={min_length}bp, min_count={min_count}): {total_filtered} kept, {length_filtered} filtered", file=sys.stderr)
        print(f"Length/count filtering removed {100 * length_filtered / len(merged_results):.1f}% of merged intervals", file=sys.stderr)
    print(f"Final reduction: {total_before - total_filtered} intervals ({100 * (total_before - total_filtered) / total_before:.1f}%)", file=sys.stderr)
    
    # Write output (use filtered results)
    print(f"Writing merged results to {output_file}...", file=sys.stderr)
    with open(output_file, 'w') as f:
        f.write(header + '\n')
        for type_str, chrom, start, end, samples, count in filtered_by_length_count:
            f.write(f"{type_str}\t{chrom}\t{start}\t{end}\t{samples}\t{count}\n")
    
    print(f"Done! Wrote {total_filtered} merged intervals.", file=sys.stderr)
    
    # Write discarded deletions if requested
    if discarded_file and discarded_results:
        print(f"Writing discarded deletions to {discarded_file}...", file=sys.stderr)
        with open(discarded_file, 'w') as f:
            f.write(header + '\n')
            for type_str, chrom, start, end, samples, count in discarded_results:
                f.write(f"{type_str}\t{chrom}\t{start}\t{end}\t{samples}\t{count}\n")
        print(f"Wrote {len(discarded_results)} discarded intervals.", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Merge overlapping deletion calls in TSV files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.tsv output.tsv
  %(prog)s outputs/tsv/d4.tsv outputs/tsv/d4_merged.tsv
        """
    )
    parser.add_argument('input', help='Input TSV file')
    parser.add_argument('output', help='Output TSV file')
    parser.add_argument('--collisions', help='Collision TSV file (optional, for filtering)')
    parser.add_argument('--discarded', help='Output file for discarded deletions (optional)')
    parser.add_argument('--min-length', type=int, default=2000, help='Minimum deletion length in bp (default: 2000)')
    parser.add_argument('--min-count', type=int, default=0, help='Minimum supporting call count (default: 0)')
    
    args = parser.parse_args()
    
    try:
        process_tsv(args.input, args.output, args.collisions, args.discarded, args.min_length, args.min_count)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
