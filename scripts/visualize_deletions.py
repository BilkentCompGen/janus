#!/usr/bin/env python3
"""
Visualize chromosome with ground truth deletions marked.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import numpy as np
from io import StringIO

# Chromosome 19 length (GRCh38)
CHR19_LENGTH = 58_617_616


def load_ground_truth(tsv_file, parental=None):
    """Load ground truth deletions."""
    with open(tsv_file, 'r') as f:
        lines = f.readlines()
    
    if lines and lines[0].startswith('#'):
        lines[0] = lines[0][1:]
    
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')
    df = df[df['copy'] == -1].copy()
    
    if parental is not None:
        df = df[df['parental'] == parental].copy()
    
    return df


def load_detected(tsv_file):
    """Load detected deletions."""
    df = pd.read_csv(tsv_file, sep='\t')
    df = df[df['type'] == 'DEL'].copy()
    return df


def main():
    parser = argparse.ArgumentParser(description='Visualize deletions on chromosome')
    parser.add_argument('--ground-truth', default='simulation/ground_truth_cnvs.tsv',
                        help='Ground truth TSV file')
    parser.add_argument('--detected', help='Detected deletions TSV (optional)')
    parser.add_argument('--parental', type=int, choices=[0, 1], default=None,
                        help='Filter to parental haplotype')
    parser.add_argument('--carrier', type=int, help='Filter to specific carrier')
    parser.add_argument('--output', '-o', default='output/images/deletions_map.png',
                        help='Output image file')
    parser.add_argument('--show-all-carriers', action='store_true',
                        help='Show all carriers in separate tracks')
    
    args = parser.parse_args()
    
    # Load ground truth
    gt_df = load_ground_truth(args.ground_truth, parental=args.parental)
    
    if args.carrier:
        gt_df = gt_df[gt_df['carrier'] == args.carrier]
    
    # Get unique carriers (include all 1-10 if showing detected)
    carriers = sorted(gt_df['carrier'].unique())
    if args.detected and args.show_all_carriers:
        carriers = list(range(1, 11))  # All carriers 1-10
    
    # Color map for carriers (consistent colors for 1-10)
    cmap = plt.cm.tab10
    carrier_colors = {i: cmap((i-1) % 10) for i in range(1, 11)}
    
    # Create figure
    if args.show_all_carriers:
        n_tracks = len(carriers) + 1  # +1 for chromosome backbone
        fig, ax = plt.subplots(figsize=(16, 2 + n_tracks * 0.8))
    else:
        fig, ax = plt.subplots(figsize=(16, 4))
    
    # Draw chromosome backbone
    chr_y = 0
    chr_height = 0.3
    backbone = mpatches.FancyBboxPatch(
        (0, chr_y - chr_height/2), CHR19_LENGTH, chr_height,
        boxstyle="round,pad=0.01,rounding_size=1000000",
        facecolor='lightgray', edgecolor='black', linewidth=1
    )
    ax.add_patch(backbone)
    
    # Add chromosome label
    ax.text(-CHR19_LENGTH * 0.02, chr_y, 'chr19', fontsize=12, va='center', ha='right', fontweight='bold')
    
    if args.show_all_carriers:
        # Draw each carrier as a separate track
        for i, carrier in enumerate(carriers):
            track_y = -(i + 1) * 0.8
            carrier_dels = gt_df[gt_df['carrier'] == carrier]
            
            # Track label - show count
            n_gt = len(carrier_dels)
            ax.text(-CHR19_LENGTH * 0.02, track_y, f'Carrier {carrier} ({n_gt} GT)', 
                    fontsize=9, va='center', ha='right')
            
            # Draw ground truth deletions
            for _, row in carrier_dels.iterrows():
                rect = mpatches.Rectangle(
                    (row['start'], track_y - 0.2), 
                    row['end'] - row['start'], 0.4,
                    facecolor=carrier_colors.get(carrier, 'gray'), 
                    edgecolor='black', linewidth=0.5, alpha=0.7
                )
                ax.add_patch(rect)
        
        y_min = -(len(carriers) + 1) * 0.8
    else:
        # Draw all deletions on chromosome with offset by carrier
        for i, carrier in enumerate(carriers):
            carrier_dels = gt_df[gt_df['carrier'] == carrier]
            offset = 0.4 + i * 0.15
            
            for _, row in carrier_dels.iterrows():
                rect = mpatches.Rectangle(
                    (row['start'], chr_y + offset), 
                    row['end'] - row['start'], 0.12,
                    facecolor=carrier_colors[carrier], 
                    edgecolor='black', linewidth=0.3, alpha=0.8
                )
                ax.add_patch(rect)
        
        y_min = -0.5
    
    # Add detected deletions if provided
    if args.detected:
        det_df = load_detected(args.detected)
        
        if args.show_all_carriers:
            # Show detected deletions per carrier, side-by-side with ground truth
            for i, carrier in enumerate(carriers):
                track_y = -(i + 1) * 0.8
                carrier_det = det_df[det_df['tumor_samples'].astype(str) == str(carrier)]
                
                # Draw detected deletions slightly offset (darker, hatched)
                for _, row in carrier_det.iterrows():
                    rect = mpatches.Rectangle(
                        (row['start'], track_y - 0.35), 
                        row['end'] - row['start'], 0.3,
                        facecolor='none', edgecolor='darkred', linewidth=1.5, 
                        linestyle='--', alpha=0.9
                    )
                    ax.add_patch(rect)
        else:
            if args.carrier:
                det_df = det_df[det_df['tumor_samples'].astype(str) == str(args.carrier)]
            
            # Draw detected as a separate track below
            det_y = y_min - 0.8
            ax.text(-CHR19_LENGTH * 0.02, det_y, 'Detected', fontsize=10, va='center', ha='right', color='red')
            
            for _, row in det_df.iterrows():
                rect = mpatches.Rectangle(
                    (row['start'], det_y - 0.25), 
                    row['end'] - row['start'], 0.5,
                    facecolor='red', edgecolor='darkred', linewidth=0.5, alpha=0.5
                )
                ax.add_patch(rect)
            
            y_min = det_y - 0.8
    
    # Formatting
    ax.set_xlim(-CHR19_LENGTH * 0.08, CHR19_LENGTH * 1.02)
    ax.set_ylim(y_min, 1.5)
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_title(f'Ground Truth Deletions on chr19' + 
                 (f' (parental={args.parental})' if args.parental is not None else ''), 
                 fontsize=14, fontweight='bold')
    
    # Format x-axis with Mb
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x/1e6:.0f} Mb'))
    
    # Remove y-axis
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Legend
    legend_patches = [mpatches.Patch(color=carrier_colors[c], label=f'Carrier {c}') 
                      for c in carriers]
    ax.legend(handles=legend_patches, loc='upper right', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=150, bbox_inches='tight')
    print(f"Saved to {args.output}")
    plt.close()


if __name__ == '__main__':
    main()
