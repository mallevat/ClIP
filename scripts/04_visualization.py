#!/usr/bin/env python3
"""
================================================================================
CLIP Pipeline - Step 4: Visualization
================================================================================

Generates publication-ready figures for CLIP analysis:
    1. Rolling Gaussian heatmaps by zone
    2. Chord network diagrams showing cell-cell communication
    3. Cross-reference panel (L-R signal + cell-cell co-localization)

Output:
    - Heatmap_{zone}.pdf/png/svg
    - ChordNetwork.pdf/png/svg
    - CrossReference_{LR_pair}.pdf/png/svg

================================================================================
"""

import pandas as pd
import numpy as np
import scanpy as sc
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe
import matplotlib.colors as mcolors
import os
import argparse
import warnings
warnings.filterwarnings('ignore')

# Set up matplotlib for publication-quality figures
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'

# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_N_EVAL_POINTS = 150
DEFAULT_SIGMA = 0.05

# Zone colors
ZONE_COLORS = {
    'Tumor': '#1976D2',      # Blue
    'Interface': '#7B1FA2',  # Purple
    'Lymphoid': '#FF6B35',   # Orange
}

ZONE_BOUNDS = {
    'Tumor': (0.0, 0.36),
    'Interface': (0.36, 0.725),
    'Lymphoid': (0.725, 1.0)
}

# Cell type colors
CELLTYPE_COLORS = {
    'Epithelial.cells': '#A3D3E3',
    'Mast.cells': '#6E540C',
    'Dendritic.cells': '#E6800B',
    'Endothelial.cells': '#A3DFE3',
    'B_Plasma.cells': '#AD0003',
    'Malignant.cells': '#3160AD',
    'T.cells': '#F200D2',
    'Macrophages': '#CA6BCF',
    'Fibroblasts': '#FCD200',
    'Myocytes': '#1ECBD4',
}

CELLTYPE_SHORT = {
    'Epithelial.cells': 'Epi',
    'Mast.cells': 'Mast',
    'Dendritic.cells': 'DC',
    'Endothelial.cells': 'Endo',
    'B_Plasma.cells': 'B/Plasma',
    'Malignant.cells': 'Malig',
    'T.cells': 'T cell',
    'Macrophages': 'Mac',
    'Fibroblasts': 'Fibro',
    'Myocytes': 'Myo',
}

CELL_ABBREV = {
    'T.cells': 'T', 'Macrophages': 'Mac', 'Dendritic.cells': 'DC',
    'B_Plasma.cells': 'B', 'Fibroblasts': 'Fib', 'Malignant.cells': 'Mal',
    'Endothelial.cells': 'Endo', 'Epithelial.cells': 'Epi',
    'Mast.cells': 'Mast', 'Myocytes': 'Myo',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='CLIP Step 4: Visualization'
    )
    parser.add_argument('--clip_results', type=str, required=True,
                        help='CLIP_scores.csv from Step 3')
    parser.add_argument('--commot_dir', type=str, required=True,
                        help='Directory containing COMMOT results')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for figures')
    parser.add_argument('--samples', nargs='+', required=True,
                        help='Sample names')
    parser.add_argument('--n_eval_points', type=int, default=DEFAULT_N_EVAL_POINTS,
                        help=f'Trajectory evaluation points (default: {DEFAULT_N_EVAL_POINTS})')
    parser.add_argument('--sigma', type=float, default=DEFAULT_SIGMA,
                        help=f'Rolling Gaussian bandwidth (default: {DEFAULT_SIGMA})')
    parser.add_argument('--formats', nargs='+', default=['pdf', 'png', 'svg'],
                        help='Output formats (default: pdf png svg)')
    return parser.parse_args()


def abbrev_cell(cell):
    if pd.isna(cell):
        return '?'
    return CELL_ABBREV.get(cell, cell[:4])


def compute_rolling_gaussian(lr_pair, sample_data, eval_points, sigma):
    """Compute rolling Gaussian-smoothed signal."""
    all_traj, all_signal = [], []
    
    for sample, data in sample_data.items():
        col_name = f's-{lr_pair}'
        if col_name not in data['sender'].columns:
            continue
        all_traj.extend(data['trajectory'])
        all_signal.extend(data['sender'][col_name].values)
    
    if len(all_traj) == 0:
        return np.full(len(eval_points), np.nan)
    
    traj = np.array(all_traj)
    signal = np.array(all_signal)
    
    smoothed = np.zeros(len(eval_points))
    for i, ep in enumerate(eval_points):
        weights = np.exp(-0.5 * ((traj - ep) / sigma) ** 2)
        if weights.sum() > 0:
            smoothed[i] = np.sum(weights * signal) / np.sum(weights)
    
    return smoothed


def make_zone_cmap(zone_color):
    """Create colormap from white to zone color."""
    return LinearSegmentedColormap.from_list('zone_cmap', ['white', zone_color])


def draw_gradient_curve(ax, x1, y1, x2, y2, color1, color2, linewidth, alpha=0.8, n_segments=50):
    """Draw a curved line with gradient color from source to target."""
    center = (0.5, 0.5)
    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    ctrl_x = center[0] + (mid_x - center[0]) * 0.3
    ctrl_y = center[1] + (mid_y - center[1]) * 0.3
    
    # Generate Bezier curve points
    t = np.linspace(0, 1, n_segments)
    x_points = (1-t)**2 * x1 + 2*(1-t)*t * ctrl_x + t**2 * x2
    y_points = (1-t)**2 * y1 + 2*(1-t)*t * ctrl_y + t**2 * y2
    
    # Convert colors to RGB
    c1 = np.array(mcolors.to_rgb(color1))
    c2 = np.array(mcolors.to_rgb(color2))
    
    # Draw segments with interpolated colors
    for i in range(len(t) - 1):
        frac = i / (len(t) - 1)
        color = c1 * (1 - frac) + c2 * frac
        ax.plot([x_points[i], x_points[i+1]],
                [y_points[i], y_points[i+1]],
                color=color, linewidth=linewidth, alpha=alpha,
                solid_capstyle='round', zorder=1)


def draw_chord_diagram(ax, zone_df, zone_name, zone_color, global_max):
    """Draw chord diagram with gradient colored edges."""
    
    edge_counts = defaultdict(int)
    for _, row in zone_df.iterrows():
        src, tgt = row['Source_1'], row['Target_1']
        if pd.notna(src) and pd.notna(tgt):
            edge_counts[(src, tgt)] += 1
    
    if not edge_counts:
        ax.text(0.5, 0.5, 'No interactions', ha='center', va='center', fontsize=14)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
        return
    
    all_cts = set()
    ct_counts = defaultdict(int)
    for (src, tgt), count in edge_counts.items():
        all_cts.add(src); all_cts.add(tgt)
        ct_counts[src] += count; ct_counts[tgt] += count
    
    all_cts = sorted(all_cts, key=lambda x: ct_counts[x], reverse=True)
    n_cts = len(all_cts)
    
    radius = 0.38
    center = (0.5, 0.5)
    
    angles = np.linspace(90, 90 - 360, n_cts, endpoint=False)
    positions = {}
    for i, ct in enumerate(all_cts):
        angle_rad = np.radians(angles[i])
        x = center[0] + radius * np.cos(angle_rad)
        y = center[1] + radius * np.sin(angle_rad)
        positions[ct] = (x, y, angles[i])
    
    # Background circle
    bg_circle = Circle(center, radius + 0.08, facecolor='#F8F9FA', edgecolor='none', zorder=0)
    ax.add_patch(bg_circle)
    
    # Sort edges by count (draw smaller first)
    sorted_edges = sorted(edge_counts.items(), key=lambda x: x[1])
    
    for (src, tgt), count in sorted_edges:
        x1, y1, a1 = positions[src]
        x2, y2, a2 = positions[tgt]
        
        lw = 3 + 10 * (count / global_max)
        alpha = 0.5 + 0.4 * (count / global_max)
        
        src_color = CELLTYPE_COLORS.get(src, '#888888')
        tgt_color = CELLTYPE_COLORS.get(tgt, '#888888')
        
        if src == tgt:
            # Self-loop
            angle_rad = np.radians(a1)
            loop_r = 0.07
            cx = x1 + loop_r * 1.5 * np.cos(angle_rad)
            cy = y1 + loop_r * 1.5 * np.sin(angle_rad)
            
            loop = Circle((cx, cy), loop_r, fill=False,
                          color=src_color, linewidth=lw, alpha=alpha, zorder=1)
            ax.add_patch(loop)
            
            ax.text(cx + loop_r * np.cos(angle_rad) * 1.8,
                   cy + loop_r * np.sin(angle_rad) * 1.8,
                   str(count), fontsize=10, ha='center', va='center',
                   fontweight='bold', color=src_color, zorder=5,
                   path_effects=[pe.withStroke(linewidth=2, foreground='white')])
        else:
            draw_gradient_curve(ax, x1, y1, x2, y2, src_color, tgt_color,
                               linewidth=lw, alpha=alpha)
    
    # Draw nodes
    node_radius = 0.06
    
    for ct in all_cts:
        x, y, angle = positions[ct]
        color = CELLTYPE_COLORS.get(ct, '#888888')
        
        involvement = ct_counts[ct]
        size_mult = 0.85 + 0.35 * (involvement / max(ct_counts.values()))
        
        # Glow
        glow = Circle((x, y), node_radius * size_mult * 1.4,
                      facecolor=color, alpha=0.25, edgecolor='none', zorder=2)
        ax.add_patch(glow)
        
        # Node
        node = Circle((x, y), node_radius * size_mult,
                     facecolor=color, edgecolor='white', linewidth=2.5, zorder=3)
        ax.add_patch(node)
        
        # Label
        label_r = radius + 0.15
        angle_rad = np.radians(angle)
        lx = center[0] + label_r * np.cos(angle_rad)
        ly = center[1] + label_r * np.sin(angle_rad)
        
        if -90 < angle < 90:
            ha = 'left'; rotation = angle
        else:
            ha = 'right'; rotation = angle + 180
        
        label = CELLTYPE_SHORT.get(ct, ct)
        ax.text(lx, ly, label, fontsize=11, ha=ha, va='center',
               fontweight='bold', rotation=rotation, rotation_mode='anchor',
               path_effects=[pe.withStroke(linewidth=3, foreground='white')])
    
    # Center label
    ax.text(center[0], center[1], zone_name, fontsize=18, fontweight='bold',
           ha='center', va='center', color=zone_color,
           path_effects=[pe.withStroke(linewidth=5, foreground='white')])
    
    ax.text(center[0], center[1] - 0.09, f'{len(zone_df)} pairs',
           fontsize=12, ha='center', va='center', color='#666666')
    
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_aspect('equal'); ax.axis('off')


def generate_heatmaps(results_df, lr_signals, output_dir, formats, n_eval_points):
    """Generate zone-specific heatmaps."""
    
    zone_cmaps = {
        'Tumor': make_zone_cmap(ZONE_COLORS['Tumor']),
        'Interface': make_zone_cmap(ZONE_COLORS['Interface']),
        'Lymphoid': make_zone_cmap(ZONE_COLORS['Lymphoid']),
    }
    
    for zone in ['Tumor', 'Interface', 'Lymphoid']:
        zone_results = results_df[results_df['Zone'] == zone].copy()
        
        if len(zone_results) == 0:
            continue
        
        zone_results = zone_results.sort_values('LR_Peak')
        lr_pairs = zone_results['LR_Pair'].tolist()
        n_pairs = len(lr_pairs)
        
        # Build signal matrix
        signal_matrix = np.zeros((n_pairs, n_eval_points))
        for i, lr_pair in enumerate(lr_pairs):
            if lr_pair in lr_signals:
                signal_matrix[i, :] = lr_signals[lr_pair]
        
        # Row-normalize
        for i in range(n_pairs):
            row = signal_matrix[i, :]
            row_min, row_max = np.nanmin(row), np.nanmax(row)
            if row_max > row_min:
                signal_matrix[i, :] = (row - row_min) / (row_max - row_min)
            else:
                signal_matrix[i, :] = 0
        
        # Create figure
        fig_height = max(6, n_pairs * 0.38)
        fig, ax = plt.subplots(figsize=(14, fig_height))
        
        cmap = zone_cmaps[zone]
        
        im = ax.imshow(signal_matrix, aspect='auto', cmap=cmap,
                       extent=[0, 1, n_pairs, 0], vmin=0, vmax=1)
        
        # Zone boundaries
        for z, (z_start, z_end) in ZONE_BOUNDS.items():
            ax.axvline(x=z_start, color='gray', linestyle='--', linewidth=1.5, alpha=0.6)
            ax.axvline(x=z_end, color='gray', linestyle='--', linewidth=1.5, alpha=0.6)
        
        # Y-axis labels
        y_labels_left = []
        y_labels_right = []
        
        for _, row in zone_results.iterrows():
            y_labels_left.append(f"{row['LR_Pair']} (pk={row['LR_Peak']:.2f})")
            src = abbrev_cell(row['Source_1'])
            tgt = abbrev_cell(row['Target_1'])
            y_labels_right.append(f"{src}→{tgt}")
        
        ax.set_yticks(np.arange(n_pairs) + 0.5)
        ax.set_yticklabels(y_labels_left, fontsize=9)
        
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(np.arange(n_pairs) + 0.5)
        ax2.set_yticklabels(y_labels_right, fontsize=9, color=ZONE_COLORS[zone])
        
        ax.set_xlabel('Trajectory (Tumor → Interface → Lymphoid)', fontsize=12)
        ax.set_xticks([0.18, 0.54, 0.86])
        ax.set_xticklabels(['Tumor', 'Interface', 'Lymphoid'], fontsize=11)
        
        ax.set_title(f'{zone} Zone L-R Pairs (n={n_pairs})',
                     fontsize=14, fontweight='bold', color=ZONE_COLORS[zone])
        
        cbar = plt.colorbar(im, ax=ax, shrink=0.5, pad=0.15)
        cbar.set_label('Normalized Signal', fontsize=10)
        
        plt.tight_layout()
        
        for fmt in formats:
            plt.savefig(os.path.join(output_dir, f'Heatmap_{zone}.{fmt}'),
                        dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"   Saved: Heatmap_{zone}.{'/'.join(formats)}")


def generate_chord_network(results_df, output_dir, formats):
    """Generate three-panel chord network diagram."""
    
    # Calculate global max
    global_max_count = 0
    for zone in ['Tumor', 'Interface', 'Lymphoid']:
        zone_res = results_df[results_df['Zone'] == zone]
        edge_counts = defaultdict(int)
        for _, row in zone_res.iterrows():
            src, tgt = row['Source_1'], row['Target_1']
            if pd.notna(src) and pd.notna(tgt):
                edge_counts[(src, tgt)] += 1
        if edge_counts:
            global_max_count = max(global_max_count, max(edge_counts.values()))
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(21, 8))
    
    for ax_idx, zone in enumerate(['Tumor', 'Interface', 'Lymphoid']):
        zone_res = results_df[results_df['Zone'] == zone]
        draw_chord_diagram(axes[ax_idx], zone_res, zone, ZONE_COLORS[zone], global_max_count)
    
    fig.suptitle('Cell-Cell Communication Networks by Spatial Zone',
                 fontsize=20, fontweight='bold', y=1.02)
    
    # Legend
    ct_elements = [mpatches.Patch(facecolor=color, edgecolor='black', linewidth=0.5,
                                  label=CELLTYPE_SHORT.get(ct, ct))
                   for ct, color in CELLTYPE_COLORS.items()
                   if ct in results_df['Source_1'].values or ct in results_df['Target_1'].values]
    
    lw_legend_counts = [c for c in [1, 3, 6, 10] if c <= global_max_count]
    if global_max_count not in lw_legend_counts:
        lw_legend_counts.append(global_max_count)
    
    lw_elements = [Line2D([0], [0], color='#666666', linewidth=3 + 10*(c/global_max_count),
                          label=f'{c} L-R pairs') for c in sorted(lw_legend_counts)]
    
    all_elements = ct_elements + [mpatches.Patch(facecolor='none', edgecolor='none', label='')] + lw_elements
    
    fig.legend(handles=all_elements, loc='lower center', bbox_to_anchor=(0.5, -0.06),
              ncol=7, fontsize=10, frameon=True, fancybox=True,
              title='Cell Types                                                    Line Width = # L-R Pairs',
              title_fontsize=11)
    
    plt.tight_layout()
    
    for fmt in formats:
        plt.savefig(os.path.join(output_dir, f'ChordNetwork.{fmt}'),
                    bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"   Saved: ChordNetwork.{'/'.join(formats)}")


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*70)
    print("CLIP STEP 4: VISUALIZATION")
    print("="*70)
    
    # Load CLIP results
    print("\n1. Loading CLIP results...", flush=True)
    results_df = pd.read_csv(args.clip_results)
    print(f"   Loaded {len(results_df)} L-R pairs")
    
    # Load COMMOT data for heatmaps
    print("\n2. Loading COMMOT data...", flush=True)
    eval_points = np.linspace(0, 1, args.n_eval_points)
    
    sample_data = {}
    for sample in args.samples:
        h5ad_path = os.path.join(args.commot_dir, f'{sample}_commot_results.h5ad')
        if os.path.exists(h5ad_path):
            adata = sc.read_h5ad(h5ad_path)
            sample_data[sample] = {
                'trajectory': adata.obs['trajectory_score'].values,
                'sender': adata.obsm['commot-cellchat-sum-sender']
            }
            print(f"   {sample}...")
    
    # Compute L-R signals
    print("\n3. Computing rolling Gaussian signals...", flush=True)
    lr_signals = {}
    for lr_pair in results_df['LR_Pair'].unique():
        lr_signals[lr_pair] = compute_rolling_gaussian(
            lr_pair, sample_data, eval_points, args.sigma
        )
    
    # Generate heatmaps
    print("\n4. Generating heatmaps...", flush=True)
    generate_heatmaps(results_df, lr_signals, args.output_dir, args.formats, args.n_eval_points)
    
    # Generate chord network
    print("\n5. Generating chord network...", flush=True)
    generate_chord_network(results_df, args.output_dir, args.formats)
    
    print("\n" + "="*70)
    print("STEP 4 COMPLETE!")
    print(f"Figures saved to: {args.output_dir}")
    print("="*70)


if __name__ == '__main__':
    main()
