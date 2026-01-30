#!/usr/bin/env python3
"""
================================================================================
CLIP Pipeline - Step 3: CLIP Scoring
================================================================================

CLIP: Cross-referenced Ligand-receptor Interaction Peak scoring

Assigns source and target cell types to each L-R pair by integrating:
    1. Expression probability from scRNA-seq
    2. Spatial peak matching along the tissue trajectory

CLIP Score = √(E_L × E_R) × exp(−d²/2σ²)

Where:
    E_L: Fraction of source cells expressing the ligand gene
    E_R: Fraction of target cells expressing the receptor gene
    d: Distance between L-R signal peak and cell-cell co-localization peak
    σ: Gaussian decay bandwidth (default: 0.08)

Input:
    - COMMOT results (h5ad files)
    - scRNA-seq expression percentages
    - Cell-cell interaction density data
    
Output:
    - CLIP_scores.csv: All L-R pairs with cell type assignments and scores

================================================================================
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from collections import defaultdict
import os
import argparse
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_N_EVAL_POINTS = 150    # Evaluation points along trajectory
DEFAULT_SIGMA = 0.05           # Rolling Gaussian bandwidth
DEFAULT_PEAK_SIGMA = 0.08      # Peak matching Gaussian decay
DEFAULT_MIN_EXPRESSION = 0.02  # Minimum expression threshold (2%)

# Zone boundaries
ZONE_BOUNDS = {
    'Tumor': (0.0, 0.36),
    'Interface': (0.36, 0.725),
    'Lymphoid': (0.725, 1.0)
}

# Cell type abbreviations for display
CELL_ABBREV = {
    'T.cells': 'T', 'Macrophages': 'Mac', 'Dendritic.cells': 'DC',
    'B_Plasma.cells': 'B', 'Fibroblasts': 'Fib', 'Malignant.cells': 'Mal',
    'Endothelial.cells': 'Endo', 'Epithelial.cells': 'Epi',
    'Mast.cells': 'Mast', 'Myocytes': 'Myo',
}

# Cell type name mapping (CARD → scRNA-seq)
CELL_MAPPING = {
    'T.cells': 'T cells', 'Macrophages': 'Macrophages',
    'Dendritic.cells': 'Dendritic cells', 'B_Plasma.cells': 'B_Plasma cells',
    'Fibroblasts': 'Fibroblasts', 'Malignant.cells': 'Malignant cells',
    'Endothelial.cells': 'Endothelial cells', 'Epithelial.cells': 'Epithelial cells',
    'Mast.cells': 'Mast cells', 'Myocytes': 'Myocytes',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description='CLIP Step 3: CLIP Scoring'
    )
    parser.add_argument('--commot_dir', type=str, required=True,
                        help='Directory containing COMMOT results')
    parser.add_argument('--scrna_file', type=str, required=True,
                        help='CSV file with scRNA-seq expression percentages')
    parser.add_argument('--interaction_file', type=str, required=True,
                        help='CSV file with cell-cell interaction densities')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--samples', nargs='+', required=True,
                        help='Sample names to process')
    parser.add_argument('--signal_stats_file', type=str, default=None,
                        help='Optional: CSV with signal statistics for filtering')
    parser.add_argument('--n_eval_points', type=int, default=DEFAULT_N_EVAL_POINTS,
                        help=f'Trajectory evaluation points (default: {DEFAULT_N_EVAL_POINTS})')
    parser.add_argument('--sigma', type=float, default=DEFAULT_SIGMA,
                        help=f'Rolling Gaussian bandwidth (default: {DEFAULT_SIGMA})')
    parser.add_argument('--peak_sigma', type=float, default=DEFAULT_PEAK_SIGMA,
                        help=f'Peak matching decay sigma (default: {DEFAULT_PEAK_SIGMA})')
    parser.add_argument('--min_expression', type=float, default=DEFAULT_MIN_EXPRESSION,
                        help=f'Minimum expression threshold (default: {DEFAULT_MIN_EXPRESSION})')
    return parser.parse_args()


def abbrev_cell(cell):
    """Get abbreviated cell type name."""
    if pd.isna(cell):
        return '?'
    return CELL_ABBREV.get(cell, cell[:4])


def get_zone_from_peak(peak_loc):
    """Assign zone based on peak location."""
    for zone, (z_start, z_end) in ZONE_BOUNDS.items():
        if z_start <= peak_loc < z_end:
            return zone
    return 'Lymphoid'


def parse_lr_pair(lr_pair):
    """Parse L-R pair string into ligand and receptor components."""
    parts = lr_pair.split('-')
    if len(parts) < 2:
        return None, []
    ligand = parts[0]
    receptors = '-'.join(parts[1:]).split('_')
    return ligand, receptors


def compute_rolling_gaussian(lr_pair, sample_data, eval_points, sigma):
    """
    Compute rolling Gaussian-smoothed signal along trajectory.
    
    Returns smoothed signal array and peak position.
    """
    all_traj, all_signal = [], []
    
    for sample, data in sample_data.items():
        col_name = f's-{lr_pair}'
        if col_name not in data['sender'].columns:
            continue
        all_traj.extend(data['trajectory'])
        all_signal.extend(data['sender'][col_name].values)
    
    if len(all_traj) == 0:
        return np.full(len(eval_points), np.nan), np.nan
    
    traj = np.array(all_traj)
    signal = np.array(all_signal)
    
    smoothed = np.zeros(len(eval_points))
    for i, ep in enumerate(eval_points):
        weights = np.exp(-0.5 * ((traj - ep) / sigma) ** 2)
        if weights.sum() > 0:
            smoothed[i] = np.sum(weights * signal) / np.sum(weights)
    
    # Find peak
    if np.nanmax(smoothed) > 0:
        peak_loc = eval_points[np.nanargmax(smoothed)]
    else:
        peak_loc = np.nan
    
    return smoothed, peak_loc


def compute_cellcell_peaks(cellcell_df, eval_points):
    """Compute peak positions for all cell-cell pairs."""
    trajectory_bins = cellcell_df.index.values.astype(float)
    pair_cols = [c for c in cellcell_df.columns if '|' in c and not c.startswith('count_')]
    
    cellcell_peaks = {}
    cellcell_interp = {}
    
    for pair in pair_cols:
        density = cellcell_df[pair].values
        interp = np.interp(eval_points, trajectory_bins, density)
        cellcell_interp[pair] = interp
        
        if np.nanmax(interp) > 0:
            cellcell_peaks[pair] = eval_points[np.nanargmax(interp)]
        else:
            cellcell_peaks[pair] = np.nan
    
    return cellcell_peaks, cellcell_interp


def clip_score(ligand_expr, receptor_expr, lr_peak, cellcell_peak, peak_sigma):
    """
    Calculate CLIP score for a source→target cell type pair.
    
    CLIP Score = √(E_L × E_R) × exp(−d²/2σ²)
    
    Parameters:
        ligand_expr: Fraction of source cells expressing ligand (0-1)
        receptor_expr: Fraction of target cells expressing receptor (0-1)
        lr_peak: Trajectory position of L-R signal peak
        cellcell_peak: Trajectory position of cell-cell co-localization peak
        peak_sigma: Gaussian decay bandwidth
    
    Returns:
        Combined CLIP score
    """
    # Expression score (geometric mean)
    expr_score = np.sqrt(ligand_expr * receptor_expr)
    
    # Peak matching score (Gaussian decay)
    if np.isnan(cellcell_peak):
        peak_match = 0.1  # Low default for unknown co-localization
    else:
        peak_distance = abs(cellcell_peak - lr_peak)
        peak_match = np.exp(-0.5 * (peak_distance / peak_sigma) ** 2)
    
    return expr_score * peak_match


def score_lr_pair(lr_pair, lr_peak, pct_df, cellcell_peaks, peak_sigma, min_expression):
    """
    Score all candidate source→target pairs for a given L-R pair.
    
    Returns list of scored pairs sorted by CLIP score (descending).
    """
    ligand, receptors = parse_lr_pair(lr_pair)
    if ligand is None:
        return []
    
    all_celltypes = list(pct_df.columns)
    
    # Get expression for all cell types
    ligand_expr = {}
    for ct in all_celltypes:
        if ligand in pct_df.index:
            ligand_expr[ct] = pct_df.loc[ligand, ct] / 100  # Convert to 0-1
        else:
            ligand_expr[ct] = 0
    
    receptor_expr = {}
    for ct in all_celltypes:
        max_pct = 0
        for r in receptors:
            if r in pct_df.index:
                max_pct = max(max_pct, pct_df.loc[r, ct])
        receptor_expr[ct] = max_pct / 100
    
    # Score all cell-cell pairs
    pair_scores = []
    
    for source_ct in all_celltypes:
        for target_ct in all_celltypes:
            expr_source = ligand_expr.get(source_ct, 0)
            expr_target = receptor_expr.get(target_ct, 0)
            
            # Skip if either doesn't express (< threshold)
            if expr_source < min_expression or expr_target < min_expression:
                continue
            
            # Get cell-cell co-localization peak
            source_mapped = CELL_MAPPING.get(source_ct, source_ct)
            target_mapped = CELL_MAPPING.get(target_ct, target_ct)
            cellcell_name = f"{source_mapped}|{target_mapped}"
            
            cellcell_peak = cellcell_peaks.get(cellcell_name, np.nan)
            
            if np.isnan(cellcell_peak):
                # Try reverse direction
                cellcell_name_rev = f"{target_mapped}|{source_mapped}"
                cellcell_peak = cellcell_peaks.get(cellcell_name_rev, np.nan)
            
            # Calculate CLIP score
            combined_score = clip_score(
                expr_source, expr_target, lr_peak, cellcell_peak, peak_sigma
            )
            
            peak_distance = abs(cellcell_peak - lr_peak) if not np.isnan(cellcell_peak) else np.nan
            peak_match = np.exp(-0.5 * (peak_distance / peak_sigma) ** 2) if not np.isnan(peak_distance) else 0.1
            
            pair_scores.append({
                'Source': source_ct,
                'Target': target_ct,
                'Expr_Source': expr_source,
                'Expr_Target': expr_target,
                'Expr_Score': np.sqrt(expr_source * expr_target),
                'CellCell_Peak': cellcell_peak,
                'Peak_Distance': peak_distance,
                'Peak_Match': peak_match,
                'CLIP_Score': combined_score,
            })
    
    # Sort by CLIP score
    pair_scores = sorted(pair_scores, key=lambda x: -x['CLIP_Score'])
    
    return pair_scores


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*70)
    print("CLIP STEP 3: CLIP SCORING")
    print("Cross-referenced Ligand-receptor Interaction Peak")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  COMMOT directory: {args.commot_dir}")
    print(f"  scRNA-seq file: {args.scrna_file}")
    print(f"  Interaction file: {args.interaction_file}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Evaluation points: {args.n_eval_points}")
    print(f"  Rolling Gaussian σ: {args.sigma}")
    print(f"  Peak matching σ: {args.peak_sigma}")
    print(f"  Min expression: {args.min_expression*100}%")
    
    # Load data
    print("\n1. Loading data...", flush=True)
    
    # scRNA-seq reference
    pct_df = pd.read_csv(args.scrna_file, index_col=0)
    print(f"   scRNA-seq: {len(pct_df)} genes × {len(pct_df.columns)} cell types")
    
    # Cell-cell interaction data
    cellcell_df = pd.read_csv(args.interaction_file, index_col=0)
    print(f"   Cell-cell interactions loaded")
    
    # Optional signal filtering
    valid_lr_pairs = None
    if args.signal_stats_file and os.path.exists(args.signal_stats_file):
        signal_stats = pd.read_csv(args.signal_stats_file)
        valid_lr_pairs = signal_stats[
            (signal_stats['Total_Signal'] >= 50) &
            (signal_stats['N_nonzero'] >= 20)
        ]['LR_Pair'].tolist()
        print(f"   Valid L-R pairs (from signal filter): {len(valid_lr_pairs)}")
    
    # Load COMMOT data
    print("\n2. Loading COMMOT data...", flush=True)
    
    eval_points = np.linspace(0, 1, args.n_eval_points)
    
    sample_data = {}
    for sample in args.samples:
        h5ad_path = os.path.join(args.commot_dir, f'{sample}_commot_results.h5ad')
        if not os.path.exists(h5ad_path):
            print(f"   {sample}: NOT FOUND")
            continue
        print(f"   {sample}...", flush=True)
        adata = sc.read_h5ad(h5ad_path)
        sample_data[sample] = {
            'trajectory': adata.obs['trajectory_score'].values,
            'sender': adata.obsm['commot-cellchat-sum-sender']
        }
    
    # Get all L-R pairs from COMMOT results
    all_lr_pairs = set()
    for sample, data in sample_data.items():
        for col in data['sender'].columns:
            if col.startswith('s-') and '-total' not in col:
                all_lr_pairs.add(col.replace('s-', ''))
    
    if valid_lr_pairs:
        all_lr_pairs = all_lr_pairs.intersection(set(valid_lr_pairs))
    
    print(f"   Total L-R pairs to score: {len(all_lr_pairs)}")
    
    # Compute L-R signal peaks
    print("\n3. Computing L-R signal peaks...", flush=True)
    
    lr_signals = {}
    lr_peaks = {}
    
    for lr_pair in all_lr_pairs:
        if lr_pair == 'total-total':
            continue
        signal, peak = compute_rolling_gaussian(lr_pair, sample_data, eval_points, args.sigma)
        if not np.isnan(peak):
            lr_signals[lr_pair] = signal
            lr_peaks[lr_pair] = peak
    
    print(f"   L-R pairs with valid peaks: {len(lr_peaks)}")
    
    # Compute cell-cell co-localization peaks
    print("\n4. Computing cell-cell co-localization peaks...", flush=True)
    
    cellcell_peaks, _ = compute_cellcell_peaks(cellcell_df, eval_points)
    valid_peaks = len([p for p in cellcell_peaks.values() if not np.isnan(p)])
    print(f"   Cell-cell pairs with peaks: {valid_peaks}")
    
    # Run CLIP scoring
    print("\n5. Running CLIP scoring...", flush=True)
    
    results = []
    
    for i, lr_pair in enumerate(lr_peaks.keys()):
        if (i + 1) % 100 == 0:
            print(f"   Processed {i + 1}/{len(lr_peaks)} L-R pairs...", flush=True)
        
        lr_peak = lr_peaks[lr_pair]
        zone = get_zone_from_peak(lr_peak)
        
        ligand, receptors = parse_lr_pair(lr_pair)
        if ligand is None:
            continue
        
        # Score all candidate pairs
        pair_scores = score_lr_pair(
            lr_pair, lr_peak, pct_df, cellcell_peaks, args.peak_sigma, args.min_expression
        )
        
        top3 = pair_scores[:3] if pair_scores else []
        
        if top3:
            best = top3[0]
            results.append({
                'LR_Pair': lr_pair,
                'Zone': zone,
                'LR_Peak': lr_peak,
                'Ligand': ligand,
                'Receptors': '/'.join(receptors),
                # Best match
                'Source_1': best['Source'],
                'Target_1': best['Target'],
                'Expr_Source_1': best['Expr_Source'],
                'Expr_Target_1': best['Expr_Target'],
                'Expr_Score_1': best['Expr_Score'],
                'Peak_Distance_1': best['Peak_Distance'],
                'Peak_Match_1': best['Peak_Match'],
                'CLIP_Score_1': best['CLIP_Score'],
                # 2nd best
                'Source_2': top3[1]['Source'] if len(top3) > 1 else None,
                'Target_2': top3[1]['Target'] if len(top3) > 1 else None,
                'CLIP_Score_2': top3[1]['CLIP_Score'] if len(top3) > 1 else None,
                # 3rd best
                'Source_3': top3[2]['Source'] if len(top3) > 2 else None,
                'Target_3': top3[2]['Target'] if len(top3) > 2 else None,
                'CLIP_Score_3': top3[2]['CLIP_Score'] if len(top3) > 2 else None,
            })
    
    results_df = pd.DataFrame(results)
    
    # Save results
    output_file = os.path.join(args.output_dir, 'CLIP_scores.csv')
    results_df.to_csv(output_file, index=False)
    print(f"\nSaved: {output_file} ({len(results_df)} pairs)")
    
    # Print summary by zone
    print("\n" + "="*70)
    print("RESULTS BY ZONE")
    print("="*70)
    
    print(f"\nZone distribution:")
    print(results_df['Zone'].value_counts())
    
    for zone in ['Tumor', 'Interface', 'Lymphoid']:
        zone_df = results_df[results_df['Zone'] == zone].sort_values('CLIP_Score_1', ascending=False)
        
        print(f"\n{zone.upper()} ({len(zone_df)} pairs)")
        print("-"*100)
        print(f"{'LR_Pair':<22} {'Peak':>6} {'Source→Target':<18} {'Expr_S':>7} {'Expr_T':>7} {'PkDist':>7} {'CLIP':>10}")
        print("-"*100)
        
        for _, row in zone_df.head(15).iterrows():
            src = abbrev_cell(row['Source_1'])
            tgt = abbrev_cell(row['Target_1'])
            pair = f"{src}→{tgt}"
            pk_dist = f"{row['Peak_Distance_1']:.3f}" if not pd.isna(row['Peak_Distance_1']) else "N/A"
            print(f"{row['LR_Pair']:<22} {row['LR_Peak']:>6.3f} {pair:<18} {row['Expr_Source_1']:>7.3f} {row['Expr_Target_1']:>7.3f} {pk_dist:>7} {row['CLIP_Score_1']:>10.4f}")
    
    # Verification: CXCL10-CXCR3
    print("\n" + "="*70)
    print("VERIFICATION: CXCL10-CXCR3")
    print("="*70)
    
    cxcl10 = results_df[results_df['LR_Pair'] == 'CXCL10-CXCR3']
    if len(cxcl10) > 0:
        row = cxcl10.iloc[0]
        print(f"\nL-R Peak: {row['LR_Peak']:.3f} → Zone: {row['Zone']}")
        print(f"\nBest match: {row['Source_1']} → {row['Target_1']}")
        print(f"  Ligand expression (CXCL10): {row['Expr_Source_1']*100:.1f}%")
        print(f"  Receptor expression (CXCR3): {row['Expr_Target_1']*100:.1f}%")
        print(f"  Expression score: {row['Expr_Score_1']:.3f}")
        print(f"  Peak distance: {row['Peak_Distance_1']:.3f}")
        print(f"  Peak match: {row['Peak_Match_1']:.3f}")
        print(f"  CLIP Score: {row['CLIP_Score_1']:.4f}")
        
        if row['Source_2']:
            print(f"\n2nd best: {row['Source_2']} → {row['Target_2']} (CLIP: {row['CLIP_Score_2']:.4f})")
        if row['Source_3']:
            print(f"3rd best: {row['Source_3']} → {row['Target_3']} (CLIP: {row['CLIP_Score_3']:.4f})")
    else:
        print("CXCL10-CXCR3 not found in results")
    
    print("\n" + "="*70)
    print("STEP 3 COMPLETE!")
    print("="*70)


if __name__ == '__main__':
    main()
