#!/usr/bin/env python3
"""
================================================================================
CLIP Pipeline - Step 1: COMMOT Ligand-Receptor Analysis
================================================================================

Performs spatially-constrained ligand-receptor communication inference using 
COMMOT (COMMunication analysis by Optimal Transport) with CellChatDB.

Input:
    - HD Visium h5ad files with spatial coordinates
    - Trajectory scores (ONTraC-derived)
    
Output:
    - Per-sample h5ad files with COMMOT sender/receiver scores
    - LR_scores_by_zone.csv

Reference:
    Cang Z, et al. Screening cell-cell communication in spatial transcriptomics 
    via collective optimal transport. Nat Methods 20, 218-228 (2023).
================================================================================
"""

import pandas as pd
import numpy as np
import scanpy as sc
import commot as ct
import os
import argparse
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Default parameters
DEFAULT_DIS_THR = 50        # Distance threshold in µm (~6 spots at 8µm resolution)
DEFAULT_N_SUBSAMPLE = 30000  # Subsample size for computational efficiency

def parse_args():
    parser = argparse.ArgumentParser(
        description='CLIP Step 1: COMMOT Ligand-Receptor Analysis'
    )
    parser.add_argument('--input_dir', type=str, required=True,
                        help='Directory containing h5ad files')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for COMMOT results')
    parser.add_argument('--trajectory_file', type=str, required=True,
                        help='CSV file with trajectory scores')
    parser.add_argument('--samples', nargs='+', required=True,
                        help='Sample names to process')
    parser.add_argument('--dis_thr', type=int, default=DEFAULT_DIS_THR,
                        help=f'Spatial distance threshold in µm (default: {DEFAULT_DIS_THR})')
    parser.add_argument('--n_subsample', type=int, default=DEFAULT_N_SUBSAMPLE,
                        help=f'Number of spots to subsample (default: {DEFAULT_N_SUBSAMPLE})')
    parser.add_argument('--h5ad_pattern', type=str, default='{sample}_8um.h5ad',
                        help='Pattern for h5ad filenames (default: {sample}_8um.h5ad)')
    return parser.parse_args()


def load_trajectory_scores(trajectory_file):
    """Load trajectory scores from CSV file."""
    print(f"Loading trajectory scores from {trajectory_file}...")
    gradient_df = pd.read_csv(trajectory_file)
    
    # Compute combined score if not present
    if 'combined_score' not in gradient_df.columns:
        if 'NT_oriented' in gradient_df.columns and 'spatial_gradient' in gradient_df.columns:
            gradient_df['combined_score'] = (
                gradient_df['NT_oriented'] + gradient_df['spatial_gradient']
            ) / 2
        elif 'trajectory_score' in gradient_df.columns:
            gradient_df['combined_score'] = gradient_df['trajectory_score']
        else:
            raise ValueError("Could not find trajectory score columns")
    
    return gradient_df


def get_cellchat_database():
    """Load CellChat ligand-receptor database for human."""
    print("Loading CellChat ligand-receptor database...")
    df_ligrec = ct.pp.ligand_receptor_database(database='CellChat', species='human')
    print(f"  Total L-R pairs in database: {len(df_ligrec)}")
    return df_ligrec


def filter_lr_pairs(df_ligrec, genes_present):
    """Filter L-R pairs to those with genes present in the data."""
    
    def check_genes_present(row):
        ligand = row[0]
        receptor = row[1]
        # Handle heteromeric receptors (e.g., "TGFBR1_TGFBR2")
        ligand_genes = ligand.split('_')
        receptor_genes = receptor.split('_')
        return (all(g in genes_present for g in ligand_genes) and 
                all(g in genes_present for g in receptor_genes))
    
    df_filtered = df_ligrec[[check_genes_present(row) for _, row in df_ligrec.iterrows()]]
    return df_filtered


def subsample_stratified(sample_gradient_valid, n_subsample, n_bins=10):
    """Stratified subsampling across trajectory bins."""
    sample_gradient_valid['traj_bin'] = pd.cut(
        sample_gradient_valid['combined_score'], bins=n_bins
    )
    sampled = sample_gradient_valid.groupby('traj_bin', group_keys=False).apply(
        lambda x: x.sample(min(len(x), n_subsample // n_bins), random_state=42)
    )
    return sampled['original_spot_id'].tolist()


def run_commot_single_sample(
    sample, 
    input_dir, 
    output_dir, 
    gradient_df, 
    df_ligrec,
    h5ad_pattern,
    dis_thr=DEFAULT_DIS_THR,
    n_subsample=DEFAULT_N_SUBSAMPLE
):
    """Run COMMOT analysis on a single sample."""
    
    print(f"\n{'='*60}")
    print(f"Processing {sample}...")
    print(f"{'='*60}")
    
    # Load h5ad
    h5ad_filename = h5ad_pattern.format(sample=sample)
    h5ad_path = os.path.join(input_dir, h5ad_filename)
    
    if not os.path.exists(h5ad_path):
        print(f"  WARNING: File not found - {h5ad_path}")
        return None
    
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded: {adata.shape[0]} spots, {adata.shape[1]} genes")
    
    # Get trajectory scores for this sample
    sample_gradient = gradient_df[gradient_df['Sample'] == sample].copy()
    sample_gradient['original_spot_id'] = sample_gradient['Cell_ID'].str.replace(
        f'{sample}_', '', regex=False
    )
    
    # Filter to spots with trajectory scores
    valid_spots = [s for s in sample_gradient['original_spot_id'] if s in adata.obs.index]
    print(f"  Spots with trajectory scores: {len(valid_spots)}")
    
    if len(valid_spots) == 0:
        print(f"  WARNING: No valid spots found for {sample}")
        return None
    
    # Subsample for computational efficiency
    sample_gradient_valid = sample_gradient[
        sample_gradient['original_spot_id'].isin(valid_spots)
    ]
    
    if len(valid_spots) > n_subsample:
        print(f"  Subsampling to {n_subsample} spots (stratified by trajectory)...")
        valid_spots = subsample_stratified(sample_gradient_valid, n_subsample)
        print(f"  Subsampled spots: {len(valid_spots)}")
    
    # Subset adata
    adata_sub = adata[valid_spots].copy()
    
    # Add trajectory scores to adata
    score_lookup = dict(zip(
        sample_gradient['original_spot_id'], 
        sample_gradient['combined_score']
    ))
    adata_sub.obs['trajectory_score'] = [score_lookup[s] for s in adata_sub.obs.index]
    
    # Basic preprocessing for COMMOT
    print("  Preprocessing for COMMOT...")
    if adata_sub.X.max() > 100:
        sc.pp.normalize_total(adata_sub, target_sum=1e4)
        sc.pp.log1p(adata_sub)
    
    # Ensure spatial coordinates are present
    if 'spatial' not in adata_sub.obsm:
        if 'x' in adata_sub.obs.columns and 'y' in adata_sub.obs.columns:
            adata_sub.obsm['spatial'] = adata_sub.obs[['x', 'y']].values
        else:
            print(f"  WARNING: No spatial coordinates found for {sample}")
            return None
    
    # Filter L-R pairs to genes present in data
    genes_present = set(adata_sub.var_names)
    df_ligrec_filtered = filter_lr_pairs(df_ligrec, genes_present)
    print(f"  L-R pairs with genes present: {len(df_ligrec_filtered)}")
    
    if len(df_ligrec_filtered) == 0:
        print(f"  WARNING: No L-R pairs found for {sample}")
        return None
    
    # Run COMMOT
    print(f"  Running COMMOT (dis_thr={dis_thr}µm)...")
    try:
        ct.tl.spatial_communication(
            adata_sub,
            database_name='cellchat',
            df_ligrec=df_ligrec_filtered,
            dis_thr=dis_thr,
            heteromeric=True,
            pathway_sum=True
        )
        print("  COMMOT completed successfully!")
    except Exception as e:
        print(f"  ERROR in COMMOT: {e}")
        return None
    
    # Save sample-specific results
    sample_output = os.path.join(output_dir, f'{sample}_commot_results.h5ad')
    adata_sub.write(sample_output)
    print(f"  Saved: {sample_output}")
    
    return adata_sub


def extract_lr_scores(adata_sub, sample):
    """Extract L-R scores by zone from COMMOT results."""
    
    # Define zone boundaries
    TUMOR_END = 0.36
    INTERFACE_END = 0.725
    
    # Define zones
    adata_sub.obs['zone'] = pd.cut(
        adata_sub.obs['trajectory_score'],
        bins=[0, TUMOR_END, INTERFACE_END, 1.0],
        labels=['Tumor', 'Interface', 'Lymphoid']
    )
    
    lr_scores = []
    
    sender_key = 'commot-cellchat-sum-sender'
    receiver_key = 'commot-cellchat-sum-receiver'
    
    if sender_key in adata_sub.obsm and receiver_key in adata_sub.obsm:
        sender_df = adata_sub.obsm[sender_key]
        receiver_df = adata_sub.obsm[receiver_key]
        
        for col in sender_df.columns:
            if col.startswith('s-') and '-total' not in col:
                lr_name = col.replace('s-', '')
                
                s_signal = sender_df[col].values
                r_col = f'r-{lr_name}'
                r_signal = receiver_df[r_col].values if r_col in receiver_df.columns else np.zeros_like(s_signal)
                
                combined_signal = s_signal + r_signal
                
                for zone in ['Tumor', 'Interface', 'Lymphoid']:
                    zone_mask = adata_sub.obs['zone'] == zone
                    if zone_mask.sum() > 0:
                        lr_scores.append({
                            'Sample': sample,
                            'LR_pair': lr_name,
                            'Zone': zone,
                            'Mean_Signal': np.mean(combined_signal[zone_mask]),
                            'Median_Signal': np.median(combined_signal[zone_mask]),
                            'N_spots': zone_mask.sum(),
                            'Pct_nonzero': (combined_signal[zone_mask] > 0).mean() * 100
                        })
    
    return lr_scores


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*70)
    print("CLIP STEP 1: COMMOT LIGAND-RECEPTOR ANALYSIS")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  Input directory: {args.input_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Samples: {args.samples}")
    print(f"  Distance threshold: {args.dis_thr}µm")
    print(f"  Subsample size: {args.n_subsample}")
    
    # Load data
    gradient_df = load_trajectory_scores(args.trajectory_file)
    df_ligrec = get_cellchat_database()
    
    # Process each sample
    all_lr_scores = []
    
    for sample in args.samples:
        adata_sub = run_commot_single_sample(
            sample=sample,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            gradient_df=gradient_df,
            df_ligrec=df_ligrec,
            h5ad_pattern=args.h5ad_pattern,
            dis_thr=args.dis_thr,
            n_subsample=args.n_subsample
        )
        
        if adata_sub is not None:
            scores = extract_lr_scores(adata_sub, sample)
            all_lr_scores.extend(scores)
    
    # Save aggregated results
    if len(all_lr_scores) > 0:
        lr_df = pd.DataFrame(all_lr_scores)
        output_file = os.path.join(args.output_dir, 'LR_scores_by_zone.csv')
        lr_df.to_csv(output_file, index=False)
        print(f"\nSaved: {output_file}")
    else:
        print("\nWARNING: No L-R scores collected!")
    
    print("\n" + "="*70)
    print("STEP 1 COMPLETE!")
    print("="*70)


if __name__ == '__main__':
    main()
