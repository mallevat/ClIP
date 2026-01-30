#!/usr/bin/env python3
"""
================================================================================
CLIP Pipeline - Step 2: Zone Mapping
================================================================================

Maps ligand-receptor signals across the spatial trajectory and assigns each
L-R pair to a spatial zone (Tumor, Interface, or Lymphoid).

Method: "High Interface Signal" approach
    - Interface: Any L-R pair with meaningful Interface signal (≥ threshold)
    - Tumor: Peak in Tumor AND low Interface signal
    - Lymphoid: Peak in Lymphoid AND low Interface signal

Input:
    - COMMOT results (h5ad files from Step 1)
    
Output:
    - LR_all_zones.csv: All L-R pairs with zone assignments
    - LR_{zone}.csv: Zone-specific L-R pair lists
    - summary_stats.csv: Summary statistics

================================================================================
"""

import pandas as pd
import numpy as np
import scanpy as sc
from collections import defaultdict
import os
import argparse
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

# Default zone boundaries (based on ONTraC niche centers)
# Tumor center: 0.07, Interface center: 0.65, Lymphoid center: 0.80
DEFAULT_TUMOR_END = 0.36        # Midpoint between 0.07 and 0.65
DEFAULT_INTERFACE_END = 0.725   # Midpoint between 0.65 and 0.80

# Signal thresholds
DEFAULT_MIN_TOTAL_SIGNAL = 0.005
DEFAULT_MIN_INTERFACE_SIGNAL = 0.005


def parse_args():
    parser = argparse.ArgumentParser(
        description='CLIP Step 2: Zone Mapping'
    )
    parser.add_argument('--commot_dir', type=str, required=True,
                        help='Directory containing COMMOT results')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory')
    parser.add_argument('--samples', nargs='+', required=True,
                        help='Sample names to process')
    parser.add_argument('--tumor_end', type=float, default=DEFAULT_TUMOR_END,
                        help=f'Tumor/Interface boundary (default: {DEFAULT_TUMOR_END})')
    parser.add_argument('--interface_end', type=float, default=DEFAULT_INTERFACE_END,
                        help=f'Interface/Lymphoid boundary (default: {DEFAULT_INTERFACE_END})')
    parser.add_argument('--min_total_signal', type=float, default=DEFAULT_MIN_TOTAL_SIGNAL,
                        help=f'Minimum total signal threshold (default: {DEFAULT_MIN_TOTAL_SIGNAL})')
    parser.add_argument('--min_interface_signal', type=float, default=DEFAULT_MIN_INTERFACE_SIGNAL,
                        help=f'Minimum Interface signal threshold (default: {DEFAULT_MIN_INTERFACE_SIGNAL})')
    return parser.parse_args()


def collect_lr_signals(commot_dir, samples, tumor_end, interface_end):
    """Collect L-R signals from all samples across zones."""
    
    lr_zone_signals = defaultdict(lambda: {'Tumor': [], 'Interface': [], 'Lymphoid': []})
    
    for sample in samples:
        h5ad_path = os.path.join(commot_dir, f'{sample}_commot_results.h5ad')
        
        if not os.path.exists(h5ad_path):
            print(f"  {sample}: NOT FOUND - {h5ad_path}")
            continue
        
        print(f"  Loading {sample}...", flush=True)
        adata = sc.read_h5ad(h5ad_path)
        
        # Get trajectory scores
        trajectory = adata.obs['trajectory_score'].values
        
        # Define zone masks
        tumor_mask = trajectory < tumor_end
        interface_mask = (trajectory >= tumor_end) & (trajectory < interface_end)
        lymphoid_mask = trajectory >= interface_end
        
        print(f"    Spots: Tumor={tumor_mask.sum():,}, Interface={interface_mask.sum():,}, Lymphoid={lymphoid_mask.sum():,}")
        
        # Get sender signals
        sender_df = adata.obsm['commot-cellchat-sum-sender']
        lr_cols = [c for c in sender_df.columns if c.startswith('s-') and '-total' not in c]
        
        for col in lr_cols:
            lr_pair = col.replace('s-', '')
            signal = sender_df[col].values
            
            # Calculate mean signal per zone
            tumor_signal = np.mean(signal[tumor_mask]) if tumor_mask.sum() > 0 else 0
            interface_signal = np.mean(signal[interface_mask]) if interface_mask.sum() > 0 else 0
            lymphoid_signal = np.mean(signal[lymphoid_mask]) if lymphoid_mask.sum() > 0 else 0
            
            lr_zone_signals[lr_pair]['Tumor'].append(tumor_signal)
            lr_zone_signals[lr_pair]['Interface'].append(interface_signal)
            lr_zone_signals[lr_pair]['Lymphoid'].append(lymphoid_signal)
    
    return lr_zone_signals


def calculate_zone_statistics(lr_zone_signals):
    """Calculate zone statistics for each L-R pair."""
    
    lr_stats = []
    
    for lr_pair, zone_signals in lr_zone_signals.items():
        # Average across samples
        tumor_mean = np.mean(zone_signals['Tumor'])
        interface_mean = np.mean(zone_signals['Interface'])
        lymphoid_mean = np.mean(zone_signals['Lymphoid'])
        total_signal = tumor_mean + interface_mean + lymphoid_mean
        
        # Standard deviation across samples
        tumor_std = np.std(zone_signals['Tumor'])
        interface_std = np.std(zone_signals['Interface'])
        lymphoid_std = np.std(zone_signals['Lymphoid'])
        
        # Count samples with signal
        n_samples = len(zone_signals['Tumor'])
        
        lr_stats.append({
            'LR_Pair': lr_pair,
            'Tumor': tumor_mean,
            'Interface': interface_mean,
            'Lymphoid': lymphoid_mean,
            'Total': total_signal,
            'Tumor_std': tumor_std,
            'Interface_std': interface_std,
            'Lymphoid_std': lymphoid_std,
            'N_Samples': n_samples,
        })
    
    return pd.DataFrame(lr_stats)


def assign_zones(df, min_interface_signal):
    """
    Assign zones using "High Interface Signal" approach.
    
    Zone assignment:
    - Interface: High Interface signal (regardless of peak)
    - Tumor: Peak in Tumor AND low Interface
    - Lymphoid: Peak in Lymphoid AND low Interface
    """
    
    def assign_zone(row):
        t, i, l = row['Tumor'], row['Interface'], row['Lymphoid']
        
        if i >= min_interface_signal:
            return 'Interface'
        elif t >= l:
            return 'Tumor'
        else:
            return 'Lymphoid'
    
    def get_pattern(row):
        """Describe the signal pattern."""
        t, i, l = row['Tumor'], row['Interface'], row['Lymphoid']
        threshold = 0.002
        
        t_active = t > threshold
        i_active = i > threshold
        l_active = l > threshold
        
        if t_active and i_active and l_active:
            return 'All zones'
        elif t_active and i_active:
            return 'Tumor→Interface'
        elif i_active and l_active:
            return 'Interface→Lymphoid'
        elif t_active and l_active:
            return 'Tumor & Lymphoid'
        elif i_active:
            return 'Interface only'
        elif t_active:
            return 'Tumor only'
        elif l_active:
            return 'Lymphoid only'
        else:
            return 'Low signal'
    
    df['Zone'] = df.apply(assign_zone, axis=1)
    df['Pattern'] = df.apply(get_pattern, axis=1)
    df['Interface_Pct'] = (df['Interface'] / df['Total'] * 100).round(1)
    
    return df


def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*70)
    print("CLIP STEP 2: ZONE MAPPING")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"  COMMOT directory: {args.commot_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Zone boundaries: Tumor <{args.tumor_end}, Interface {args.tumor_end}-{args.interface_end}, Lymphoid >{args.interface_end}")
    print(f"  Min total signal: {args.min_total_signal}")
    print(f"  Min Interface signal: {args.min_interface_signal}")
    
    # Collect L-R signals
    print("\n" + "="*70)
    print("Loading COMMOT results and calculating zone signals...")
    print("="*70)
    
    lr_zone_signals = collect_lr_signals(
        args.commot_dir, args.samples, args.tumor_end, args.interface_end
    )
    print(f"\n  Total L-R pairs across samples: {len(lr_zone_signals)}")
    
    # Calculate statistics
    print("\n" + "="*70)
    print("Calculating zone statistics...")
    print("="*70)
    
    df = calculate_zone_statistics(lr_zone_signals)
    
    print(f"\n  Total L-R pairs: {len(df)}")
    print(f"\n  Signal distribution (Total):")
    print(f"    Min: {df['Total'].min():.6f}")
    print(f"    25%: {df['Total'].quantile(0.25):.6f}")
    print(f"    50%: {df['Total'].quantile(0.50):.6f}")
    print(f"    75%: {df['Total'].quantile(0.75):.6f}")
    print(f"    Max: {df['Total'].max():.6f}")
    
    # Filter by total signal
    df_filtered = df[df['Total'] >= args.min_total_signal].copy()
    print(f"\n  L-R pairs with Total >= {args.min_total_signal}: {len(df_filtered)}")
    
    # Assign zones
    print("\n" + "="*70)
    print("Assigning zones (High Interface Signal approach)...")
    print("="*70)
    
    df_filtered = assign_zones(df_filtered, args.min_interface_signal)
    
    print("\nZone distribution:")
    print(df_filtered['Zone'].value_counts())
    
    print("\nPattern distribution:")
    print(df_filtered['Pattern'].value_counts())
    
    # Show top L-R pairs per zone
    for zone in ['Tumor', 'Interface', 'Lymphoid']:
        zone_df = df_filtered[df_filtered['Zone'] == zone].copy()
        zone_df = zone_df.sort_values(zone, ascending=False)
        
        print(f"\n{'='*70}")
        print(f"{zone.upper()} - Top 15 (n={len(zone_df)})")
        print("="*70)
        print(f"{'L-R Pair':<35} {'Tumor':>8} {'Interface':>10} {'Lymphoid':>10} {'Int%':>6}")
        print("-"*75)
        
        for _, row in zone_df.head(15).iterrows():
            print(f"{row['LR_Pair']:<35} {row['Tumor']:>8.4f} {row['Interface']:>10.4f} {row['Lymphoid']:>10.4f} {row['Interface_Pct']:>5.1f}%")
    
    # Save results
    print("\n" + "="*70)
    print("Saving results...")
    print("="*70)
    
    # Save all filtered L-R pairs
    df_all = df_filtered.sort_values(['Zone', 'Interface'], ascending=[True, False])
    df_all.to_csv(os.path.join(args.output_dir, 'LR_all_zones.csv'), index=False)
    print(f"  Saved: LR_all_zones.csv ({len(df_all)} pairs)")
    
    # Save zone-specific files
    for zone in ['Tumor', 'Interface', 'Lymphoid']:
        zone_df = df_filtered[df_filtered['Zone'] == zone].sort_values(zone, ascending=False)
        zone_df.to_csv(os.path.join(args.output_dir, f'LR_{zone}.csv'), index=False)
        print(f"  Saved: LR_{zone}.csv ({len(zone_df)} pairs)")
    
    # Save summary stats
    summary = {
        'Total_LR_pairs': len(df),
        'Filtered_LR_pairs': len(df_filtered),
        'Tumor_pairs': len(df_filtered[df_filtered['Zone'] == 'Tumor']),
        'Interface_pairs': len(df_filtered[df_filtered['Zone'] == 'Interface']),
        'Lymphoid_pairs': len(df_filtered[df_filtered['Zone'] == 'Lymphoid']),
        'Min_total_signal': args.min_total_signal,
        'Min_interface_signal': args.min_interface_signal,
        'Tumor_boundary': args.tumor_end,
        'Interface_boundary': args.interface_end,
    }
    pd.DataFrame([summary]).to_csv(os.path.join(args.output_dir, 'summary_stats.csv'), index=False)
    print(f"  Saved: summary_stats.csv")
    
    print("\n" + "="*70)
    print("STEP 2 COMPLETE!")
    print("="*70)


if __name__ == '__main__':
    main()
