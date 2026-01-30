#!/usr/bin/env python3
"""
================================================================================
CLIP Example Usage
================================================================================

This script demonstrates how to run the complete CLIP pipeline.

Prerequisites:
    1. HD Visium h5ad files with spatial coordinates
    2. ONTraC trajectory scores
    3. scRNA-seq reference expression data
    4. CARD cell type deconvolution results

================================================================================
"""

import os
import subprocess

# =============================================================================
# CONFIGURATION - Edit these paths for your data
# =============================================================================

# Input data paths
INPUT_DIR = '/path/to/your/h5ad_files/'
TRAJECTORY_FILE = '/path/to/your/trajectory_scores.csv'
SCRNA_FILE = '/path/to/your/scrnaseq_pct_expression.csv'
INTERACTION_FILE = '/path/to/your/cellcell_interactions.csv'

# Output directories
OUTPUT_BASE = '/path/to/output/'
COMMOT_DIR = os.path.join(OUTPUT_BASE, 'commot_results')
ZONE_DIR = os.path.join(OUTPUT_BASE, 'zone_mapping')
CLIP_DIR = os.path.join(OUTPUT_BASE, 'clip_scores')
FIG_DIR = os.path.join(OUTPUT_BASE, 'figures')

# Sample names
SAMPLES = ['MA-1', 'MA-4', 'MA-7', 'MA-9', 'MA-11', 'MA-12']

# h5ad file pattern
H5AD_PATTERN = 'Sample_11236-{sample}_8um_with_niches.h5ad'

# =============================================================================
# STEP 1: Run COMMOT
# =============================================================================
print("\n" + "="*70)
print("STEP 1: COMMOT Analysis")
print("="*70)

cmd = [
    'python', 'scripts/01_run_commot.py',
    '--input_dir', INPUT_DIR,
    '--output_dir', COMMOT_DIR,
    '--trajectory_file', TRAJECTORY_FILE,
    '--samples', *SAMPLES,
    '--dis_thr', '50',
    '--n_subsample', '30000',
    '--h5ad_pattern', H5AD_PATTERN,
]
print(f"\nRunning: {' '.join(cmd)}\n")
# subprocess.run(cmd, check=True)

# =============================================================================
# STEP 2: Zone Mapping
# =============================================================================
print("\n" + "="*70)
print("STEP 2: Zone Mapping")
print("="*70)

cmd = [
    'python', 'scripts/02_zone_mapping.py',
    '--commot_dir', COMMOT_DIR,
    '--output_dir', ZONE_DIR,
    '--samples', *SAMPLES,
    '--tumor_end', '0.36',
    '--interface_end', '0.725',
]
print(f"\nRunning: {' '.join(cmd)}\n")
# subprocess.run(cmd, check=True)

# =============================================================================
# STEP 3: CLIP Scoring
# =============================================================================
print("\n" + "="*70)
print("STEP 3: CLIP Scoring")
print("="*70)

cmd = [
    'python', 'scripts/03_clip_scoring.py',
    '--commot_dir', COMMOT_DIR,
    '--scrna_file', SCRNA_FILE,
    '--interaction_file', INTERACTION_FILE,
    '--output_dir', CLIP_DIR,
    '--samples', *SAMPLES,
    '--peak_sigma', '0.08',
    '--min_expression', '0.02',
]
print(f"\nRunning: {' '.join(cmd)}\n")
# subprocess.run(cmd, check=True)

# =============================================================================
# STEP 4: Visualization
# =============================================================================
print("\n" + "="*70)
print("STEP 4: Visualization")
print("="*70)

cmd = [
    'python', 'scripts/04_visualization.py',
    '--clip_results', os.path.join(CLIP_DIR, 'CLIP_scores.csv'),
    '--commot_dir', COMMOT_DIR,
    '--output_dir', FIG_DIR,
    '--samples', *SAMPLES,
    '--formats', 'pdf', 'png', 'svg',
]
print(f"\nRunning: {' '.join(cmd)}\n")
# subprocess.run(cmd, check=True)

# =============================================================================
# Done!
# =============================================================================
print("\n" + "="*70)
print("CLIP PIPELINE COMPLETE!")
print("="*70)
print(f"\nOutputs:")
print(f"  COMMOT results: {COMMOT_DIR}")
print(f"  Zone mapping: {ZONE_DIR}")
print(f"  CLIP scores: {CLIP_DIR}")
print(f"  Figures: {FIG_DIR}")
