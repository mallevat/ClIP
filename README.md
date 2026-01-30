# CLIP: Cross-referenced Ligand-receptor Interaction Peak Scoring

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**CLIP** is a computational method for assigning source and target cell types to ligand-receptor (L-R) pairs by integrating single-cell expression data with spatial co-localization patterns along tissue trajectories.

## Overview

Standard approaches assign cell types to L-R pairs based solely on expression levels (e.g., "which cell type expresses the most CXCL10?"). However, the highest-expressing cell types may not be spatially co-localized at the site of L-R activity.

**CLIP solves this** by cross-referencing:
1. **Expression probability** from scRNA-seq (which cell types *can* send/receive the signal?)
2. **Spatial peak matching** from spatial transcriptomics (which cell types *actually co-localize* where the signal peaks?)

![CLIP Workflow](docs/clip_workflow.png)

## The CLIP Score

For each candidate source→target cell type pair, CLIP computes:

```
CLIP Score = √(E_L × E_R) × exp(−d²/2σ²)
```

Where:
- **E_L**: Fraction of source cells expressing the ligand gene (from scRNA-seq)
- **E_R**: Fraction of target cells expressing the receptor gene (from scRNA-seq)
- **d**: Distance between L-R signal peak and cell-cell co-localization peak along the trajectory
- **σ**: Gaussian decay bandwidth (default: 0.08)

## Why CLIP Works

**Example: CXCL10-CXCR3**

| Cell Pair | E_L | E_R | d | CLIP Score |
|-----------|-----|-----|---|------------|
| Mac→DC | 27% | 49% | 0.29 | 0.001 |
| Mal→DC | 21% | 49% | 0.00 | **0.319** |

Despite macrophages expressing more CXCL10 (27% vs 21%), **Mal→DC scores 300× higher** because malignant-DC co-localization peaks at the same trajectory position as the L-R signal (d=0.00), while macrophage-DC co-localization peaks in the Lymphoid zone (d=0.29).

**Biological interpretation**: Tumor cells at the Interface are the relevant CXCL10 source, not the higher-expressing macrophages in the Lymphoid zone.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/CLIP.git
cd CLIP

# Create conda environment
conda create -n clip_env python=3.10
conda activate clip_env

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

- scanpy >= 1.9.0
- pandas >= 1.5.0
- numpy >= 1.23.0
- scipy >= 1.9.0
- matplotlib >= 3.6.0
- seaborn >= 0.12.0
- commot >= 0.0.3

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│  Step 1: COMMOT Analysis                                        │
│  - Spatially-constrained L-R communication inference            │
│  - Input: HD Visium h5ad files + CellChatDB                     │
│  - Output: Per-spot sender/receiver scores                      │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  Step 2: Zone Mapping                                           │
│  - Assign L-R pairs to Tumor/Interface/Lymphoid zones           │
│  - Input: COMMOT results + trajectory scores                    │
│  - Output: Zone-annotated L-R pairs                             │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  Step 3: CLIP Scoring                                           │
│  - Cross-reference L-R peaks with cell-cell co-localization     │
│  - Input: L-R signals + scRNA-seq expression + CARD proportions │
│  - Output: Cell type assignments with CLIP scores               │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  Step 4: Visualization                                          │
│  - Rolling Gaussian heatmaps                                    │
│  - Chord network diagrams                                       │
│  - Cross-reference panels                                       │
└─────────────────────────────────────────────────────────────────┘
```

## Quick Start

```python
from clip import CLIPScorer

# Initialize scorer
scorer = CLIPScorer(
    commot_dir='path/to/commot_results/',
    scrna_expression='path/to/scrnaseq_expression.csv',
    trajectory_scores='path/to/trajectory_scores.csv',
    card_proportions='path/to/card_proportions.csv'
)

# Run CLIP scoring
results = scorer.score_all_lr_pairs()

# Get top assignments for a specific L-R pair
scorer.get_assignment('CXCL10-CXCR3')
# Output: {'Source': 'Malignant.cells', 'Target': 'Dendritic.cells', 'Score': 0.319}
```

## Usage

### Step 1: Run COMMOT

```bash
python scripts/01_run_commot.py \
    --input_dir /path/to/h5ad_files/ \
    --output_dir /path/to/commot_output/ \
    --samples MA-1 MA-4 MA-7 MA-9 MA-11 MA-12 \
    --dis_thr 50
```

### Step 2: Zone Mapping

```bash
python scripts/02_zone_mapping.py \
    --commot_dir /path/to/commot_output/ \
    --trajectory_file /path/to/trajectory_scores.csv \
    --output_dir /path/to/zone_output/
```

### Step 3: CLIP Scoring

```bash
python scripts/03_clip_scoring.py \
    --commot_dir /path/to/commot_output/ \
    --scrna_file /path/to/scrnaseq_pct_expression.csv \
    --interaction_file /path/to/cellcell_interactions.csv \
    --output_dir /path/to/clip_output/ \
    --peak_sigma 0.08
```

### Step 4: Visualization

```bash
python scripts/04_visualization.py \
    --clip_results /path/to/clip_output/CLIP_scores.csv \
    --output_dir /path/to/figures/
```

## Configuration

Edit `config/parameters.yaml` to customize:

```yaml
# COMMOT parameters
commot:
  dis_thr: 50              # Spatial distance threshold (µm)
  heteromeric: true        # Allow heteromeric receptors
  n_subsample: 30000       # Spots per sample

# Zone boundaries
zones:
  tumor_end: 0.36          # Tumor/Interface boundary
  interface_end: 0.725     # Interface/Lymphoid boundary

# CLIP scoring
clip:
  sigma: 0.05              # Rolling Gaussian bandwidth
  peak_sigma: 0.08         # Peak matching decay
  n_eval_points: 150       # Trajectory evaluation points
  min_expression: 0.02     # Minimum expression threshold (2%)

# Signal filtering
filtering:
  min_total_signal: 0.005
  min_interface_signal: 0.005
```

## Input Data Requirements

### 1. Spatial Transcriptomics (HD Visium)
- AnnData h5ad files with:
  - `.X`: Gene expression matrix
  - `.obsm['spatial']`: Spatial coordinates
  - `.obs['trajectory_score']`: ONTraC trajectory scores (or computed separately)

### 2. scRNA-seq Reference
- CSV file with percentage of cells expressing each gene per cell type:
```
Gene,Malignant.cells,T.cells,Macrophages,Dendritic.cells,...
CXCL10,21.3,5.2,26.6,18.4,...
CXCR3,2.1,32.4,8.7,49.4,...
```

### 3. Cell Type Deconvolution (CARD)
- Per-spot cell type proportions from CARD deconvolution

### 4. Trajectory Scores
- ONTraC-derived niche trajectory scores (0 = Tumor, 1 = Lymphoid)

## Output Files

```
output/
├── CLIP_scores.csv              # Full results with all candidate pairs
├── CLIP_top_assignments.csv     # Top assignment per L-R pair
├── LR_by_zone/
│   ├── Tumor_LR_pairs.csv
│   ├── Interface_LR_pairs.csv
│   └── Lymphoid_LR_pairs.csv
└── figures/
    ├── Heatmap_Tumor.pdf
    ├── Heatmap_Interface.pdf
    ├── Heatmap_Lymphoid.pdf
    ├── ChordNetwork.pdf
    └── CrossReference_Panel.pdf
```

## Citation

If you use CLIP in your research, please cite:

```bibtex
@article{allevato2026clip,
  title={CLIP: Cross-referenced Ligand-receptor Interaction Peak scoring for 
         spatially-resolved cell-cell communication analysis},
  author={Allevato, Michael and others},
  journal={Nature Cancer},
  year={2026}
}
```

## Methods References

1. **COMMOT**: Cang Z, et al. Screening cell–cell communication in spatial transcriptomics via collective optimal transport. *Nat Methods* 20, 218-228 (2023).

2. **CellChatDB**: Jin S, et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun* 12, 1088 (2021).

3. **ONTraC**: Wang W, et al. ONTraC characterizes spatially continuous variations of tissue microenvironment. *Genome Biol* 26, 117 (2025).

4. **CARD**: Ma Y & Zhou X. Spatially informed cell-type deconvolution for spatial transcriptomics. *Nat Biotechnol* 40, 1349-1359 (2022).

5. **HNSCC scRNA-seq**: Choi JH, et al. Single-cell transcriptome profiling of the stepwise progression of head and neck cancer. *Nat Commun* 14, 1055 (2023).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- Michael Allevato -mallevat@umich.edu
- GitHub Issues: [https://github.com/yourusername/CLIP/issues](https://github.com/yourusername/CLIP/issues)
