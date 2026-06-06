# GLP-1R Molecular Dynamics Analysis Pipeline

## Overview

This repository contains a complete molecular dynamics (MD) analysis pipeline for GLP-1 receptor (GLP-1R) membrane protein systems simulated over 500 ns using GROMACS.

The pipeline is designed for reproducible structural, dynamical, and energetic analysis of GPCR–ligand systems and supports comparative evaluation across multiple ligands.

---

## Features

The workflow includes:

* Trajectory preprocessing (PBC correction, centering, and fitting)
* Root Mean Square Deviation (RMSD) analysis
* Root Mean Square Fluctuation (RMSF) analysis
* Radius of gyration (Rg)
* Solvent Accessible Surface Area (SASA)
* Hydrogen bond analysis
* Minimum distance and COM distance calculations
* Contact analysis and occupancy
* Contact maps (MDMAT)
* Principal Component Analysis (PCA)
* Free Energy Landscape (FEL)
* GROMOS clustering
* Dynamic Cross-Correlation Matrix (DCCM)
* Ligand RMSF and B-factor mapping
* Ramachandran analysis

---

## Requirements

* GROMACS (≥ 2022 recommended)
* Python (≥ 3.8)
* MDAnalysis
* NumPy
* Matplotlib

Install Python dependencies:

```bash
pip install numpy matplotlib MDAnalysis
```

---

## Input Files

The pipeline expects:

* `md_0_500.tpr` → GROMACS run input file
* `md_0_500.xtc` → production trajectory

---

## Usage

Run the full pipeline:

```bash
bash md_analysis_full.sh
```

All outputs are stored in:

```text
full_md_analysis/
```

---

## Key Outputs

| Analysis   | Output File                            |
| ---------- | -------------------------------------- |
| RMSD       | `rmsd_backbone.xvg`, `ligand_rmsd.xvg` |
| RMSF       | `rmsf_residue.xvg`                     |
| SASA       | `protein_sasa.xvg`, `ligand_sasa.xvg`  |
| H-bonds    | `hbonds.xvg`                           |
| PCA        | `pca_proj.xvg`, `eigenval.xvg`         |
| FEL        | `free_energy.xpm`                      |
| Clustering | `clusters.pdb`, `cluster.log`          |
| DCCM       | `dccm_matrix.npy`                      |

---

## Methodology Notes

* Trajectories are processed using PBC removal and structural alignment prior to analysis.
* PCA is performed on Cα atoms to capture dominant collective motions.
* DCCM is computed from covariance of Cα positional fluctuations and normalized to produce correlation values between -1 and +1.
* Clustering is performed using the GROMOS method with a 0.3 nm cutoff.

---

## Applications

This pipeline is suitable for:

* GPCR ligand stability analysis
* Drug discovery and lead optimization
* Comparative MD of multiple ligands
* Membrane protein conformational studies

---

## Author

Developed for GLP-1R computational studies and GPCR ligand dynamics analysis.

---

## License

MIT License
