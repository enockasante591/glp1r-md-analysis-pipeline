# GLP-1R Molecular Dynamics Analysis Pipeline

## Overview

This repository contains a comprehensive molecular dynamics (MD) analysis pipeline for **GLP-1 receptor (GLP-1R) membrane protein systems** simulated over **500 ns** using **GROMACS**.

The pipeline provides a fully automated and reproducible workflow for preprocessing trajectories, structural and dynamical analyses, conformational sampling, free energy characterization, kinetic modeling, and protein–ligand interaction analysis. It is designed for comparative studies of GPCR–ligand systems and is suitable for large-scale molecular dynamics investigations.

---

## Features

The pipeline performs the following analyses:

### Trajectory Preprocessing
- Periodic boundary condition (PBC) correction
- System centering
- Structural fitting/alignment
- Automatic index generation
- Final structure extraction
- Periodic trajectory snapshots every 100 ns (0–500 ns)

### Structural Stability
- Root Mean Square Deviation (RMSD)
  - Protein backbone
  - Ligand
- Root Mean Square Fluctuation (RMSF)
- Radius of Gyration (Rg)
- Solvent Accessible Surface Area (SASA)

### Protein–Ligand Interactions
- Hydrogen bond analysis
- Minimum distance calculations
- Center-of-mass (COM) distance
- Contact occupancy analysis
- Contact maps
- Ligand atomic RMSF
- Ligand B-factor mapping

### Conformational Dynamics
- Principal Component Analysis (PCA)
- Free Energy Landscape (FEL)
- GROMOS clustering
- Dynamic Cross-Correlation Matrix (DCCM)
- Time-lagged Independent Component Analysis (tICA)

### Structural Characterization
- DSSP secondary structure analysis
- Ramachandran plot analysis

---

## Requirements

### Software

- GROMACS (2023 or later recommended)
- Python 3.8+

### Python Packages

- NumPy
- Matplotlib
- MDAnalysis
- MDTraj
- deeptime

Install the required Python packages:

```bash
pip install numpy matplotlib MDAnalysis mdtraj deeptime
```

> **Note:** For some Linux distributions, installation may require:

```bash
pip install mdtraj deeptime --break-system-packages
```

---

## Input Files

The pipeline expects the following input files:

```
md_0_500.tpr
md_0_500.xtc
```

where

- **md_0_500.tpr** – GROMACS portable run input
- **md_0_500.xtc** – Production MD trajectory

---

## Usage

Run the complete workflow using:

```bash
bash glp1r_md_pipeline.sh
```

All analysis results will be written automatically to:

```
full_md_analysis/
```

---

## Output Structure

```
full_md_analysis/
│
├── snapshots/
│   ├── snapshot_0ns.pdb
│   ├── snapshot_100ns.pdb
│   ├── ...
│   └── snapshot_500ns.pdb
│
├── rmsd_backbone.xvg
├── ligand_rmsd.xvg
├── rmsf_residue.xvg
├── gyration.xvg
├── protein_sasa.xvg
├── ligand_sasa.xvg
├── hbonds.xvg
├── contact_map.xpm
├── contact_occupancy.xvg
├── pca_proj.xvg
├── eigenval.xvg
├── free_energy.xpm
├── cluster.log
├── clusters.pdb
├── dccm_matrix.npy
├── dccm_matrix.png
├── receptor_tica_landscape.png
├── tica_projection.npy
├── dssp.dat
├── dssp_count.xvg
├── rama.xvg
└── pipeline.log
```

---

## Key Outputs

| Analysis | Output |
|----------|--------|
| RMSD | `rmsd_backbone.xvg`, `ligand_rmsd.xvg` |
| RMSF | `rmsf_residue.xvg` |
| Radius of Gyration | `gyration.xvg` |
| SASA | `protein_sasa.xvg`, `ligand_sasa.xvg` |
| Hydrogen Bonds | `hbonds.xvg` |
| Contact Analysis | `contact_occupancy.xvg`, `contact_map.xpm` |
| PCA | `pca_proj.xvg`, `eigenval.xvg` |
| Free Energy Landscape | `free_energy.xpm` |
| Clustering | `clusters.pdb`, `cluster.log` |
| DCCM | `dccm_matrix.npy`, `dccm_matrix.png` |
| tICA | `tica_projection.npy`, `receptor_tica_landscape.png` |
| DSSP | `dssp.dat`, `dssp_count.xvg` |
| Ramachandran Analysis | `rama.xvg` |
| Ligand Flexibility | `ligand_rmsf_atoms.xvg`, `ligand_bfac.pdb` |

---

## Methodology

- Trajectories undergo periodic boundary correction, centering, and structural alignment before analysis.
- Principal Component Analysis (PCA) is performed using **Cα atoms** to capture dominant collective motions.
- Dynamic Cross-Correlation Matrices (DCCM) are computed from Cα positional fluctuations and normalized to produce correlation coefficients ranging from **−1 to +1**.
- Time-lagged Independent Component Analysis (tICA) is employed to identify the slowest kinetic motions governing receptor conformational changes.
- Free Energy Landscapes (FEL) are generated from PCA projections using the GROMACS `sham` module.
- Conformational clustering is performed using the **GROMOS** algorithm with a **0.3 nm** RMSD cutoff.
- Secondary structure evolution is analyzed using DSSP throughout the simulation.

---

## Applications

This pipeline is suitable for:

- GPCR molecular dynamics studies
- Protein–ligand stability analysis
- Membrane protein simulations
- Drug discovery and lead optimization
- Comparative analysis of multiple ligands
- Conformational landscape characterization
- Dynamic correlation analysis
- Kinetic modeling of receptor activation

---

## Compatibility

The workflow has been developed for modern versions of **GROMACS (2023+)**. Since newer releases replaced the legacy hydrogen bond interface, the pipeline uses **`gmx hbond-legacy`** for hydrogen bond calculations.

---

## Author

Developed for **GLP-1 receptor computational studies**, molecular dynamics simulations, and GPCR ligand dynamics analysis.

---

## License

This project is licensed under the **MIT License**.
