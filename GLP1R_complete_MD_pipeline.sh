#!/bin/bash
# ============================================================
# GLP-1R FULL MD ANALYSIS PIPELINE
# Membrane GPCR - 500 ns trajectory
# ============================================================

set -euo pipefail
IFS=$'\n\t'

TPR="md_0_500.tpr"
XTC="md_0_500.xtc"
OUTDIR="full_md_analysis"

DUMP_TIME=500000
CLUSTER_CUTOFF=0.3
DCCM_STRIDE=10

mkdir -p "$OUTDIR"
LOG="$OUTDIR/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

skip_if_done() {
    [[ -f "$1" ]]
}

echo "======================================================"
echo "GLP-1R FULL MD ANALYSIS PIPELINE"
date
echo "======================================================"

# ------------------------------------------------------
# STEP 1 PBC
# ------------------------------------------------------
if ! skip_if_done "$OUTDIR/nojump.xtc"; then
echo "0" | gmx trjconv -s "$TPR" -f "$XTC" \
-o "$OUTDIR/nojump.xtc" -pbc nojump -tu ns -quiet
fi

# ------------------------------------------------------
# STEP 2 INDEX
# ------------------------------------------------------
if ! skip_if_done "$OUTDIR/index.ndx"; then
gmx make_ndx -f "$TPR" -o "$OUTDIR/index.ndx" << EOF
1 | 13
name 19 Protein_Ligand
14 | 15
name 20 Membrane
1 | 13 | 14 | 15
name 21 Prot_Lig_Lipids
q
EOF
fi

# ------------------------------------------------------
# STEP 3 CENTER
# ------------------------------------------------------
if ! skip_if_done "$OUTDIR/center.xtc"; then
printf "1\n0\n" | gmx trjconv \
-s "$TPR" \
-f "$OUTDIR/nojump.xtc" \
-n "$OUTDIR/index.ndx" \
-o "$OUTDIR/center.xtc" \
-center -pbc mol -ur compact -tu ns -quiet
fi

# ------------------------------------------------------
# STEP 4 FIT
# ------------------------------------------------------
if ! skip_if_done "$OUTDIR/fit.xtc"; then
printf "19\n0\n" | gmx trjconv \
-s "$TPR" \
-f "$OUTDIR/center.xtc" \
-n "$OUTDIR/index.ndx" \
-o "$OUTDIR/fit.xtc" \
-fit rot+trans -quiet
fi

# ------------------------------------------------------
# STEP 5 SNAPSHOT
# ------------------------------------------------------
if ! skip_if_done "$OUTDIR/500ns_snapshot.pdb"; then
echo "19" | gmx trjconv \
-s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-o "$OUTDIR/500ns_snapshot.pdb" \
-dump "$DUMP_TIME" -quiet
fi

# RMSD
if ! skip_if_done "$OUTDIR/rmsd_backbone.xvg"; then
printf "1\n1\n" | gmx rms -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" -o "$OUTDIR/rmsd_backbone.xvg" -tu ns
fi

if ! skip_if_done "$OUTDIR/ligand_rmsd.xvg"; then
printf "1\n13\n" | gmx rms -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" -o "$OUTDIR/ligand_rmsd.xvg" -tu ns
fi

# RMSF
if ! skip_if_done "$OUTDIR/rmsf_residue.xvg"; then
echo "1" | gmx rmsf -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" -o "$OUTDIR/rmsf_residue.xvg" -res
fi

# Rg
if ! skip_if_done "$OUTDIR/gyration.xvg"; then
echo "1" | gmx gyrate -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" -o "$OUTDIR/gyration.xvg"
fi

# SASA
if ! skip_if_done "$OUTDIR/protein_sasa.xvg"; then
echo "1" | gmx sasa -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" -o "$OUTDIR/protein_sasa.xvg"
fi

if ! skip_if_done "$OUTDIR/ligand_sasa.xvg"; then
echo "13" | gmx sasa -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" -o "$OUTDIR/ligand_sasa.xvg"
fi

# HBONDS
if ! skip_if_done "$OUTDIR/hbonds.xvg"; then
printf "1\n13\n" | gmx hbond \
-s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-num "$OUTDIR/hbonds.xvg" \
-hbm "$OUTDIR/hbond_matrix.xpm" \
-hbn "$OUTDIR/hbond_index.ndx"
fi

# DISTANCES
if ! skip_if_done "$OUTDIR/com_distance.xvg"; then
gmx distance -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-select 'com of group "Protein" plus com of group "UNK"' \
-oall "$OUTDIR/com_distance.xvg"
fi

if ! skip_if_done "$OUTDIR/mindist.xvg"; then
printf "1\n13\n" | gmx mindist \
-s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-od "$OUTDIR/mindist.xvg"
fi

# CONTACT OCCUPANCY
if ! skip_if_done "$OUTDIR/contact_occupancy.xvg"; then
gmx select -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-select 'group "Protein" and within 0.4 of group "UNK"' \
-os "$OUTDIR/contact_occupancy.xvg"
fi

# CONTACT MAP
if ! skip_if_done "$OUTDIR/contact_map.xpm"; then
echo "1" | gmx mdmat -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" -mean "$OUTDIR/contact_map.xpm"
fi

# =========================================================
# LIGAND RMSF + B-FACTOR
# =========================================================
if ! skip_if_done "$OUTDIR/ligand_bfac.pdb"; then

echo "13" | gmx rmsf \
    -s "$TPR" \
    -f "$OUTDIR/fit.xtc" \
    -n "$OUTDIR/index.ndx" \
    -o "$OUTDIR/ligand_rmsf_atoms.xvg" \
    -oq "$OUTDIR/ligand_bfac.pdb"
fi

# RAMA
if ! skip_if_done "$OUTDIR/rama.xvg"; then
gmx rama -s "$TPR" -f "$OUTDIR/fit.xtc" -o "$OUTDIR/rama.xvg"
fi

# PCA
if ! skip_if_done "$OUTDIR/eigenval.xvg"; then
printf "4\n4\n" | gmx covar \
-s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-o "$OUTDIR/eigenval.xvg" \
-v "$OUTDIR/eigenvec.trr" \
-av "$OUTDIR/average.pdb"
fi

if ! skip_if_done "$OUTDIR/pca_2d.xvg"; then
printf "4\n4\n" | gmx anaeig \
-v "$OUTDIR/eigenvec.trr" \
-s "$TPR" \
-f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-first 1 -last 2 \
-2d "$OUTDIR/pca_2d.xvg" \
-proj "$OUTDIR/pca_proj.xvg"
fi

# FEL
if ! skip_if_done "$OUTDIR/free_energy.xpm"; then
gmx sham -f "$OUTDIR/pca_2d.xvg" \
-ls "$OUTDIR/free_energy.xpm" \
-g "$OUTDIR/sham.log" -notime
fi

## CLUSTER
if ! skip_if_done "$OUTDIR/cluster.log"; then
printf "4\n4\n" | gmx cluster \
    -s "$TPR" -f "$OUTDIR/fit.xtc" \
    -n "$OUTDIR/index.ndx" \
    -method gromos \
    -cutoff "$CLUSTER_CUTOFF" \
    -cl "$OUTDIR/clusters.pdb" \
    -g "$OUTDIR/cluster.log" \
    -skip 10
fi

# =========================================================
# DCCM (Dynamic Cross-Correlation Matrix)
# =========================================================

export DCCM_STRIDE=10

if ! skip_if_done "$OUTDIR/dccm_matrix.npy"; then

python3 << PYEOF
import os
import numpy as np
import MDAnalysis as mda

stride = int(os.environ.get("DCCM_STRIDE", 10))

print("Loading trajectory...")

u = mda.Universe(
    "md_0_500.tpr",
    "full_md_analysis/fit.xtc"
)

ca = u.select_atoms("protein and name CA")
n_res = ca.n_atoms

print(f"CA atoms: {n_res}")

# ---------------------------------------------------------
# STEP 1: Collect CA coordinates (stride-controlled)
# ---------------------------------------------------------
frames = []

for ts in u.trajectory[::stride]:
    frames.append(ca.positions.copy())

frames = np.array(frames)   # shape: (frames, residues, 3)

n_frames = frames.shape[0]

print(f"Frames used: {n_frames}")

# ---------------------------------------------------------
# STEP 2: Reshape to (time, 3N)
# ---------------------------------------------------------
X = frames.reshape(n_frames, n_res * 3)

# remove mean (fluctuations only)
X -= X.mean(axis=0)

# ---------------------------------------------------------
# STEP 3: covariance matrix
# ---------------------------------------------------------
cov = np.cov(X.T)

# reshape into residue blocks
cov = cov.reshape(n_res, 3, n_res, 3)

# trace over xyz components → correlation between residues
dccm = np.trace(cov, axis1=1, axis2=3)

# ---------------------------------------------------------
# STEP 4: normalization (true DCCM form)
# ---------------------------------------------------------
diag = np.sqrt(np.diag(dccm))
outer = np.outer(diag, diag)

with np.errstate(divide='ignore', invalid='ignore'):
    dccm = np.where(outer > 0, dccm / outer, 0.0)

# enforce physical constraints
np.fill_diagonal(dccm, 1.0)

# ---------------------------------------------------------
# STEP 5: save matrix
# ---------------------------------------------------------
out_file = "full_md_analysis/dccm_matrix.npy"
np.save(out_file, dccm)

print("DCCM saved:", out_file)

# ---------------------------------------------------------
# STEP 6 (optional but useful): quick heatmap
# ---------------------------------------------------------
import matplotlib.pyplot as plt

plt.figure(figsize=(6,5))
plt.imshow(dccm, cmap="RdBu_r", vmin=-1, vmax=1, origin="lower")
plt.colorbar(label="Correlation")
plt.title("DCCM (Cα)")
plt.tight_layout()

plt.savefig("full_md_analysis/dccm_matrix.png", dpi=300)
plt.close()

print("DCCM heatmap saved.")

PYEOF

fi

echo "Analysis completed."
