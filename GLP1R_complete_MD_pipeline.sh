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

DUMP_TIME=500000          # ps (500 ns) - final full-trajectory snapshot
SNAPSHOT_INTERVAL=100000  # ps (100 ns) - interval for periodic snapshots
CLUSTER_CUTOFF=0.3
DCCM_STRIDE=10
TICA_STRIDE=10             # subsample every Nth frame when loading (XTC decompression
                            # of every frame in a large file is the slow part, not the math)
TICA_LAGTIME=10            # in STRIDED frames — i.e. lag = TICA_LAGTIME * TICA_STRIDE
                            # original-trajectory frames (see Step TICA for details)

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
# STEP 5 SNAPSHOT (final frame, 500 ns)
# ------------------------------------------------------
if ! skip_if_done "$OUTDIR/500ns_snapshot.pdb"; then
echo "19" | gmx trjconv \
-s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-o "$OUTDIR/500ns_snapshot.pdb" \
-dump "$DUMP_TIME" -quiet
fi

# ------------------------------------------------------
# STEP 5B SNAPSHOTS EVERY 100 ns (0, 100, 200, 300, 400, 500 ns)
# ------------------------------------------------------
SNAPDIR="$OUTDIR/snapshots"
mkdir -p "$SNAPDIR"
t=0
while [ "$t" -le "$DUMP_TIME" ]; do
    ns=$(( t / 1000 ))
    SNAP_FILE="$SNAPDIR/snapshot_${ns}ns.pdb"
    if ! skip_if_done "$SNAP_FILE"; then
        echo "19" | gmx trjconv \
        -s "$TPR" -f "$OUTDIR/fit.xtc" \
        -n "$OUTDIR/index.ndx" \
        -o "$SNAP_FILE" \
        -dump "$t" -quiet
    fi
    t=$(( t + SNAPSHOT_INTERVAL ))
done

# RMSD
if ! skip_if_done "$OUTDIR/rmsd_proteinss.xvg"; then
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
# NOTE: recent GROMACS (2023+) rewrote "gmx hbond" with a new selection-based
# interface that dropped -hbm/-hbn. The old interface (interactive groups +
# these flags) is preserved as a separate tool: "gmx hbond-legacy".
if ! skip_if_done "$OUTDIR/hbonds.xvg"; then
printf "1\n13\n" | gmx hbond-legacy \
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

# =========================================================
# DSSP (Secondary Structure)
# =========================================================
if ! skip_if_done "$OUTDIR/dssp.dat"; then
gmx dssp -s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-sel 'group "Protein"' \
-o "$OUTDIR/dssp.dat" \
-num "$OUTDIR/dssp_count.xvg" \
-tu ns
fi

# PCA (Cα atoms — group 3 in the default GROMACS index; more standard for
# essential-dynamics/PCA than full backbone (N,Ca,C) or main-chain (+O),
# and matches convention widely used in GPCR literature)
if ! skip_if_done "$OUTDIR/eigenval.xvg"; then
printf "3\n3\n" | gmx covar \
-s "$TPR" -f "$OUTDIR/fit.xtc" \
-n "$OUTDIR/index.ndx" \
-o "$OUTDIR/eigenval.xvg" \
-v "$OUTDIR/eigenvec.trr" \
-av "$OUTDIR/average.pdb"
fi

if ! skip_if_done "$OUTDIR/pca_2d.xvg"; then
printf "3\n3\n" | gmx anaeig \
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

## CLUSTER (also switched to Cα for consistency with the PCA atom selection)
if ! skip_if_done "$OUTDIR/cluster.log"; then
printf "3\n3\n" | gmx cluster \
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

export DCCM_STRIDE

if ! skip_if_done "$OUTDIR/dccm_matrix.npy"; then

python3 << PYEOF
import os
import numpy as np
import MDAnalysis as mda

stride = int(os.environ.get("DCCM_STRIDE", 10))

print("Loading trajectory...")

u = mda.Universe(
    "$TPR",
    "$OUTDIR/fit.xtc"
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
out_file = "$OUTDIR/dccm_matrix.npy"
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

plt.savefig("$OUTDIR/dccm_matrix.png", dpi=300)
plt.close()

print("DCCM heatmap saved.")

PYEOF

fi

# =========================================================
# tICA (time-lagged Independent Component Analysis)
# =========================================================
# Requires: mdtraj, deeptime  (pip install mdtraj deeptime --break-system-packages)

# MDTraj cannot read .tpr files directly (unlike MDAnalysis, used in the DCCM
# step above) — it only supports pdb/gro/psf/etc. Generate a plain PDB
# structure from the tpr once, and use that as the "topology" for MDTraj.
MDTRAJ_TOP="$OUTDIR/topol_for_mdtraj.pdb"
if ! skip_if_done "$MDTRAJ_TOP"; then
    gmx editconf -f "$TPR" -o "$MDTRAJ_TOP"
fi

if ! skip_if_done "$OUTDIR/receptor_tica_landscape.png"; then

if ! python3 -c "import mdtraj, deeptime" 2>/dev/null; then
    echo "⚠️  Skipping tICA: 'mdtraj' and/or 'deeptime' not installed."
    echo "    Install with: pip install mdtraj deeptime --break-system-packages"
else

python3 << PYEOF
import mdtraj as md
from deeptime.decomposition import TICA
import matplotlib.pyplot as plt
import numpy as np

# ==============================================================================
# 1. LOAD TOPOLOGY AND TRAJECTORY
# ==============================================================================
print("Loading trajectory...")
# MDTraj is highly optimized for fast geometric featurization
# Using C-alpha atoms ensures clean tracking of global activation states
top_file = "$MDTRAJ_TOP"
traj_file = "$OUTDIR/fit.xtc"

# fit.xtc contains the FULL system (protein + ligand + membrane + water +
# ions, ~145k atoms) across tens of thousands of frames. Loading all of it
# would require tens of GB of RAM, and XTC's per-frame compression means
# every single frame must be decompressed even when only Cα atoms are kept
# — with ~59k frames in a large file this can take a long time with no
# progress output in between. atom_indices avoids the memory blowup;
# stride cuts the number of frames actually decompressed.
topology = md.load_topology(top_file)
ca_indices = topology.select("name CA")
print(f"Cα atoms selected: {len(ca_indices)}")

tica_stride = $TICA_STRIDE
traj = md.load(traj_file, top=top_file, atom_indices=ca_indices, stride=tica_stride)
print(f"Loaded {traj.n_frames} frames with {traj.n_atoms} atoms (Cα only).")

# ==============================================================================
# 2. FEATURE EXTRACTION (C-alpha coordinates)
# ==============================================================================
print("Extracting features...")
# traj now contains ONLY Cα atoms, so no further indexing is needed here.
traj.superpose(traj, 0)

# Reshape the (n_frames, n_ca_atoms, 3) array into a flat (n_frames, n_features) matrix
features = traj.xyz.reshape(traj.n_frames, -1)
print(f"Feature matrix generated with shape: {features.shape}")

# ==============================================================================
# 3. COMPUTE TIME-LAGGED INDEPENDENT COMPONENT ANALYSIS (tICA)
# ==============================================================================
print("Fitting tICA model...")
# Lag time (tau): measured in frames (steps in this trajectory).
# Choose a lag time shorter than your target slow process but longer than fast noise.
lag_frame_time = $TICA_LAGTIME

estimator = TICA(lagtime=lag_frame_time, dim=2)  # dim=2 extracts TIC1 and TIC2
tica_model = estimator.fit_transform(features)

tic1 = tica_model[:, 0]
tic2 = tica_model[:, 1]

# Save the projected coordinates for later reuse (e.g. combining with DCCM/PCA)
np.save("$OUTDIR/tica_projection.npy", tica_model)

print("tICA execution complete.")

# ==============================================================================
# 4. PLOT RECEPTOR CONFORMATIONAL KINETIC LANDSCAPE
# ==============================================================================
print("Generating tICA projection plot...")
plt.figure(figsize=(8, 6))
sc = plt.scatter(tic1, tic2, c=np.arange(len(tic1)), cmap="viridis", s=5, alpha=0.6)
cbar = plt.colorbar(sc)
cbar.set_label("Simulation Frame Index", fontsize=12)
plt.xlabel("TIC 1 (Slowest Kinetic Mode)", fontsize=12)
plt.ylabel("TIC 2 (Second Slowest Kinetic Mode)", fontsize=12)
plt.title("Receptor Activation Trajectory in tICA Space", fontsize=14, fontweight="bold")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("$OUTDIR/receptor_tica_landscape.png", dpi=300)
plt.close()

print("tICA landscape saved.")
PYEOF

fi
fi

echo "Analysis completed."
