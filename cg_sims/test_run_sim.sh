#!/bin/bash
# Version 2.0 — Minimal GROMACS MD workflow (min → nvt → npt)
# Requires: ./mdps/{min.mdp,nvt.mdp,npt.mdp} in the current run folder
# Expects:  output_initial_state.gro, topol.top, index.ndx

set -e  # stop on first error

# Load GROMACS if needed (adjust module name if on a cluster)
module load gromacs/2021.4-AVX2-GPU 2>/dev/null || echo "GROMACS already available."

# Inputs (all relative to the current run directory)
MDP_DIR="../../mdps"
STRUCTURE="output_initial_state.gro"
TOPOLOGY="topol.top"
INDEX="index.ndx"

# Quick checks (kept lightweight)
[[ -f "$STRUCTURE" ]] || { echo "❌ Missing $STRUCTURE — run COBY first."; exit 1; }
[[ -f "$TOPOLOGY"  ]] || { echo "❌ Missing $TOPOLOGY — run COBY first.";  exit 1; }
[[ -f "$INDEX"     ]] || { echo "❌ Missing $INDEX — build index.ndx first."; exit 1; }
for f in min.mdp nvt.mdp npt.mdp; do
  [[ -f "$MDP_DIR/$f" ]] || { echo "❌ Missing $MDP_DIR/$f"; exit 1; }
done

# Run sequence
declare -a STEPS=("min" "nvt" "npt" "prod")

for STEP in "${STEPS[@]}"; do
  echo
  echo "============================================"
  echo "Running $STEP ..."
  echo "============================================"

  gmx grompp \
    -f "$MDP_DIR/${STEP}.mdp" \
    -c "$STRUCTURE" \
    -p "$TOPOLOGY" \
    -n "$INDEX" \
    -o "${STEP}.tpr" \
    -r "$STRUCTURE" \
    -maxwarn 1

  gmx mdrun -deffnm "$STEP" -v

  # Set structure for next stage
  STRUCTURE="${STEP}.gro"
done

echo
echo "✅ All stages completed successfully!"
