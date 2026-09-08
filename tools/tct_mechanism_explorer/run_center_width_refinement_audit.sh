#!/usr/bin/env bash
set -euo pipefail

REPO="/home/ubuntu/work/openmc/sweep"
cd "$REPO"

SWEEP_TMP="${TCT_MPI_TMPDIR:-/tmp/tct-$USER}"
mkdir -p "$SWEEP_TMP"
export TMPDIR="$SWEEP_TMP"
export OMPI_MCA_orte_tmpdir_base="$SWEEP_TMP"
mkdir -p "$TMPDIR" "$OMPI_MCA_orte_tmpdir_base"

source "$HOME/spack/share/spack/setup-env.sh"
spack env activate m3dc1-deps

python3 tools/tct_mechanism_explorer/center_width_refinement_audit.py

echo
echo "Center-width refinement summary:"
cat validation_runs/m3dc1_tct_center_width_refinement/center_width_refinement_summary.json

