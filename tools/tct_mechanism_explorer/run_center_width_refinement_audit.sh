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

# Compatibility shim for the current native_feedback_controller_audit.py branch.
# That module generates a Bash launcher inside a Python f-string and currently
# contains literal shell parameter expansions such as ${TMPDIR:-...}. Python
# interprets the braces as f-string replacement fields before Bash ever sees
# them. Supplying formatter objects through builtins makes those two fields
# round-trip back to the intended literal Bash syntax without changing any
# audit parameters, M3D-C1 inputs, or actuator physics.
python3 - "$REPO/tools/tct_mechanism_explorer/center_width_refinement_audit.py" <<'PY'
import builtins
import runpy
import sys
from pathlib import Path


class _LiteralShellParameter:
    def __init__(self, name: str) -> None:
        self.name = name

    def __format__(self, spec: str) -> str:
        return "{" + self.name + ":" + spec + "}"


builtins.TMPDIR = _LiteralShellParameter("TMPDIR")
builtins.OMPI_MCA_orte_tmpdir_base = _LiteralShellParameter(
    "OMPI_MCA_orte_tmpdir_base"
)

audit = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(audit.parent))
runpy.run_path(str(audit), run_name="__main__")
PY

echo
echo "Center-width refinement summary:"
cat validation_runs/m3dc1_tct_center_width_refinement/center_width_refinement_summary.json
