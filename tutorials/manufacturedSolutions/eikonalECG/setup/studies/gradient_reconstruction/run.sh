#!/usr/bin/env bash
# Canonical entry for the registered "eikonal_gradient_tet" verification
# experiment (applications/scripts/driverFoam/verification_experiments.json):
#   matrix: mesh_family {tet_generic, tet_frontal} x gradient_scheme
#           {leastSquares} x N {10, 20, 40, 80}   (8 cases)
#   result: setup/results/eikonal_gradient_tet.csv
#
# gradientReconstructionOrder (applications/test/gradientReconstructionOrder/)
# exercises the gradient operator alone against an exact analytic field --
# no cardiacFoam solve, no case materialization -- so it is not something an
# existing driverFOAM tutorial entry can drive (see SPEC_FACTORIES in
# applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/tutorials/registry.py,
# none of which parametrize their workflow_dag's solve command). A genuine
# driverFOAM entry for it would need a new tutorial spec module (mirroring
# manufactured_eikonal_ecg.py) whose _workflow_dag_for swaps the "cardiacFoam"
# solve step for "gradientReconstructionOrder" -- allowed in principle, since
# AGENT_GUIDE.md documents workflow_dag steps as "any core OpenFOAM app or
# your own compiled utility" -- plus a bespoke output reader, since this
# utility's metrics (n_cells, Linf_max/mean, L2_bulk/boundary/total) are
# printed to stdout/log, not written to a postProcessing/*.dat file the way
# every cardiacFoam verifier's output is. That is new-entry-authoring work,
# not a sweep JSON against an existing entry, so per this refactor's brief
# ("do not invent new entries unless genuinely required") it has not been
# done; this thin wrapper instead reuses the working bash implementation,
# restricted to this experiment's registered (leastSquares-only) matrix.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
GRADIENT_STUDY="$CASE_DIR/setup/studies/gradientVerification"
RESULT_LOCAL="$SCRIPT_DIR/results/eikonal_gradient_tet.csv"

mkdir -p "$SCRIPT_DIR/results"
SCHEME_LABELS="leastSquares" RESULT="$RESULT_LOCAL" \
    "$GRADIENT_STUDY/run_gradient_verification.sh"

mkdir -p "$CASE_DIR/setup/results"
cp "$RESULT_LOCAL" "$CASE_DIR/setup/results/eikonal_gradient_tet.csv"

echo "Wrote $CASE_DIR/setup/results/eikonal_gradient_tet.csv"
