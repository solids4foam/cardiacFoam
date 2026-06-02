#!/bin/bash
# runComparisons.sh
#
# Helper script to run the full set of batched comparisons on the Niederer benchmark.
# This runs both the mode comparison and the substep sweep sequentially.
#
# Usage:
#   ./runComparisons.sh
#

set -euo pipefail

CASE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "================================================================="
echo " Starting Full Comparison Run"
echo "================================================================="
echo "Step 1: Running mode comparisons (cpu vs batched integrators)..."
bash "$CASE_ROOT/run_niederer_bueno_orovio_batched_comparison.sh"

echo ""
echo "================================================================="
echo "Step 2: Running substep sweep (tradeoff analysis)..."
bash "$CASE_ROOT/run_substep_sweep.sh"

echo ""
echo "================================================================="
echo " All comparisons finished successfully!"
echo " Results and metrics: comparisonResults/"
echo " Plots              : comparisonResults/plots/"
echo "================================================================="
