#!/bin/bash
# Reproduce the reported conformal tetrahedral bath-interface ladder
# (tbl-bath-bidomain-tet): the matchedSubmesh assembly with the
# distanceWeightedHarmonic interface interpolation, at N=10,20,40,80.
#
# Delegates to run_parallel_interface_sweep.sh, which runs each resolution
# under MPI (decomposePar -> mpirun -np ${NPROCS:-6} cardiacFoam -parallel ->
# reconstructPar) and writes one bathBidomainInterfaceMetrics.csv per N into
#   setup/mesh/tet/interfaceStudy/matchedSubmesh/distanceWeightedHarmonic/N*/
# which aggregate.py's `bath_tet` adapter canonicalises. N=80 (6.8M cells)
# requires the parallel path; N<=40 agree with their serial baseline to CSV
# precision (see run_parallel_equivalence.sh / parallelEquivalence/).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RESOLUTIONS="${RESOLUTIONS:-10 20 40 80}" \
ASSEMBLY="${ASSEMBLY:-matchedSubmesh}" \
METHODS="${METHODS:-distanceWeightedHarmonic}" \
    bash "$SCRIPT_DIR/run_parallel_interface_sweep.sh"
