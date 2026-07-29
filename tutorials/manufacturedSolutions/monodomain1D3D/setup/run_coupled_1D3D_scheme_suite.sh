#!/bin/bash
# Run the canonical coupled 1D-3D MMS scheme matrix.
#
# The low-level run_coupling1D3D_hex.sh runs one configuration.  This wrapper
# runs the configurations needed to compare PVJ assembly and retrograde coupling:
#
#   couplingMode:        unidirectional, bidirectional
#   pvjCouplingScheme:   explicit, implicit
#   solutionAlgorithm:   implicit
#
# Environment overrides forwarded to the low-level sweep:
#   ENDTIME, RPVJ, PVJRADIUS, COUPLED_1D3D_PAIRS
#
# Suite-specific overrides:
#   SUITE_OUTPUT_DIR           manifest directory
#   INCLUDE_TISSUE_EXPLICIT=1  add explicit tissue-solver checks with
#                              pvjCouplingScheme explicit only

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

SUITE_OUTPUT_DIR="${SUITE_OUTPUT_DIR:-$CASE_DIR/outputs/coupled1D3DSchemeSuite}"
INCLUDE_TISSUE_EXPLICIT="${INCLUDE_TISSUE_EXPLICIT:-0}"

mkdir -p "$SUITE_OUTPUT_DIR"
MANIFEST="$SUITE_OUTPUT_DIR/manifest.csv"

printf "case_id,couplingMode,pvjCouplingScheme,solutionAlgorithm,outputDir,status\n" > "$MANIFEST"

declare -a CASES=(
    "uni_pvjExplicit:unidirectional:explicit:implicit"
    "uni_pvjImplicit:unidirectional:implicit:implicit"
    "bi_pvjExplicit:bidirectional:explicit:implicit"
    "bi_pvjImplicit:bidirectional:implicit:implicit"
)

if [[ "$INCLUDE_TISSUE_EXPLICIT" == "1" ]]; then
    CASES+=(
        "uni_tissueExplicit_pvjExplicit:unidirectional:explicit:explicit"
        "bi_tissueExplicit_pvjExplicit:bidirectional:explicit:explicit"
    )
fi

for entry in "${CASES[@]}"
do
    IFS=: read -r CASE_ID COUPLING_MODE_VALUE PVJ_SCHEME_VALUE SOLUTION_ALGORITHM_VALUE <<< "$entry"

    OUTPUT_SUFFIX="_${CASE_ID}"
    OUTPUT_DIR="$CASE_DIR/outputs/coupled1D3DConvergence${OUTPUT_SUFFIX}"

    echo
    echo "=== Scheme suite case: ${CASE_ID} ==="
    echo "    couplingMode=${COUPLING_MODE_VALUE}"
    echo "    pvjCouplingScheme=${PVJ_SCHEME_VALUE}"
    echo "    solutionAlgorithm=${SOLUTION_ALGORITHM_VALUE}"
    echo "    outputDir=${OUTPUT_DIR}"

    status="ok"
    if ! COUPLING_MODE="$COUPLING_MODE_VALUE" \
         PVJ_COUPLING_SCHEME="$PVJ_SCHEME_VALUE" \
         SOLUTION_ALGORITHM="$SOLUTION_ALGORITHM_VALUE" \
         SKIP_PAPERI_AGGREGATE=1 \
         OUTPUT_SUFFIX="$OUTPUT_SUFFIX" \
         "$SCRIPT_DIR/run_coupling1D3D_hex.sh"
    then
        status="failed"
    fi

    printf "%s,%s,%s,%s,%s,%s\n" \
        "$CASE_ID" \
        "$COUPLING_MODE_VALUE" \
        "$PVJ_SCHEME_VALUE" \
        "$SOLUTION_ALGORITHM_VALUE" \
        "$OUTPUT_DIR" \
        "$status" \
        >> "$MANIFEST"

    if [[ "$status" != "ok" ]]; then
        echo "ERROR: scheme suite case ${CASE_ID} failed" >&2
        exit 1
    fi
done

echo
echo "Coupled 1D-3D scheme suite manifest written to $MANIFEST"
