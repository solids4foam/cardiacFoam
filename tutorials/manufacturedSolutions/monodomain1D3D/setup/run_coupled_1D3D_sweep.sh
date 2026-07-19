#!/bin/bash
# Joint 1D-3D convergence sweep.
# Pairs 3D mesh resolution with 1D graph so that h_3D ≈ h_1D:
#   10^3 / nodes011, 20^3 / nodes021, 40^3 / nodes041, 80^3 / nodes081
#
# Environment overrides:
#   ENDTIME   final simulation time (default 0.1)
#   RPVJ      PVJ resistance, overrides rPvj in electroProperties (default: no override)
#   COUPLING_MODE  overrides couplingMode (unidirectional or bidirectional)
#   OUTPUT_SUFFIX  appended to the output dir name (default: empty)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

ENDTIME="${ENDTIME:-0.1}"
RPVJ="${RPVJ:-}"
PVJRADIUS="${PVJRADIUS:-}"
COUPLING_MODE="${COUPLING_MODE:-}"
OUTPUT_SUFFIX="${OUTPUT_SUFFIX:-}"
OUTPUT_DIR="$CASE_DIR/outputs/coupled1D3DConvergence${OUTPUT_SUFFIX}"

_OF_BASHRC=""
if [[ -f "/Volumes/OpenFOAM-v2412/etc/bashrc" ]]; then
    _OF_BASHRC="/Volumes/OpenFOAM-v2412/etc/bashrc"
elif [[ -n "${WM_PROJECT_DIR:-}" && -f "$WM_PROJECT_DIR/etc/bashrc" ]]; then
    _OF_BASHRC="$WM_PROJECT_DIR/etc/bashrc"
fi
if [[ -n "$_OF_BASHRC" ]]; then
    set +eu
    . "$_OF_BASHRC"
    set -eu
fi
unset _OF_BASHRC

# Pairs: (N_3D graph_id)
declare -a PAIRS=(
    "10 nodes011"
    "20 nodes021"
    "40 nodes041"
    "80 nodes081"
)

mkdir -p "$OUTPUT_DIR"

# Restore backups on any exit
_restore_backups() {
    local bak="$CASE_DIR/system/controlDict.bak"
    [[ -f "$bak" ]] && mv "$bak" "$CASE_DIR/system/controlDict"
    local epbak="$CASE_DIR/constant/electroProperties.bak"
    [[ -f "$epbak" ]] && mv "$epbak" "$CASE_DIR/constant/electroProperties"
}
trap _restore_backups EXIT

for pair in "${PAIRS[@]}"
do
    read -r N_CELLS GRAPH_ID <<< "$pair"

    # Skip if pvjRadius < cell size (pvjMapper would abort)
    if [[ -n "$PVJRADIUS" ]]; then
        H=$(python3 -c "print(1.0/${N_CELLS})")
        SKIP=$(python3 -c "print('yes' if ${PVJRADIUS} < ${H} else 'no')")
        if [[ "$SKIP" == "yes" ]]; then
            echo "=== Skipping N=${N_CELLS}: pvjRadius=${PVJRADIUS} < h=${H} ==="
            continue
        fi
    fi

    echo "=== Coupled sweep: ${N_CELLS}^3 mesh / ${GRAPH_ID}  endTime=${ENDTIME}  rPvj=${RPVJ:-default}  pvjRadius=${PVJRADIUS:-default}  couplingMode=${COUPLING_MODE:-default} ==="

    # --- 3D mesh ---
    sed "s/(10 10 10)/(${N_CELLS} ${N_CELLS} ${N_CELLS})/g" \
        "$CASE_DIR/system/blockMeshDict.3D" \
        > "$CASE_DIR/system/blockMeshDict.3D.active"

    (cd "$CASE_DIR" && blockMesh -dict system/blockMeshDict.3D.active \
        > "$CASE_DIR/log.blockMesh.${N_CELLS}" 2>&1)

    # --- 1D graph ---
    "$SCRIPT_DIR/select_purkinje_graph.sh" "$GRAPH_ID" >/dev/null

    # --- clean previous fields ---
    rm -rf "$CASE_DIR"/[1-9]* "$CASE_DIR"/[0-9]*.[0-9]* 2>/dev/null || true
    rm -f  "$CASE_DIR"/postProcessing/graph_*_nodes.dat
    rm -f  "$CASE_DIR"/postProcessing/3D_*_cells_*.dat
    rm -f  "$CASE_DIR"/verification/coupled1D3DMonodomain_diagnostics.csv

    # --- dt ~ h²: anchor N=80 at DT_BASE; coarser meshes scale as (80/N)² ---
    # This matches the standalone MMS driver and ensures temporal error << h²
    DT=$(python3 -c "print(f'{1.40174e-04 * (80/${N_CELLS})**2:.8e}')")

    cp "$CASE_DIR/system/controlDict" "$CASE_DIR/system/controlDict.bak"
    sed -e "s/^deltaT.*/deltaT          ${DT};/" \
        -e "s/^endTime.*/endTime         ${ENDTIME};/" \
        "$CASE_DIR/system/controlDict.bak" \
        > "$CASE_DIR/system/controlDict"

    # --- optional electroProperties overrides ---
    if [[ -n "$RPVJ" || -n "$PVJRADIUS" || -n "$COUPLING_MODE" ]]; then
        cp "$CASE_DIR/constant/electroProperties" \
           "$CASE_DIR/constant/electroProperties.bak"
        _ep_sed_args=()
        [[ -n "$RPVJ"      ]] && _ep_sed_args+=(-e "s/rPvj[[:space:]].*[0-9];/rPvj            ${RPVJ};/")
        [[ -n "$PVJRADIUS" ]] && _ep_sed_args+=(-e "s/pvjRadius[[:space:]].*[0-9];/pvjRadius       ${PVJRADIUS};/")
        [[ -n "$COUPLING_MODE" ]] && _ep_sed_args+=(-e "s/couplingMode[[:space:]].*;/couplingMode    ${COUPLING_MODE};/")
        sed "${_ep_sed_args[@]}" \
            "$CASE_DIR/constant/electroProperties.bak" \
            > "$CASE_DIR/constant/electroProperties"
    fi

    # --- run ---
    (cd "$CASE_DIR" && cardiacFoam \
        > "$CASE_DIR/log.cardiacFoam.${N_CELLS}" 2>&1)

    mv "$CASE_DIR/system/controlDict.bak" "$CASE_DIR/system/controlDict"
    [[ -f "$CASE_DIR/constant/electroProperties.bak" ]] && \
        mv "$CASE_DIR/constant/electroProperties.bak" \
           "$CASE_DIR/constant/electroProperties"

    # --- collect outputs ---
    case_output="$OUTPUT_DIR/${N_CELLS}"
    mkdir -p "$case_output"

    cp "$CASE_DIR"/postProcessing/graph_*_nodes.dat     "$case_output/" 2>/dev/null || true
    cp "$CASE_DIR"/postProcessing/3D_*_cells_*.dat      "$case_output/" 2>/dev/null || true
    cp "$CASE_DIR"/verification/coupled1D3DMonodomain_diagnostics.csv \
        "$case_output/coupling_diagnostics.csv"          2>/dev/null || true
done

# clean temp dict
rm -f "$CASE_DIR/system/blockMeshDict.3D.active"

"$SCRIPT_DIR/post_processing_coupled_1D3D.py" --output-dir "$OUTPUT_DIR"

echo "Coupled 1D-3D sweep outputs written to $OUTPUT_DIR"

# --- Paper I: persist canonical convergence CSV (additive; does not alter the sweep above) ---
_PAPERI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
python3 "$_PAPERI_ROOT/applications/scripts/paperI_results/aggregate.py" coupling \
    --repo-root "$_PAPERI_ROOT" \
    || echo "WARN: paperI aggregate (coupling) failed; native output untouched" >&2
