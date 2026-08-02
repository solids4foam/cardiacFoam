#!/usr/bin/env bash
# Isolated least-squares gradient reconstruction on the two tet mesh families.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
        set +eu
        source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
        set -eu
    else
        echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
        exit 2
    fi
fi
set +eu
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1
set -eu

PYTHON_BIN="${PYTHON:-python3}"
RESOLUTIONS="${RESOLUTIONS:-10 20 40 80}"
RESULT="$SCRIPT_DIR/results/eikonal_gradient_tet.csv"
GENERIC_TEMPLATE="$CASE_DIR/setup/mesh/tet/box.geo.template"
FRONTAL_TEMPLATE="$CASE_DIR/setup/mesh/tet/box.geo.template.optimised"

fv_schemes_backup="$(mktemp)"
cp "$CASE_DIR/system/fvSchemes" "$fv_schemes_backup"
trap 'cp "$fv_schemes_backup" "$CASE_DIR/system/fvSchemes"; rm -f "$fv_schemes_backup"' EXIT

mkdir -p "$(dirname "$RESULT")"
echo "mesh_family,N,h,n_cells,Linf_max,Linf_mean,n_cells_Linf_gt_0_05,L2_bulk,L2_boundary,L2_total" > "$RESULT"

run_case() {
    local mesh_family="$1"
    local n="$2"
    local template="$3"
    local h log_file n_cells linf_max linf_mean n_above l2_bulk l2_boundary l2_total

    cd "$CASE_DIR"
    ./Allclean >/dev/null 2>&1
    h="$($PYTHON_BIN -c "print(1.0/$n)")"
    sed "s|__LC__|$h|" "$template" > setup/mesh/tet/box.geo
    gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1
    gmshToFoam box.msh >/dev/null 2>&1
    sed -E \
        's/^([[:space:]]*)default([[:space:]]+)Gauss linear;.*$/\1default\2leastSquares;/; s/^([[:space:]]*)default([[:space:]]+)leastSquares.*$/\1default\2leastSquares;/' \
        system/fvSchemes > system/fvSchemes.tmp
    mv system/fvSchemes.tmp system/fvSchemes

    log_file="$SCRIPT_DIR/results/${mesh_family}_N${n}.log"
    gradientReconstructionOrder > "$log_file" 2>&1

    n_cells="$(awk '/cells[[:space:]]*:/ {print $NF; exit}' "$log_file")"
    linf_max="$(awk -F= '/E_inf =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    linf_mean="$(awk -F= '/Mean E_inf =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    n_above="$(awk -F= '/Cells with error > 0.05 =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    l2_bulk="$(awk -F= '/L2 Bulk =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    l2_boundary="$(awk -F= '/L2 Bound =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    l2_total="$(awk -F= '/L2 Total =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"

    if [[ -z "$n_cells" || -z "$l2_bulk" || -z "$l2_boundary" || -z "$l2_total" ]]; then
        echo "Incomplete gradient metrics for $mesh_family N=$n; see $log_file" >&2
        exit 1
    fi
    echo "$mesh_family,$n,$h,$n_cells,$linf_max,$linf_mean,$n_above,$l2_bulk,$l2_boundary,$l2_total" >> "$RESULT"
}

for n in $RESOLUTIONS; do
    run_case generic "$n" "$GENERIC_TEMPLATE"
    run_case frontal "$n" "$FRONTAL_TEMPLATE"
done

echo "Wrote $RESULT"
