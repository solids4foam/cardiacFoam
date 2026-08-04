#!/usr/bin/env bash
# Isolated gradient reconstruction on the two tet mesh families, both
# gradSchemes (leastSquares and Gauss linear) by default.
#
# This is the underlying implementation for the REGISTERED "eikonal_gradient_tet"
# verification experiment (applications/scripts/driverFoam/verification_experiments.json),
# but its own registered matrix only covers leastSquares (the scheme that
# converges) on both mesh families -- the canonical entry point for that
# experiment is therefore ../gradient_reconstruction/run.sh, a thin wrapper
# that calls this script with SCHEME_LABELS=leastSquares and writes to this
# tutorial's normalized setup/results/eikonal_gradient_tet.csv. Run this
# script directly (as it always has) for the full leastSquares-vs-gaussLinear
# comparison used in the paper's qualitative discussion; both schemes' output
# still lands here at studies/gradientVerification/results/eikonal_gradient_tet.csv.
#
# gradientReconstructionOrder (applications/test/gradientReconstructionOrder/)
# is a standalone utility, not a cardiacFoam solve -- no driverFOAM tutorial
# entry drives it (see SPEC_FACTORIES in
# applications/scripts/driverFoam/openfoam_driver/plugins/cardiacfoam/tutorials/registry.py),
# so there is no sweep_*.json for this study; this bash loop over mesh
# family x scheme x N remains the whole implementation.
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
# Gradient schemes to sweep, given as labels. Both are run by default so that
# the bulk/boundary decomposition exists for the scheme that converges and for
# the one that does not, on the same mesh ladder. Labels are mapped to their
# fvSchemes entry by scheme_entry() below; keep them space-free.
SCHEME_LABELS="${SCHEME_LABELS:-leastSquares gaussLinear}"
# Overridable so a restricted-scope caller (../gradient_reconstruction/run.sh,
# which sets SCHEME_LABELS=leastSquares for the registered eikonal_gradient_tet
# experiment) can write to its own path instead of overwriting this directory's
# full leastSquares-vs-gaussLinear comparison file.
RESULT="${RESULT:-$SCRIPT_DIR/results/eikonal_gradient_tet.csv}"

scheme_entry() {
    case "$1" in
        leastSquares) printf 'leastSquares' ;;
        gaussLinear)  printf 'Gauss linear' ;;
        *) echo "Unknown gradient scheme label '$1'" >&2; exit 2 ;;
    esac
}
GENERIC_TEMPLATE="$CASE_DIR/setup/mesh/tet/box.geo.template"
FRONTAL_TEMPLATE="$CASE_DIR/setup/mesh/tet/box.geo.template.optimised"

fv_schemes_backup="$(mktemp)"
cp "$CASE_DIR/system/fvSchemes" "$fv_schemes_backup"
trap 'cp "$fv_schemes_backup" "$CASE_DIR/system/fvSchemes"; rm -f "$fv_schemes_backup"' EXIT

mkdir -p "$(dirname "$RESULT")"
echo "mesh_family,scheme,N,h,n_cells,Linf_max,Linf_mean,n_cells_Linf_gt_0_05,L2_bulk,L2_boundary,L2_total" > "$RESULT"

run_case() {
    local mesh_family="$1"
    local n="$2"
    local template="$3"
    local scheme_label="$4"
    local h log_file n_cells linf_max linf_mean n_above l2_bulk l2_boundary l2_total
    local scheme

    scheme="$(scheme_entry "$scheme_label")"

    cd "$CASE_DIR"
    ./Allclean >/dev/null 2>&1
    h="$($PYTHON_BIN -c "print(1.0/$n)")"
    sed "s|__LC__|$h|" "$template" > setup/mesh/tet/box.geo
    gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1
    gmshToFoam box.msh >/dev/null 2>&1
    # Rewrite the gradSchemes default to the requested scheme. The range
    # restriction matters: laplacianSchemes also carries a
    # "default Gauss linear corrected;" entry, and rewriting that would
    # silently change the Laplacian discretisation as well as the gradient.
    # Restored from the backup by the EXIT trap.
    sed -E \
        "/^gradSchemes/,/^}/ s|^([[:space:]]*)default([[:space:]]+)[^;]*;.*$|\1default\2${scheme};|" \
        system/fvSchemes > system/fvSchemes.tmp
    mv system/fvSchemes.tmp system/fvSchemes

    log_file="$SCRIPT_DIR/results/${mesh_family}_${scheme_label}_N${n}.log"
    gradientReconstructionOrder > "$log_file" 2>&1

    n_cells="$(awk '/cells[[:space:]]*:/ {print $NF; exit}' "$log_file")"
    linf_max="$(awk -F= '/E_inf =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    linf_mean="$(awk -F= '/Mean E_inf =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    n_above="$(awk -F= '/Cells with error > 0.05 =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    l2_bulk="$(awk -F= '/L2 Bulk =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    l2_boundary="$(awk -F= '/L2 Bound =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"
    l2_total="$(awk -F= '/L2 Total =/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' "$log_file")"

    if [[ -z "$n_cells" || -z "$l2_bulk" || -z "$l2_boundary" || -z "$l2_total" ]]; then
        echo "Incomplete gradient metrics for $mesh_family $scheme_label N=$n;" \
             "see $log_file" >&2
        exit 1
    fi
    echo "$mesh_family,$scheme_label,$n,$h,$n_cells,$linf_max,$linf_mean,$n_above,$l2_bulk,$l2_boundary,$l2_total" >> "$RESULT"
}

for scheme_label in $SCHEME_LABELS; do
    for n in $RESOLUTIONS; do
        run_case generic "$n" "$GENERIC_TEMPLATE" "$scheme_label"
        run_case frontal "$n" "$FRONTAL_TEMPLATE" "$scheme_label"
    done
done

echo "Wrote $RESULT"
