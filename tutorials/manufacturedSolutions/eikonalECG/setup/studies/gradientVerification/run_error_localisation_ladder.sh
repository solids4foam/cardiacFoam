#!/usr/bin/env bash
# Eikonal error-localisation ladder: the same single-case deep dive as
# run_error_localisation.sh, swept over N and reduced to one archivable CSV.
#
# Why this exists. The Results section excludes propagation as the mechanism
# limiting the eikonal rate, on the grounds that the accumulation ramp decays
# under refinement while the near-wall elevation does not. That argument needs
# the correlations at more than one level, and it needs them on disk: the
# single-case script prints them to stdout and nothing is retained, so the
# claim previously had no archived support.
#
# Why bash and not driverFOAM. The measurement is not a resolution ladder over
# a solver-reported norm, which is what the driverFOAM sweep machinery reduces.
# It needs the retained time directory's cell centres and the cellwise signed
# error field, then a positional correlation over them by a separate utility.
# driverFOAM has no stage that keeps a time directory and hands it to an
# external analysis, so this follows run_gradient_verification.sh instead: a
# bash sweep that owns its own meshing, dictionary edits and restore.
#
# Output:
#   results/eikonal_error_localisation.csv
#     N, scheme, n_cells, interior_cut, spearman_d_seed_all,
#     spearman_d_wall_all, spearman_d_seed_interior, mean_abs_err, max_abs_err
#   results/localisation_<scheme>_N<N>.log   full stdout per level
#
# Usage:
#   ./run_error_localisation_ladder.sh
#   RESOLUTIONS="20 40" ./run_error_localisation_ladder.sh
#   SCHEME=gaussLinear ./run_error_localisation_ladder.sh
set +e
if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
        source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
    else
        echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
        exit 2
    fi
fi
# Darwin/SIP strips DYLD_LIBRARY_PATH across a fresh bash exec; RunFunctions
# restores it from FOAM_LD_LIBRARY_PATH.
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
PY="${PYTHON:-python3}"
RESOLUTIONS="${RESOLUTIONS:-20 40 80}"
SCHEME="${SCHEME:-leastSquares}"
NPROCS="${NPROCS:-6}"
OUT_DIR="$SCRIPT_DIR/results"
CSV="$OUT_DIR/eikonal_error_localisation.csv"

mkdir -p "$OUT_DIR"
# Start a fresh ladder rather than appending to a previous one.
rm -f "$CSV"

cd "$CASE_DIR"

BAK_SOL="$(mktemp)"; cp system/fvSolution "$BAK_SOL"
BAK_SCH="$(mktemp)"; cp system/fvSchemes  "$BAK_SCH"
BAK_EP="$(mktemp)";  cp constant/electroProperties "$BAK_EP"
trap 'cp "$BAK_SOL" system/fvSolution; cp "$BAK_SCH" system/fvSchemes; \
      cp "$BAK_EP" constant/electroProperties; \
      rm -f "$BAK_SOL" "$BAK_SCH" "$BAK_EP"' EXIT

if [[ "$SCHEME" == "leastSquares" ]]; then GRAD="leastSquares"; else GRAD="Gauss linear"; fi

for N in $RESOLUTIONS; do
    echo "=== N=$N scheme=$SCHEME ==="
    LOG="$OUT_DIR/localisation_${SCHEME}_N${N}.log"

    ./Allclean >/dev/null 2>&1

    LC=$($PY -c "print(1.0/$N)")
    sed "s|__LC__|$LC|" setup/mesh/tet/box.geo.template > setup/mesh/tet/box.geo
    gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1
    gmshToFoam box.msh >/dev/null 2>&1
    rm -f box.msh

    cp setup/mesh/tet/fvSolution system/fvSolution

    # gradSchemes default only. laplacianSchemes carries its own
    # "Gauss linear corrected" that must not be rewritten.
    sed -E "/^gradSchemes/,/^}/ s|^([[:space:]]*)default([[:space:]]+)[^;]*;.*$|\1default\2${GRAD};|" \
        system/fvSchemes > system/fvSchemes.tmp
    mv system/fvSchemes.tmp system/fvSchemes

    # Ask the verifier for the cellwise signed error field.
    $PY - "$CASE_DIR/constant/electroProperties" <<'PYEOF' >/dev/null
import re, sys
p = sys.argv[1]
s = open(p).read()
if "writeErrorField" not in s:
    s = re.sub(r"(verificationModel\s*\{[^}]*?enabled\s+yes;)",
               r"\1\n        writeErrorField yes;", s, count=1)
    open(p, "w").write(s)
PYEOF

    set +e
    decomposePar > "$LOG" 2>&1
    mpirun --oversubscribe -np "$NPROCS" cardiacFoam -parallel >> "$LOG" 2>&1
    RC=$?
    reconstructPar >> "$LOG" 2>&1
    set -e
    if [[ $RC -ne 0 ]]; then
        echo "FAILED N=$N: cardiacFoam exit=$RC; see $LOG" >&2
        exit 1
    fi

    postProcess -func writeCellCentres -latestTime >> "$LOG" 2>&1

    LATEST="$(foamListTimes -latestTime 2>/dev/null | tail -1)"
    if [[ -z "$LATEST" || ! -d "$LATEST" ]]; then
        echo "FAILED N=$N: no time directory retained; see $LOG" >&2
        exit 1
    fi

    LOCALISATION_CSV="$CSV" LOCALISATION_N="$N" LOCALISATION_SCHEME="$SCHEME" \
        "$SCRIPT_DIR/analyse_error_localisation.py" "$CASE_DIR/$LATEST" 10 \
        | tee -a "$LOG"

    echo "  N=$N done"
done

echo
echo "Ladder complete: $CSV"
[[ -f "$CSV" ]] && cat "$CSV"
