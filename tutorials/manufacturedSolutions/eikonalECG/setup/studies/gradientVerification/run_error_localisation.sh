#!/usr/bin/env bash
# Single manufactured eikonal solve that retains its time directory, so the
# cellwise activation-time error can be correlated against position.
#
# Discriminates whether the boundary is the SOURCE of the ~1.5 order (error
# injected at the domain edge and then spread inward by the elliptic solve and
# the outward characteristic sweep) or whether the shortfall is interior
# truncation error. See analyse_error_localisation.py for the test.
#
# This is a single-case (default N=40, leastSquares) spatial deep dive, kept
# for direct/manual use; it is not converted to a driverFOAM sweep since its
# whole point is the ONE retained time directory's Cx/Cy/Cz + cellwise error
# fields, not a resolution ladder. For the N-swept, per-scheme quantitative
# counterpart -- the same writeErrorField bulk/boundary split, reduced to
# L2_bulk/L2_boundary/L2_total across N=10..80 -- see the canonical
# ../errorLocalisation/run_bulk_boundary_tet.sh, which drives the same
# manufacturedEikonalECG entry through driverFOAM.
#
# Unlike run_eikonal_tet_generic.sh this does NOT call Allclean at the end.
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
N="${N:-40}"
SCHEME="${SCHEME:-leastSquares}"

cd "$CASE_DIR"

# Restore every dictionary this script rewrites, whatever happens.
BAK_SOL="$(mktemp)"; cp system/fvSolution "$BAK_SOL"
BAK_SCH="$(mktemp)"; cp system/fvSchemes  "$BAK_SCH"
BAK_EP="$(mktemp)";  cp constant/electroProperties "$BAK_EP"
trap 'cp "$BAK_SOL" system/fvSolution; cp "$BAK_SCH" system/fvSchemes; \
      cp "$BAK_EP" constant/electroProperties; rm -f "$BAK_SOL" "$BAK_SCH" "$BAK_EP"' EXIT

./Allclean >/dev/null 2>&1

LC=$($PY -c "print(1.0/$N)")
sed "s|__LC__|$LC|" setup/mesh/tet/box.geo.template > setup/mesh/tet/box.geo
gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1
gmshToFoam box.msh >/dev/null 2>&1
rm -f box.msh

# tet fvSolution overlay, matching the convergence ladder.
cp setup/mesh/tet/fvSolution system/fvSolution

# gradSchemes default only; laplacianSchemes carries its own
# "Gauss linear corrected" that must not be touched.
if [[ "$SCHEME" == "leastSquares" ]]; then GRAD="leastSquares"; else GRAD="Gauss linear"; fi
sed -E "/^gradSchemes/,/^}/ s|^([[:space:]]*)default([[:space:]]+)[^;]*;.*$|\1default\2${GRAD};|" \
    system/fvSchemes > system/fvSchemes.tmp
mv system/fvSchemes.tmp system/fvSchemes

# Ask the verifier for the cellwise signed error field.
$PY - "$CASE_DIR/constant/electroProperties" <<'PYEOF'
import re, sys
p = sys.argv[1]
s = open(p).read()
if "writeErrorField" not in s:
    s = re.sub(
        r"(verificationModel\s*\{[^}]*?enabled\s+yes;)",
        r"\1\n        writeErrorField yes;",
        s,
        count=1,
    )
    open(p, "w").write(s)
    print("enabled writeErrorField")
else:
    print("writeErrorField already present")
PYEOF
grep -A4 "verificationModel" constant/electroProperties | sed 's/^/    /'

echo "Running N=$N scheme=$SCHEME (this retains the time directory)"
decomposePar > log.localisation 2>&1
mpirun --oversubscribe -np 6 cardiacFoam -parallel >> log.localisation 2>&1
RC=$?
reconstructPar >> log.localisation 2>&1
if [[ $RC -ne 0 ]]; then
    echo "cardiacFoam exit=$RC; see $CASE_DIR/log.localisation" >&2
    exit 1
fi

# Cell centres for the positional correlation.
postProcess -func writeCellCentres -latestTime >> log.localisation 2>&1

LATEST="$(foamListTimes -latestTime 2>/dev/null | tail -1)"
if [[ -z "$LATEST" || ! -d "$LATEST" ]]; then
    echo "Could not identify the latest time directory." >&2
    exit 1
fi

echo "Latest time directory: $CASE_DIR/$LATEST"
ls "$LATEST" | sed 's/^/    /'
echo
"$SCRIPT_DIR/analyse_error_localisation.py" "$CASE_DIR/$LATEST" 10
