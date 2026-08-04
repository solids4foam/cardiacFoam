#!/usr/bin/env bash
# Per-electrode pseudo-ECG convergence on the generic Delaunay tet family.
#
# SUPERSEDED by run_mono_per_electrode_tet.sh in this same directory, which
# drives the identical physics through driverFOAM's manufacturedFDA entry
# (sweep_tet_per_electrode.json) and reads driverFOAM's own sweepCases/
# archive instead of hand-rescuing manufacturedPseudoECGSummary.dat before
# Allclean. Kept here, still runnable directly, for manual/ad hoc use --
# see that script's header for exactly what it replaces and why.
#
# The standard sweep (setup/mesh/tet/run_monodomain_tet_generic.sh) collapses
# the verifier's per-electrode table to a single max-over-electrodes value
# before the summary file is deleted by Allclean.  The reported tetrahedral
# pseudo-ECG rate is therefore set by whichever electrode happens to be worst,
# with no way to tell which, or whether the others behave differently.
#
# The five electrodes are not equidistant from the domain:
#
#     E1 (-0.50, 0.50, 0.50)   distance to [0,1]^3 = 0.500
#     E2 ( 1.50, 0.50, 0.50)                        0.500
#     E3 ( 1.20, 0.23, 0.61)                        0.200   <-- closest
#     E4 ( 1.35, 0.74, 0.28)                        0.350
#     E5 ( 1.55, 0.41, 0.83)                        0.550
#
# The lead-field kernel is 1/|x - r_e|, so its gradient scales as 1/|x - r_e|^2
# and E3 sees a kernel roughly six times sharper than E1/E2.  If midpoint
# quadrature against a sharp kernel is what limits the observed rate, E3 should
# be the outlier and the remaining four should be better behaved.
#
# This script keeps manufacturedPseudoECGSummary.dat for every level so that
# question can be answered.  The underlying monodomain field on this family
# converges cleanly at ~2.03, so any scatter here is attributable to the
# observation functional rather than to a contaminated source field.
set +e
if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]; then
        source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null 2>&1
    else
        echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
        exit 2
    fi
fi
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
PY="${PYTHON:-python3}"
RESOLUTIONS="${RESOLUTIONS:-10 20 40 80}"
SCHEME="${SCHEME:-leastSquares}"
OUT="$SCRIPT_DIR/results/per_electrode"

cd "$CASE_DIR"
mkdir -p "$OUT"

# Restore every dictionary this script rewrites, whatever happens.
BAK_SOL="$(mktemp)"; cp system/fvSolution         "$BAK_SOL"
BAK_SCH="$(mktemp)"; cp system/fvSchemes          "$BAK_SCH"
BAK_CD="$(mktemp)";  cp system/controlDict        "$BAK_CD"
BAK_EP="$(mktemp)";  cp constant/electroProperties "$BAK_EP"
trap 'cp "$BAK_SOL" system/fvSolution; cp "$BAK_SCH" system/fvSchemes; \
      cp "$BAK_CD" system/controlDict; cp "$BAK_EP" constant/electroProperties; \
      rm -f "$BAK_SOL" "$BAK_SCH" "$BAK_CD" "$BAK_EP"' EXIT

# The pseudo-ECG analytic reference in manufacturedFDAReference.H differentiates
# the FDA manufactured shape function, so the ECG error only converges when the
# solved field is the FDA one.  Commit 48a60ae6 repurposed this case to the
# rotated-anisotropy verifier and dropped the ECG block, which is why the
# as-committed case cannot reproduce the reported pseudo-ECG ladder.  Point
# ELECTRO_PROPERTIES at a matching configuration to reproduce it -- the
# bundled configs/electroProperties.fdaDiagonalECG this used to default
# nothing at (there was no default) has been removed as redundant: it was a
# byte-for-byte copy of what the manufacturedFDA driverFOAM entry's own
# defaults (verification_model_type, conductivity, ECG electrodes/quadrature
# in core/defaults/manufactured_fda.py) already produce -- see
# run_mono_per_electrode_tet.sh, which needs no such file. Direct use of
# *this* script still requires an explicit ELECTRO_PROPERTIES override.
if [[ -n "${ELECTRO_PROPERTIES:-}" ]]; then
    if [[ ! -f "$ELECTRO_PROPERTIES" ]]; then
        echo "ELECTRO_PROPERTIES not found: $ELECTRO_PROPERTIES" >&2
        exit 2
    fi
    cp "$ELECTRO_PROPERTIES" constant/electroProperties
    echo "using electroProperties override: $ELECTRO_PROPERTIES"
fi
echo "  verifier    : $(awk '/verificationModel/{f=1} f&&/type/{print $2; exit}' constant/electroProperties)"
echo "  conductivity: $(awk '/^ *conductivity/{$1="";print;exit}' constant/electroProperties)"
echo "  ecgDomains  : $(grep -c ecgDomains constant/electroProperties)"

cp setup/mesh/tet/fvSolution system/fvSolution

dt_for_n(){ case "$1" in
    10) echo 0.00892857;; 20) echo 0.00224215;;
    40) echo 0.000560538;; 80) echo 0.000140174;;
esac; }

# gradSchemes default only. laplacianSchemes carries its own
# "Gauss linear corrected" that must not be rewritten.
if [[ "$SCHEME" == "leastSquares" ]]; then GRAD="leastSquares"; else GRAD="Gauss linear"; fi
sed -E "/^gradSchemes/,/^}/ s|^([[:space:]]*)default([[:space:]]+)[^;]*;.*$|\1default\2${GRAD};|" \
    system/fvSchemes > system/fvSchemes.tmp
mv system/fvSchemes.tmp system/fvSchemes

for N in $RESOLUTIONS; do
    ./Allclean >/dev/null 2>&1
    LC=$($PY -c "print(1.0/$N)")
    sed "s|__LC__|$LC|" setup/mesh/tet/box.geo.template > setup/mesh/tet/box.geo
    gmsh -3 setup/mesh/tet/box.geo -o box.msh -format msh2 >/dev/null 2>&1
    gmshToFoam box.msh >/dev/null 2>&1
    rm -f box.msh

    DT=$(dt_for_n "$N")
    sed -E "s/^deltaT.*/deltaT    $DT;/; s/^endTime.*/endTime    0.2;/" \
        system/controlDict > system/controlDict.tmp
    mv system/controlDict.tmp system/controlDict

    echo "=== N=$N scheme=$SCHEME ==="
    decomposePar > log.perElectrode 2>&1

    # mpirun must not be under 'set -e': its exit status is inspected below,
    # and back-to-back launches occasionally fail to bring up the PMIx server
    # ("listener thread failed to start"), which is transient and worth a retry
    # rather than aborting a multi-hour sweep.
    set +e
    RC=1
    for attempt in 1 2 3; do
        mpirun --oversubscribe -np 6 cardiacFoam -parallel >> log.perElectrode 2>&1
        RC=$?
        if [[ $RC -eq 0 ]]; then break; fi
        if grep -q "PMIx server" log.perElectrode; then
            echo "  MPI startup failed (attempt $attempt); retrying in 15 s" >&2
            sleep 15
            continue
        fi
        break
    done
    reconstructPar >> log.perElectrode 2>&1
    set -e

    SUMMARY=postProcessing/manufacturedPseudoECGSummary.dat
    if [[ $RC -ne 0 || ! -f "$SUMMARY" ]]; then
        echo "FAILED N=$N: cardiacFoam exit=$RC, summary present: $([ -f "$SUMMARY" ] && echo yes || echo no)" >&2
        echo "  see $CASE_DIR/log.perElectrode" >&2
        exit 1
    fi

    # Preserve the per-electrode table before Allclean removes it.
    cp "$SUMMARY" "$OUT/manufacturedPseudoECGSummary_N${N}.dat"
    NCELLS=$(awk '/^Grid spacing/{print $NF}' postProcessing/3D_*_cells_implicit.dat 2>/dev/null | head -1)
    echo "  kept $OUT/manufacturedPseudoECGSummary_N${N}.dat  (dx=${NCELLS:-?})"
done

echo
echo "Per-electrode summaries in $OUT"
"$SCRIPT_DIR/analyse_per_electrode.py" "$OUT" $RESOLUTIONS
