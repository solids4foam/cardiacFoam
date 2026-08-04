#!/usr/bin/env bash
# Cartesian bath-bidomain ladder WITH the interface-current diagnostics.
#
# The existing Cartesian study (setup/studies/cartesianConvergence) reports the
# regional potentials only; it never runs bathBidomainInterfaceMetrics, so the
# assembled-current-density diagnostic exists for the tetrahedral ladder alone.
# That is the gap this script fills, and it is the controlling comparison for
# the tetrahedral result.
#
# Why it discriminates. On this geometry the myocardium occupies x in [0,1] and
# the bath the two slabs either side, so both interfaces are planar with
# x-aligned normals. The manufactured conductivities are diagonal
# (constant/electroProperties: conductivityExtracellular xy=xz=yz=0), so the
# face normal is an eigenvector of the conductivity tensor and every interface
# face is exactly K-orthogonal. The two-point distanceWeightedHarmonic
# transmissibility is exact in precisely that case.
#
# So:
#   assembled current density converges at ~2 here  -> the first-order rate
#       measured on tetrahedra is a consequence of interface non-orthogonality,
#       and the DDFV forward pointer becomes a quantified claim.
#   it converges at ~1 here too                     -> the two-point flux is
#       first-order regardless of orthogonality, and non-orthogonality is NOT
#       the mechanism.
#
# Both passes are run on the same converged case: the solved fields, and then
# -exactFields, which substitutes the manufactured phiE so the flux construction
# is evaluated without solve or coupling error in its extracellular input.
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
RESOLUTIONS="${RESOLUTIONS:-10 20 40}"
NPROCS="${NPROCS:-6}"
OUT_ROOT="$SCRIPT_DIR/results/interfaceCartesian"

cd "$CASE_DIR"
mkdir -p "$OUT_ROOT"

# Restore every dictionary this script rewrites, whatever happens.
BAK_BM="$(mktemp)";  cp system/blockMeshDict.3D "$BAK_BM"
BAK_CD="$(mktemp)";  cp system/controlDict      "$BAK_CD"
BAK_BMD="$(mktemp)"; cp system/blockMeshDict    "$BAK_BMD" 2>/dev/null || true
trap 'cp "$BAK_BM" system/blockMeshDict.3D; cp "$BAK_CD" system/controlDict; \
      cp "$BAK_BMD" system/blockMeshDict 2>/dev/null || true; \
      rm -f "$BAK_BM" "$BAK_CD" "$BAK_BMD"' EXIT

# deltaT matched to the tetrahedral ladder's dt ~ h^2 convention, and the
# matching step count so writeInterval lands exactly on the final step. The
# case default is writeControl runTime with writeInterval 0.1, which at
# endTime 0.02 never fires -- the run then completes with no time directory
# written and the metrics utility silently falls back to time 0.
dt_for_n(){ case "$1" in
    10) echo 0.00892857;; 20) echo 0.00224215;;
    40) echo 0.000560538;; 80) echo 0.000140174;;
    *) echo "Unsupported N=$1" >&2; exit 2;;
esac; }
steps_for_n(){ case "$1" in
    10) echo 2;; 20) echo 9;; 40) echo 36;; 80) echo 143;;
    *) echo "Unsupported N=$1" >&2; exit 2;;
esac; }

for N in $RESOLUTIONS; do
    OUT_DIR="$OUT_ROOT/N$N"; mkdir -p "$OUT_DIR"
    echo "=== Cartesian N=$N ==="

    ./Allclean >/dev/null 2>&1

    # All three blocks carry the same (n n n) count, so one substitution does it.
    sed -E "s/\(80 80 80\)/($N $N $N)/g" system/blockMeshDict.3D > system/blockMeshDict
    blockMesh > "$OUT_DIR/log.blockMesh" 2>&1
    # myocardium / bath cellZones -- required by bathBidomainInterfaceMetrics.
    topoSet > "$OUT_DIR/log.topoSet" 2>&1

    if ! grep -q "myocardium" "$OUT_DIR/log.topoSet"; then
        echo "FAILED N=$N: topoSet did not create the myocardium zone" >&2
        exit 1
    fi

    # The bath domain reads a bodyAndOrgansConductivity field that does not
    # exist until this utility writes it; without it the solver aborts on a
    # missing 0/bodyAndOrgansConductivity as soon as the extracellular domain
    # is constructed.
    setTorsoOrganConductivityField > "$OUT_DIR/log.setConductivity" 2>&1
    if [[ ! -f 0/bodyAndOrgansConductivity ]]; then
        echo "FAILED N=$N: setTorsoOrganConductivityField wrote no field;" \
             "see $OUT_DIR/log.setConductivity" >&2
        exit 1
    fi

    foamDictionary system/controlDict -entry deltaT        -set "$(dt_for_n "$N")"    > /dev/null
    foamDictionary system/controlDict -entry endTime       -set 0.02                  > /dev/null
    foamDictionary system/controlDict -entry writeControl  -set timeStep              > /dev/null
    foamDictionary system/controlDict -entry writeInterval -set "$(steps_for_n "$N")" > /dev/null

    set +e
    if [[ "$NPROCS" -gt 1 ]]; then
        # decomposeParDict is committed with 6 subdomains; without this the
        # solver aborts with "specifies N processors but job was started with
        # M ranks" whenever NPROCS differs from it.
        foamDictionary system/decomposeParDict \
            -entry numberOfSubdomains -set "$NPROCS" > /dev/null 2>&1 || true
        decomposePar -force > "$OUT_DIR/log.decomposePar" 2>&1
        mpirun --oversubscribe -np "$NPROCS" cardiacFoam -parallel \
            > "$OUT_DIR/log.cardiacFoam" 2>&1
        RC=$?
        reconstructPar -latestTime > "$OUT_DIR/log.reconstructPar" 2>&1
    else
        cardiacFoam > "$OUT_DIR/log.cardiacFoam" 2>&1
        RC=$?
    fi
    set -e
    if [[ $RC -ne 0 ]]; then
        echo "FAILED N=$N: cardiacFoam exit=$RC; see $OUT_DIR/log.cardiacFoam" >&2
        exit 1
    fi

    # Guard: if no time directory past 0 was written, -latestTime silently
    # selects 0 and the utility reports metrics for the initial condition (or
    # dies on a missing 0/phiE). Fail here instead.
    if [[ -z "$(find . -maxdepth 1 -type d -regex '\./0\.[0-9].*' -o -maxdepth 1 -type d -regex '\./[1-9][0-9]*' | head -1)" ]]; then
        echo "FAILED N=$N: solver wrote no time directory past 0 --" \
             "check writeControl/writeInterval against endTime" >&2
        exit 1
    fi

    # FACE_ERRORS=1 adds the signed per-face dump. It is what distinguishes a
    # spread error, which a manufactured-solution or parameter inconsistency
    # would give, from an error concentrated in a handful of faces.
    FACE_FLAG=""
    if [[ "${FACE_ERRORS:-0}" == "1" ]]; then FACE_FLAG="-writeFaceErrors"; fi

    bathBidomainInterfaceMetrics -latestTime $FACE_FLAG \
        > "$OUT_DIR/log.interfaceMetrics" 2>&1
    cp postProcessing/bathBidomainInterfaceMetrics.csv "$OUT_DIR/" 2>/dev/null || {
        echo "FAILED N=$N: no interface metrics written" >&2; exit 1; }

    bathBidomainInterfaceMetrics -latestTime -exactFields $FACE_FLAG \
        > "$OUT_DIR/log.interfaceMetricsExactField" 2>&1
    cp postProcessing/bathBidomainInterfaceMetricsExactField.csv "$OUT_DIR/" 2>/dev/null || {
        echo "FAILED N=$N: no exact-field metrics written" >&2; exit 1; }
    if [[ -n "$FACE_FLAG" ]]; then
        cp postProcessing/bathBidomainFaceErrors*.csv "$OUT_DIR/" 2>/dev/null || true
    fi

    echo "  N=$N done -> $OUT_DIR"
done

echo
echo "Cartesian interface metrics in $OUT_ROOT"
