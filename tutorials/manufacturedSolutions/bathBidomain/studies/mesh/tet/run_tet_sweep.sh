#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$CASE_DIR"

if [[ -z "${WM_PROJECT_DIR:-}" ]]
then
    if [[ -f /Volumes/OpenFOAM-v2412/etc/bashrc ]]
    then
        set +eu
        source /Volumes/OpenFOAM-v2412/etc/bashrc
        # Darwin/SIP strips DYLD_LIBRARY_PATH across a fresh bash exec; RunFunctions
        # restores it from FOAM_LD_LIBRARY_PATH (see its own "Darwin workaround" block).
        # Without this, cardiacFoam aborts with "Library not loaded: @rpath/libOpenFOAM.dylib".
        source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1
        set -eu
    else
        echo "OpenFOAM is not sourced. Source the v2412 etc/bashrc first." >&2
        exit 2
    fi
fi

dt_for_n()
{
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        80) echo 0.000140174 ;;
        *) echo "No deltaT configured for N=$1" >&2; exit 2 ;;
    esac
}

RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40}"
ENDTIME="${ENDTIME:-0.02}"
read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"

# Merged case (bathBidomain): system/controlDict, constant/electroProperties,
# and system/fvSchemes are the hex defaults (dimension "1D", no
# bathPotentialDomain interface-assembly entries, Gauss linear gradient) --
# activate this case's own tet overlay for all three and restore byte-for-byte
# on exit.
CONTROL_BACKUP="$(mktemp)"
cp system/controlDict "$CONTROL_BACKUP"
ELECTRO_BACKUP="$(mktemp)"
cp constant/electroProperties "$ELECTRO_BACKUP"
cp studies/mesh/tet/electroProperties constant/electroProperties
SCHEMES_BACKUP="$(mktemp)"
cp system/fvSchemes "$SCHEMES_BACKUP"
cp studies/mesh/tet/fvSchemes system/fvSchemes
restore_control_dict()
{
    cp "$CONTROL_BACKUP" system/controlDict
    rm -f "$CONTROL_BACKUP"
    cp "$ELECTRO_BACKUP" constant/electroProperties
    rm -f "$ELECTRO_BACKUP"
    cp "$SCHEMES_BACKUP" system/fvSchemes
    rm -f "$SCHEMES_BACKUP"
}
trap restore_control_dict EXIT

for N in "${RESOLUTIONS[@]}"
do
    echo "=== bath-bidomain tet N=$N ==="
    OUT_DIR="studies/mesh/tet/results/N${N}"

    if [[ "${RESUME:-1}" == "1" \
       && -f "$OUT_DIR/summary.dat" \
       && -f "$OUT_DIR/log.cardiacFoam" \
       && -f "$OUT_DIR/log.checkMesh" ]] \
       && grep -q '^End$' "$OUT_DIR/log.cardiacFoam" \
       && grep -q 'Mesh OK' "$OUT_DIR/log.checkMesh"
    then
        echo "reusing completed N=$N result"
        continue
    fi

    bash studies/mesh/tet/run_mesh_gate.sh "$N"

    rm -rf 0 postProcessing [0-9]*

    DT="$(dt_for_n "$N")"
    foamDictionary system/controlDict -entry deltaT -set "$DT"
    foamDictionary system/controlDict -entry endTime -set "$ENDTIME"
    foamDictionary system/controlDict -entry writeControl -set runTime
    foamDictionary system/controlDict -entry writeInterval -set 1e9

    setTorsoOrganConductivityField > log.setTorsoOrganConductivityField 2>&1
    cardiacFoam > log.cardiacFoam 2>&1

    DAT="$(find postProcessing -maxdepth 1 -name 'bathBidomain_3D_*_cells_implicit.dat' -print -quit)"
    if [[ -z "$DAT" ]]
    then
        echo "Missing manufactured error summary for N=$N" >&2
        exit 1
    fi

    cp "$DAT" "$OUT_DIR/summary.dat"
    cp log.setTorsoOrganConductivityField log.cardiacFoam "$OUT_DIR/"

    if grep -Eiq 'FOAM FATAL|Caught signal|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' log.cardiacFoam
    then
        echo "Solver failure signature found for N=$N" >&2
        exit 1
    fi

    echo "completed N=$N"
done

restore_control_dict
trap - EXIT

python3 studies/mesh/tet/summarize_tet.py studies/mesh/tet/results \
    --resolutions "${RESOLUTIONS[@]}" \
    --out studies/mesh/tet/results/summary.csv
