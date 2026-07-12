#!/bin/bash
set -euo pipefail

CASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$CASE_DIR"

set +eu
source /Volumes/OpenFOAM-v2412/etc/bashrc
set -eu

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

CONTROL_BACKUP="$(mktemp)"
cp system/controlDict "$CONTROL_BACKUP"
restore_control_dict()
{
    cp "$CONTROL_BACKUP" system/controlDict
    rm -f "$CONTROL_BACKUP"
}
trap restore_control_dict EXIT

for N in "${RESOLUTIONS[@]}"
do
    echo "=== bath-bidomain tet N=$N ==="
    OUT_DIR="setup/results/N${N}"

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

    bash setup/run_mesh_gate.sh "$N"

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

python3 setup/summarize_tet.py setup/results \
    --resolutions "${RESOLUTIONS[@]}" \
    --out setup/results/summary.csv
