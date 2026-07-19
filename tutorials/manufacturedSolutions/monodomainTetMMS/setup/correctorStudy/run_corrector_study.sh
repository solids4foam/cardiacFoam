#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
PYTHON="${PYTHON:-python3}"
RESOLUTIONS_STR="${RESOLUTIONS:-10}"
OUTER_COUNTS_STR="${OUTER_COUNTS:-1 2 3 4 8}"
NONORTH_COUNTS_STR="${NONORTH_COUNTS:-0 1 2}"
NONORTH_FIXED_OUTER="${NONORTH_FIXED_OUTER:-1}"
REPEATS="${REPEATS:-1}"
STEPS="${STEPS:-4}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
WORK_ROOT="${WORK_ROOT:-$(mktemp -d /tmp/cardiacfoam-corrector-study-XXXXXX)}"
KEEP_WORK="${KEEP_WORK:-0}"

read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"
read -r -a OUTER_COUNTS <<< "$OUTER_COUNTS_STR"
read -r -a NONORTH_COUNTS <<< "$NONORTH_COUNTS_STR"

if [[ ! -f "$OPENFOAM_BASHRC" ]]; then
    echo "OpenFOAM bashrc not found: $OPENFOAM_BASHRC" >&2
    exit 2
fi

set +eu
source "$OPENFOAM_BASHRC" >/dev/null
set -eu

for executable in cardiacFoam gmsh gmshToFoam; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "Required executable not found after sourcing OpenFOAM: $executable" >&2
        exit 2
    fi
done

cleanup()
{
    if [[ "$KEEP_WORK" != "1" ]]; then
        rm -rf "$WORK_ROOT"
    else
        echo "Keeping work directory: $WORK_ROOT"
    fi
}
trap cleanup EXIT

rm -rf "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

dt_for_n()
{
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        *) echo "No timestep configured for N=$1" >&2; exit 2 ;;
    esac
}

set_dictionary_controls()
{
    local fv_solution="$1"
    local control_dict="$2"
    local outer="$3"
    local nonorth="$4"
    local dt="$5"
    local end_time="$6"

    "$PYTHON" - "$fv_solution" "$control_dict" "$outer" "$nonorth" "$dt" "$end_time" <<'PY'
from pathlib import Path
import re
import sys

fv_path, control_path, outer, nonorth, dt, end_time = sys.argv[1:]
fv = Path(fv_path).read_text()

def replace_required(text, key, value):
    pattern = rf"(?m)^(\s*{re.escape(key)}\s+)\d+(\s*;)"
    updated, count = re.subn(pattern, rf"\g<1>{value}\g<2>", text)
    if count != 1:
        raise SystemExit(f"Expected exactly one {key} entry, found {count}")
    return updated

fv = replace_required(fv, "nOuterCorrectors", outer)
nonorth_pattern = r"(?m)^(\s*nNonOrthogonalCorrectors\s+)\d+(\s*;)"
if re.search(nonorth_pattern, fv):
    fv = re.sub(nonorth_pattern, rf"\g<1>{nonorth}\g<2>", fv)
else:
    outer_line = r"(?m)^(\s*nOuterCorrectors\s+\d+\s*;\s*)$"
    fv, count = re.subn(
        outer_line,
        rf"\1\n    nNonOrthogonalCorrectors {nonorth};",
        fv,
    )
    if count != 1:
        raise SystemExit("Could not insert nNonOrthogonalCorrectors")
Path(fv_path).write_text(fv)

control = Path(control_path).read_text()

def replace_scalar(text, key, value):
    pattern = rf"(?m)^(\s*{re.escape(key)}\s+)[^;]+(;)"
    updated, count = re.subn(pattern, rf"\g<1>{value}\g<2>", text)
    if count != 1:
        raise SystemExit(f"Expected exactly one {key} entry, found {count}")
    return updated

for key, value in (
    ("deltaT", dt),
    ("endTime", end_time),
    ("writeControl", "timeStep"),
    ("writeInterval", "1"),
    ("writePrecision", "16"),
):
    control = replace_scalar(control, key, value)
Path(control_path).write_text(control)
PY
}

clean_run_state()
{
    local run_case="$1"
    local time_dir
    for time_dir in "$run_case"/[0-9]*; do
        if [[ -d "$time_dir" ]]; then
            rm -rf "$time_dir"
        fi
    done
    rm -rf "$run_case/postProcessing" "$run_case/processor"*
    rm -f "$run_case"/log.*
}

run_one()
{
    local run_case="$1"
    local n="$2"
    local mode="$3"
    local outer="$4"
    local nonorth="$5"
    local repeat="$6"
    local dt="$7"
    local end_time="$8"
    local run_id="N${n}_${mode}_o${outer}_no${nonorth}_r${repeat}"
    local out_dir="$RESULTS_DIR/$run_id"

    clean_run_state "$run_case"
    set_dictionary_controls \
        "$run_case/system/fvSolution" \
        "$run_case/system/controlDict" \
        "$outer" "$nonorth" "$dt" "$end_time"

    mkdir -p "$out_dir"
    printf 'run_id=%s\nmode=%s\nresolution=%s\nn_outer=%s\nn_nonorth=%s\nrepeat=%s\ndelta_t=%s\nend_time=%s\n' \
        "$run_id" "$mode" "$n" "$outer" "$nonorth" "$repeat" "$dt" "$end_time" \
        > "$out_dir/metadata.env"

    echo "=== $run_id ==="
    local harness_start_ns
    local harness_end_ns
    local harness_wall_time
    harness_start_ns="$($PYTHON -c 'import time; print(time.perf_counter_ns())')"
    (
        cd "$run_case"
        cardiacFoam > log.cardiacFoam 2>&1
    )
    harness_end_ns="$($PYTHON -c 'import time; print(time.perf_counter_ns())')"
    harness_wall_time="$($PYTHON -c "print((int('$harness_end_ns') - int('$harness_start_ns')) / 1e9)")"
    printf 'harness_wall_time_s=%s\n' "$harness_wall_time" >> "$out_dir/metadata.env"

    cp "$run_case/log.cardiacFoam" "$out_dir/"
    cp "$run_case/system/fvSolution" "$out_dir/"
    cp "$run_case/system/controlDict" "$out_dir/"

    local dat
    dat="$(find "$run_case/postProcessing" -maxdepth 1 -name '3D_*_cells_implicit.dat' -print | head -1)"
    if [[ -z "$dat" ]]; then
        echo "No MMS summary produced for $run_id" >&2
        exit 3
    fi
    cp "$dat" "$out_dir/summary.dat"

    local latest_time
    latest_time="$(
        for time_dir in "$run_case"/[0-9]*; do
            if [[ -d "$time_dir" ]]; then
                basename "$time_dir"
            fi
        done | sort -g | tail -1
    )"
    if [[ -z "$latest_time" || ! -f "$run_case/$latest_time/Vm" ]]; then
        echo "No final Vm field found for $run_id" >&2
        exit 3
    fi
    cp "$run_case/$latest_time/Vm" "$out_dir/Vm"
    printf 'latest_time=%s\n' "$latest_time" >> "$out_dir/metadata.env"
}

for n in "${RESOLUTIONS[@]}"; do
    dt="$(dt_for_n "$n")"
    if [[ -n "${ENDTIME:-}" ]]; then
        end_time="$ENDTIME"
    else
        end_time="$($PYTHON -c "print(float('$dt') * int('$STEPS'))")"
    fi

    run_case="$WORK_ROOT/N$n"
    mkdir -p "$run_case"
    rsync -a \
        --exclude 'setup/results/' \
        --exclude 'setup/correctorStudy/' \
        --exclude 'constant/polyMesh/' \
        --exclude '[0-9]*/' \
        "$CASE_DIR/" "$run_case/"

    lc="$($PYTHON -c "print(1.0/float('$n'))")"
    sed "s/__LC__/$lc/" "$CASE_DIR/setup/box.geo.template" > "$run_case/setup/box.geo"
    (
        cd "$run_case"
        gmsh -3 setup/box.geo -o box.msh -format msh2 > log.gmsh 2>&1
        gmshToFoam box.msh > log.gmshToFoam 2>&1
        rm -f box.msh
        checkMesh > log.checkMesh 2>&1 || true
    )
    cp "$run_case/log.checkMesh" "$RESULTS_DIR/N${n}_mesh.log"

    for repeat in $(seq 1 "$REPEATS"); do
        for outer in "${OUTER_COUNTS[@]}"; do
            run_one "$run_case" "$n" outer "$outer" 0 "$repeat" "$dt" "$end_time"
        done
        for nonorth in "${NONORTH_COUNTS[@]}"; do
            run_one "$run_case" "$n" nonorth "$NONORTH_FIXED_OUTER" "$nonorth" "$repeat" "$dt" "$end_time"
        done
    done
done

"$PYTHON" "$SCRIPT_DIR/summarize_corrector_study.py" "$RESULTS_DIR"
echo "Corrector study complete: $RESULTS_DIR"
