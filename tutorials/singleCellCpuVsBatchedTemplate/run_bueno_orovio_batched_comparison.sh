#!/bin/bash

set -euo pipefail

CASE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ROOT="${RUN_ROOT:-$CASE_ROOT/comparisonRuns}"
RESULT_ROOT="${RESULT_ROOT:-$CASE_ROOT/comparisonResults}"
N_RUNS="${N_RUNS:-3}"
MODES="${MODES:-cpu batched_euler batched_rl batched_soa}"
CONTROL_END_TIME="${CONTROL_END_TIME:-0.03}"
OPENFOAM_BASHRC="${CF_OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"

if [ -f "$OPENFOAM_BASHRC" ]; then
    # shellcheck disable=SC1090
    set +eu
    source "$OPENFOAM_BASHRC"
    set -eu
fi

config_for_mode()
{
    case "$1" in
        cpu)
            echo "electroProperties.cpu"
            ;;
        batched_euler)
            echo "electroProperties.batched_euler"
            ;;
        batched_rl)
            echo "electroProperties.batched_rl"
            ;;
        batched_soa)
            echo "electroProperties.batched_soa"
            ;;
        *)
            echo "Unknown mode '$1'" >&2
            return 1
            ;;
    esac
}

set_control_end_time()
{
    local control_dict="$1"
    local end_time="$2"

    python3 - "$control_dict" "$end_time" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
end_time = sys.argv[2]
text = path.read_text()
text, count = re.subn(r"(endTime\s+)[^;]+;", rf"\g<1>{end_time};", text)
if count != 1:
    raise SystemExit(f"Could not replace endTime in {path}")
path.write_text(text)
PY
}

prepare_case()
{
    local mode="$1"
    local run_id="$2"
    local config_name
    local case_dir

    config_name="$(config_for_mode "$mode")"
    case_dir="$RUN_ROOT/$mode/run$(printf '%02d' "$run_id")"

    rm -rf "$case_dir"
    mkdir -p "$case_dir"

    cp -R "$CASE_ROOT/constant" "$case_dir/"
    cp -R "$CASE_ROOT/system" "$case_dir/"
    cp "$CASE_ROOT/constant/$config_name" "$case_dir/constant/electroProperties"
    set_control_end_time "$case_dir/system/controlDict" "$CONTROL_END_TIME"
}

run_case()
{
    local mode="$1"
    local run_id="$2"
    local case_dir="$RUN_ROOT/$mode/run$(printf '%02d' "$run_id")"
    local start
    local end
    local elapsed

    start="$(python3 -c 'import time; print(time.perf_counter())')"

    (
        cd "$case_dir"
        cardiacFoam > log.cardiacFoam 2>&1
    )

    end="$(python3 -c 'import time; print(time.perf_counter())')"
    elapsed="$(python3 -c "print(${end} - ${start})")"
    printf 'real %.9f\n' "$elapsed" > "$case_dir/time.txt"
}

mkdir -p "$RUN_ROOT" "$RESULT_ROOT"

for mode in $MODES; do
    rm -rf "$RUN_ROOT/$mode"
    mkdir -p "$RUN_ROOT/$mode"

    for run_id in $(seq 1 "$N_RUNS"); do
        echo "Running $mode run $run_id/$N_RUNS"
        prepare_case "$mode" "$run_id"
        run_case "$mode" "$run_id"
    done
done

python3 "$CASE_ROOT/compare_bueno_orovio_batched.py" \
    "$RUN_ROOT" \
    "$RESULT_ROOT/comparison_metrics.json" \
    $MODES

python3 "$CASE_ROOT/plot_bueno_orovio_batched.py" \
    "$RUN_ROOT" \
    "$RESULT_ROOT/plots" \
    $MODES

echo "Comparison metrics written to $RESULT_ROOT/comparison_metrics.json"
echo "Comparison plots written to $RESULT_ROOT/plots"
