#!/bin/bash
# run_substep_sweep.sh
#
# Runs cpu (reference) + euler and soa at substep counts 50, 20, 10, 5.
# The sub100 data points are pulled from the first comparison run
# (comparisonResults/comparison_metrics.json) by plot_substep_sweep.py,
# so we skip re-running them here.
#
# Layout:
#   comparisonRuns/<mode>/run01/     (cpu, euler_sub050, euler_sub020, ...)
#   comparisonResults/substep_sweep_metrics.json
#   comparisonResults/plots/substep_*.png
#
# Env-var overrides:
#   SUBSTEPS              substep counts to sweep   (default: "50 20 10 5")
#   COMPARISON_METRICS    first comparison JSON to merge sub100 data from
#                         (default: comparisonResults/comparison_metrics.json)
#   N_RUNS                repeat runs per mode       (default: 1)
#   CF_OPENFOAM_BASHRC
#   RUN_ROOT / RESULT_ROOT

set -euo pipefail

CASE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ROOT="${RUN_ROOT:-$CASE_ROOT/comparisonRuns}"
RESULT_ROOT="${RESULT_ROOT:-$CASE_ROOT/comparisonResults}"
N_RUNS="${N_RUNS:-1}"
SUBSTEPS="${SUBSTEPS:-50 20 10 5}"
COMPARISON_METRICS="${COMPARISON_METRICS:-$RESULT_ROOT/comparison_metrics.json}"
OPENFOAM_BASHRC="${CF_OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
N_SUBDOMAINS="${N_SUBDOMAINS:-6}"

if [ -f "$OPENFOAM_BASHRC" ]; then
    set +eu; source "$OPENFOAM_BASHRC"; set -eu
fi

# ---------------------------------------------------------------------------
# Build the mode list: cpu + euler_subXXX + soa_subXXX (no sub100 here)
# ---------------------------------------------------------------------------
build_mode_list()
{
    local modes="cpu"
    for s in $SUBSTEPS; do
        modes="$modes euler_sub$(printf '%03d' "$s")"
    done
    for s in $SUBSTEPS; do
        modes="$modes soa_sub$(printf '%03d' "$s")"
    done
    echo "$modes"
}

# ---------------------------------------------------------------------------
# Map mode → electroProperties file
# ---------------------------------------------------------------------------
config_for_mode()
{
    case "$1" in
        cpu)           echo "electroProperties.cpu" ;;
        euler_sub050)  echo "electroProperties.batched_euler_sub050" ;;
        euler_sub020)  echo "electroProperties.batched_euler_sub020" ;;
        euler_sub010)  echo "electroProperties.batched_euler_sub010" ;;
        euler_sub005)  echo "electroProperties.batched_euler_sub005" ;;
        soa_sub050)    echo "electroProperties.batched_soa_sub050" ;;
        soa_sub020)    echo "electroProperties.batched_soa_sub020" ;;
        soa_sub010)    echo "electroProperties.batched_soa_sub010" ;;
        soa_sub005)    echo "electroProperties.batched_soa_sub005" ;;
        *)
            echo "Unknown mode '$1'" >&2; return 1 ;;
    esac
}

# ---------------------------------------------------------------------------
# Stage a case directory
# ---------------------------------------------------------------------------
prepare_case()
{
    local mode="$1" run_id="$2"
    local config_name case_dir

    config_name="$(config_for_mode "$mode")"
    case_dir="$RUN_ROOT/$mode/run$(printf '%02d' "$run_id")"

    rm -rf "$case_dir"
    mkdir -p "$case_dir"
    cp -R "$CASE_ROOT/constant" "$case_dir/"
    cp -R "$CASE_ROOT/system"   "$case_dir/"
    cp "$CASE_ROOT/constant/$config_name" "$case_dir/constant/electroProperties"

    echo "  Prepared $case_dir"
}

# ---------------------------------------------------------------------------
# Run one case — times cardiacFoam only
# ---------------------------------------------------------------------------
run_case()
{
    local mode="$1" run_id="$2"
    local case_dir="$RUN_ROOT/$mode/run$(printf '%02d' "$run_id")"
    local start end elapsed

    (
        cd "$case_dir"
        blockMesh      > log.blockMesh      2>&1
        decomposePar   > log.decomposePar   2>&1

        start="$(python3 -c 'import time; print(time.perf_counter())')"
        mpirun -np "$N_SUBDOMAINS" cardiacFoam -parallel > log.cardiacFoam 2>&1
        end="$(python3 -c 'import time; print(time.perf_counter())')"
        elapsed="$(python3 -c "print(${end} - ${start})")"
        printf 'real %.9f\n' "$elapsed" > time.txt

        reconstructPar > log.reconstructPar 2>&1
        postProcess -func Niedererpoints -latestTime > log.postProcess_points 2>&1
        postProcess -func Niedererlines  -latestTime > log.postProcess_lines  2>&1

        echo "  cardiacFoam wall time: ${elapsed} s"
    )
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
MODES="$(build_mode_list)"
echo "Substep sweep modes: $MODES"
echo "sub100 data will be pulled from: $COMPARISON_METRICS"

mkdir -p "$RUN_ROOT" "$RESULT_ROOT"

for mode in $MODES; do
    rm -rf "$RUN_ROOT/$mode"
    mkdir -p "$RUN_ROOT/$mode"

    for run_id in $(seq 1 "$N_RUNS"); do
        echo ""
        echo "===================================================="
        echo "  Mode: $mode  |  Run: $run_id / $N_RUNS"
        echo "===================================================="
        prepare_case "$mode" "$run_id"
        run_case     "$mode" "$run_id"
    done
done

# ---------------------------------------------------------------------------
# Compare (sweep modes only — cpu is reference, sub100 not needed here)
# ---------------------------------------------------------------------------
echo ""
echo "Running comparison ..."
python3 "$CASE_ROOT/compare_niederer_bueno_orovio_batched.py" \
    "$RUN_ROOT" \
    "$RESULT_ROOT/substep_sweep_metrics.json" \
    $MODES

# ---------------------------------------------------------------------------
# Substep tradeoff plots — merges sub100 from first comparison automatically
# ---------------------------------------------------------------------------
echo ""
echo "Generating substep tradeoff plots ..."
python3 "$CASE_ROOT/plot_substep_sweep.py" \
    "$RESULT_ROOT/substep_sweep_metrics.json" \
    "$COMPARISON_METRICS" \
    "$RESULT_ROOT/plots" \
    $SUBSTEPS

echo ""
echo "Sweep metrics  : $RESULT_ROOT/substep_sweep_metrics.json"
echo "Plots          : $RESULT_ROOT/plots/"
