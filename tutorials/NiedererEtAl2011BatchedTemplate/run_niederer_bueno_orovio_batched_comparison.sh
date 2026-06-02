#!/bin/bash
# run_niederer_bueno_orovio_batched_comparison.sh
#
# Runs the NiedererEtAl2011 Niederer benchmark for each ionic model mode
# (cpu, batched_euler, batched_rl, batched_soa) and collects
# wall-clock time for the cardiacFoam solve step only.
#
# Layout produced:
#   comparisonRuns/<mode>/run01/     ← full OpenFOAM case (constant/ system/ log.* postProcessing/)
#   comparisonResults/comparison_metrics.json
#   comparisonResults/plots/
#
# Key env-var overrides:
#   N_RUNS=1              number of repeat runs per mode (default 1; tissue runs are slow)
#   MODES="cpu ..."       space-separated list of modes
#   CF_OPENFOAM_BASHRC    path to OpenFOAM bashrc
#   RUN_ROOT              override comparisonRuns location
#   RESULT_ROOT           override comparisonResults location

set -euo pipefail

CASE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_ROOT="${RUN_ROOT:-$CASE_ROOT/comparisonRuns}"
RESULT_ROOT="${RESULT_ROOT:-$CASE_ROOT/comparisonResults}"
N_RUNS="${N_RUNS:-1}"
MODES="${MODES:-cpu batched_euler batched_rl batched_soa}"
OPENFOAM_BASHRC="${CF_OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
N_SUBDOMAINS="${N_SUBDOMAINS:-6}"

# ---------------------------------------------------------------------------
# Source OpenFOAM environment if available
# ---------------------------------------------------------------------------
if [ -f "$OPENFOAM_BASHRC" ]; then
    # shellcheck disable=SC1090
    set +eu
    source "$OPENFOAM_BASHRC"
    set -eu
fi

# ---------------------------------------------------------------------------
# Map mode name → electroProperties filename
# ---------------------------------------------------------------------------
config_for_mode()
{
    case "$1" in
        cpu)            echo "electroProperties.cpu" ;;
        batched_euler)  echo "electroProperties.batched_euler" ;;
        batched_rl)     echo "electroProperties.batched_rl" ;;
        batched_soa)    echo "electroProperties.batched_soa" ;;
        *)
            echo "Unknown mode '$1'" >&2
            return 1
            ;;
    esac
}

# ---------------------------------------------------------------------------
# Stage a fresh case directory for one run
# ---------------------------------------------------------------------------
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

    # Copy constant/ and system/ from the tutorial root
    cp -R "$CASE_ROOT/constant" "$case_dir/"
    cp -R "$CASE_ROOT/system"   "$case_dir/"

    # Install the selected electroProperties
    cp "$CASE_ROOT/constant/$config_name" "$case_dir/constant/electroProperties"

    echo "Prepared $case_dir (mode=$mode, config=$config_name)"
}

# ---------------------------------------------------------------------------
# Run blockMesh + decomposePar, then time cardiacFoam, then reconstructPar
# and postProcess.  Wall time covers cardiacFoam only.
# ---------------------------------------------------------------------------
run_case()
{
    local mode="$1"
    local run_id="$2"
    local case_dir="$RUN_ROOT/$mode/run$(printf '%02d' "$run_id")"
    local start end elapsed

    (
        cd "$case_dir"

        # --- mesh (not timed) ---
        echo "[${mode}/run${run_id}] Running blockMesh ..."
        blockMesh > log.blockMesh 2>&1

        # --- decompose (not timed) ---
        echo "[${mode}/run${run_id}] Running decomposePar ..."
        decomposePar > log.decomposePar 2>&1

        # --- time cardiacFoam only ---
        echo "[${mode}/run${run_id}] Running cardiacFoam (parallel, N=${N_SUBDOMAINS}) ..."
        start="$(python3 -c 'import time; print(time.perf_counter())')"
        mpirun -np "$N_SUBDOMAINS" cardiacFoam -parallel > log.cardiacFoam 2>&1
        end="$(python3 -c 'import time; print(time.perf_counter())')"
        elapsed="$(python3 -c "print(${end} - ${start})")"
        printf 'real %.9f\n' "$elapsed" > time.txt
        echo "[${mode}/run${run_id}] cardiacFoam wall time: ${elapsed} s"

        # --- reconstruct (not timed) ---
        echo "[${mode}/run${run_id}] Running reconstructPar ..."
        reconstructPar > log.reconstructPar 2>&1

        # --- post-process: extract activationTime at probe locations ---
        echo "[${mode}/run${run_id}] Running postProcess (Niedererpoints) ..."
        postProcess -func Niedererpoints -latestTime > log.postProcess_points 2>&1
        echo "[${mode}/run${run_id}] Running postProcess (Niedererlines) ..."
        postProcess -func Niedererlines  -latestTime > log.postProcess_lines  2>&1
    )
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
mkdir -p "$RUN_ROOT" "$RESULT_ROOT"

for mode in $MODES; do
    rm -rf "$RUN_ROOT/$mode"
    mkdir -p "$RUN_ROOT/$mode"

    for run_id in $(seq 1 "$N_RUNS"); do
        echo ""
        echo "========================================================"
        echo "  Mode: $mode  |  Run: $run_id / $N_RUNS"
        echo "========================================================"
        prepare_case "$mode" "$run_id"
        run_case     "$mode" "$run_id"
    done
done

# ---------------------------------------------------------------------------
# Compare and plot
# ---------------------------------------------------------------------------
echo ""
echo "All runs complete. Running comparison ..."
python3 "$CASE_ROOT/compare_niederer_bueno_orovio_batched.py" \
    "$RUN_ROOT" \
    "$RESULT_ROOT/comparison_metrics.json" \
    $MODES

echo ""
echo "Generating plots ..."
python3 "$CASE_ROOT/plot_niederer_bueno_orovio_batched.py" \
    "$RUN_ROOT" \
    "$RESULT_ROOT/plots" \
    $MODES

echo ""
echo "Comparison metrics : $RESULT_ROOT/comparison_metrics.json"
echo "Plots              : $RESULT_ROOT/plots/"
