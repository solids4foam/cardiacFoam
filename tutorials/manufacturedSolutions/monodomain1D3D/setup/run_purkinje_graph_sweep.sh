#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
OUTPUT_DIR="$CASE_DIR/outputs/1dGraphConvergence"

if [[ -n "${WM_PROJECT_DIR:-}" && -f "$WM_PROJECT_DIR/etc/bashrc" ]]
then
    # Re-source inside the script so macOS dyld paths survive script launch.
    # shellcheck source=/dev/null
    set +eu
    . "$WM_PROJECT_DIR/etc/bashrc"
    set -eu
fi

N_STEPS="${N_STEPS:-1427}"
DELTA_T="${DELTA_T:-0.000140174}"

if [[ "$#" -gt 0 ]]
then
    GRAPH_IDS=("$@")
else
    GRAPH_IDS=(nodes003 nodes011 nodes021 nodes041 nodes081 nodes161)
fi

if [[ ! -d "$CASE_DIR/constant/polyMesh" ]]
then
    (cd "$CASE_DIR" && blockMesh -dict system/blockMeshDict.3D > log.blockMesh)
fi

mkdir -p "$OUTPUT_DIR"

for graph_id in "${GRAPH_IDS[@]}"
do
    echo "Running Purkinje graph sweep case: ${graph_id}"

    "$SCRIPT_DIR/select_purkinje_graph.sh" "$graph_id" >/dev/null

    case_output="$OUTPUT_DIR/$graph_id"
    mkdir -p "$case_output"

    (
        cd "$CASE_DIR"
        rm -f postProcessing/graph_*_nodes.dat
        endTime=$(awk -v n="$N_STEPS" -v dt="$DELTA_T" 'BEGIN {print n * dt}')
        foamDictionary system/controlDict -entry deltaT -set "$DELTA_T" > /dev/null 2>&1
        foamDictionary system/controlDict -entry endTime -set "$endTime" > /dev/null 2>&1
        runPurkinjeGraph -case . > "$case_output/log.runPurkinjeGraph"
    )

    cp "$CASE_DIR/postProcessing/purkinjeNetwork.dat" \
        "$case_output/purkinjeNetwork.dat"

    cp "$CASE_DIR"/postProcessing/graph_*_nodes.dat "$case_output/"

    if [[ -d "$CASE_DIR/postProcessing/purkinjeNetworkVTK" ]]
    then
        mkdir -p "$case_output/purkinjeNetworkVTK"
        cp "$CASE_DIR"/postProcessing/purkinjeNetworkVTK/* \
            "$case_output/purkinjeNetworkVTK/"
    fi
done

"$SCRIPT_DIR/post_processing_purkinje_graph.py" --output-dir "$OUTPUT_DIR"

echo "Wrote graph-only sweep outputs to $OUTPUT_DIR"
