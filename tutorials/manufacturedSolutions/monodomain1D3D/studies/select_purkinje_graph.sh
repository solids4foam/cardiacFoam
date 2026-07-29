#!/bin/bash

set -euo pipefail

if [[ $# -ne 1 ]]
then
    echo "Usage: $0 nodes003|nodes011|nodes021|nodes041|nodes081|nodes161" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GRAPH_ID="$1"
GRAPH_FILE="$CASE_DIR/constant/purkinjeGraph.${GRAPH_ID}"

if [[ ! -f "$GRAPH_FILE" ]]
then
    echo "Missing graph file: $GRAPH_FILE" >&2
    exit 1
fi

cp "$GRAPH_FILE" "$CASE_DIR/constant/purkinjeGraph"
echo "Selected constant/purkinjeGraph.${GRAPH_ID}"
