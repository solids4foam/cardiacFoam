#!/bin/bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_case.sh --case-dir <path> [options]

Options:
  --case-dir <path>         Case directory (required)
  --dimension <1D|2D|3D>    Run blockMesh with system/blockMeshDict.<dimension>
  --parallel                Pass "parallel" to ./Allrun
  --touch-case-foam         Touch case.foam after run
  --openfoam-bashrc <path>  OpenFOAM bashrc path
  -h, --help                Show this help
EOF
}

CASE_DIR=""
DIMENSION=""
ALLRUN_ARGS=()
TOUCH_CASE_FOAM=false
DEFAULT_OPENFOAM_BASHRC=""
if [[ -n "${WM_PROJECT_DIR:-}" && -f "${WM_PROJECT_DIR}/etc/bashrc" ]]; then
    DEFAULT_OPENFOAM_BASHRC="${WM_PROJECT_DIR}/etc/bashrc"
elif [[ -f "/Volumes/OpenFOAM-v2412/etc/bashrc" ]]; then
    DEFAULT_OPENFOAM_BASHRC="/Volumes/OpenFOAM-v2412/etc/bashrc"
fi
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-$DEFAULT_OPENFOAM_BASHRC}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --case-dir)
            CASE_DIR="$2"
            shift 2
            ;;
        --dimension)
            DIMENSION="$2"
            shift 2
            ;;
        --parallel)
            ALLRUN_ARGS+=("parallel")
            shift
            ;;
        --touch-case-foam)
            TOUCH_CASE_FOAM=true
            shift
            ;;
        --openfoam-bashrc)
            OPENFOAM_BASHRC="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "$CASE_DIR" ]]; then
    echo "Missing required argument: --case-dir" >&2
    usage >&2
    exit 1
fi

if [[ -z "$OPENFOAM_BASHRC" || ! -f "$OPENFOAM_BASHRC" ]]; then
    echo "OpenFOAM bashrc not found: $OPENFOAM_BASHRC" >&2
    echo "Set OPENFOAM_BASHRC, pass --openfoam-bashrc, or source OpenFOAM before running." >&2
    exit 1
fi

cd "$CASE_DIR"
# OpenFOAM bashrc relies on undefined vars and non-zero helper branches.
set +e +u
source "$OPENFOAM_BASHRC"
set -e -u

echo "Cleaning case"
if ! ./Allclean; then
    echo "Warning: Allclean returned non-zero, continuing." >&2
fi

if [[ -n "$DIMENSION" ]]; then
    DICT="system/blockMeshDict.$DIMENSION"
    echo "Running blockMesh with: $DICT (output -> log.blockMesh)"
    blockMesh -dict "$DICT" > log.blockMesh 2>&1
fi

echo "Running case"
if [[ ${#ALLRUN_ARGS[@]} -gt 0 ]]; then
    ./Allrun "${ALLRUN_ARGS[@]}"
else
    ./Allrun
fi

if $TOUCH_CASE_FOAM; then
    echo "Touching case.foam"
    touch case.foam
fi
