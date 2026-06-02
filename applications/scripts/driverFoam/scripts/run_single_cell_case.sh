#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: run_single_cell_case.sh --case-dir <path> [--openfoam-bashrc <path>]" >&2
}

CASE_DIR=""
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
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$CASE_DIR" ]]; then
    usage
    exit 1
fi

if [[ -z "$OPENFOAM_BASHRC" || ! -f "$OPENFOAM_BASHRC" ]]; then
    echo "OpenFOAM bashrc not found: $OPENFOAM_BASHRC" >&2
    exit 1
fi

cd "$CASE_DIR"

set +e +u
source "$OPENFOAM_BASHRC"
set -e -u

rm -f log.cardiacFoam
cardiacFoam > log.cardiacFoam 2>&1
