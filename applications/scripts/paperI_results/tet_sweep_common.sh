#!/usr/bin/env bash
# Shared helpers for the tetrahedral-mesh MMS convergence sweeps, sourced by
# each case's setup/run_*.sh. This library installs NO top-level `trap ... EXIT`
# — the calling script owns its traps (see mms_controldict_backup/restore) so a
# case-specific trap can never be silently overwritten by the library.

# Source the OpenFOAM environment only if it is not already active.
mms_env_bootstrap() {
    if [ -z "${WM_PROJECT_DIR:-}" ]; then
        # shellcheck disable=SC1091
        source /Volumes/OpenFOAM-v2412/etc/bashrc
    fi
}

# Effective mesh spacing lc = 1/N, formatted WITHOUT the %.17g float artifact.
# `awk 'BEGIN{printf "%.17g",1.0/n}'` emits e.g. 0.025000000000000001 for N=40,
# which then leaks into the generated .geo; python's repr is the clean 0.025.
mms_lc() {
    python3 -c "print(1.0/$1)"
}

# Per-N time step for the transient MMS cases (the steady eikonal case ignores it).
dt_for_n() {
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        80) echo 0.000140174 ;;
        *) echo "no dt configured for N=$1" >&2; return 1 ;;
    esac
}

# Transactional controlDict edit: back up the dict, then restore on demand.
# Sourced functions must never `exit`; the CALLER installs
# `trap 'mms_controldict_restore' EXIT` so the case dir is left pristine.
_MMS_CONTROLDICT_SRC=""
_MMS_CONTROLDICT_BAK=""
mms_controldict_backup() {
    _MMS_CONTROLDICT_SRC="$1"
    _MMS_CONTROLDICT_BAK="$(mktemp)"
    cp "$_MMS_CONTROLDICT_SRC" "$_MMS_CONTROLDICT_BAK"
}
mms_controldict_restore() {
    if [ -n "$_MMS_CONTROLDICT_BAK" ] && [ -f "$_MMS_CONTROLDICT_BAK" ]; then
        cp "$_MMS_CONTROLDICT_BAK" "$_MMS_CONTROLDICT_SRC"
        rm -f "$_MMS_CONTROLDICT_BAK"
        _MMS_CONTROLDICT_BAK=""
    fi
}
