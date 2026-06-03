#!/bin/bash

# Resolve and export SOLIDS4FOAM_INST_DIR, preferring a system solids4foam
# install and falling back to the bundled submodule if necessary

# Guard against repeated sourcing (e.g. via AllwmakeParseArguments)
if [ -n "$_SOLIDS4FOAM_RESOLVED" ]; then
    return 0 2>/dev/null || exit 0
fi
export _SOLIDS4FOAM_RESOLVED=1

# Directory containing this script (absolute). The repo root is its parent;
# resolve it with cd+pwd so the exported paths never contain a literal '..'
# (a '..' segment makes $(abspath ...) in Make/options collapse the repo dir
# name and mismatch CARDIACFOAM_LIGHTWEIGHT_ROOT, breaking the build).
_thisDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_repoRoot="$(cd "$_thisDir/.." && pwd)"
_bundledSolids4Foam="$_repoRoot/modules/solids4foam"
_bundledPhysicsModel="$_repoRoot/modules/physicsModel"
_solids4FoamHeader="src/solids4FoamModels/physicsModel/physicsModel.H"
# A solids4foam tree is "built" once its lnInclude has been generated; test for
# a representative header there to avoid selecting an un-built source tree.
_s4fLnHeader="src/solids4FoamModels/lnInclude/physicsModel.H"

useLightweightPhysicsModel()
{
    local reason="$1"

    SOLIDS4FOAM_INST_DIR="$_bundledPhysicsModel"
    export USE_LIGHTWEIGHT_PHYSICSMODEL=1

    echo "$reason"
    echo "Using lightweight electrophysiology compilation."
    echo "SOLIDS4FOAM_INST_DIR=$SOLIDS4FOAM_INST_DIR"
    echo

    wmake libso "$SOLIDS4FOAM_INST_DIR/src/solids4FoamModels"

    export SOLIDS4FOAM_INST_DIR
}

# Force lightweight mode if requested
if [ "${FORCE_LIGHTWEIGHT_PHYSICSMODEL:-0}" = "1" ]
then
    useLightweightPhysicsModel "FORCE_LIGHTWEIGHT_PHYSICSMODEL=1: skipping solids4foam."
    return 0 2>/dev/null || exit 0
fi

if [ -n "$SOLIDS4FOAM_INST_DIR" ] && [ -d "$SOLIDS4FOAM_INST_DIR" ]
then
    if [ -f "$SOLIDS4FOAM_INST_DIR/$_solids4FoamHeader" ]
    then
        echo "Using compiled solids4foam installation."
        echo "SOLIDS4FOAM_INST_DIR=$SOLIDS4FOAM_INST_DIR"
        echo
        export USE_LIGHTWEIGHT_PHYSICSMODEL=0
    else
        useLightweightPhysicsModel \
            "SOLIDS4FOAM_INST_DIR is set, but solids4foam is not compiled or is missing $_solids4FoamHeader."
    fi
    echo
else
    # SOLIDS4FOAM_INST_DIR not set: auto-discover a *built* solids4foam so it
    # doesn't have to be exported every time. A candidate counts only when its
    # lnInclude is populated (i.e. actually compiled) — this avoids silently
    # selecting an un-built source tree (e.g. a freshly checked-out submodule),
    # which otherwise fails the build with missing headers (solidModel.H, ...).
    _s4fFound=""
    for _cand in "$HOME/solids4foam" "$WM_PROJECT_USER_DIR/solids4foam" "$_bundledSolids4Foam"
    do
        if [ -n "$_cand" ] && [ -f "$_cand/$_s4fLnHeader" ]
        then
            _s4fFound="$_cand"
            break
        fi
    done

    if [ -n "$_s4fFound" ]
    then
        echo
        echo "Auto-detected built solids4foam at: $_s4fFound"
        SOLIDS4FOAM_INST_DIR="$_s4fFound"
        export USE_LIGHTWEIGHT_PHYSICSMODEL=0
    elif [ -f "$_bundledSolids4Foam/$_solids4FoamHeader" ]
    then
        useLightweightPhysicsModel \
            "Bundled solids4foam is present but not built (empty lnInclude). Build it with 'cd modules/solids4foam && ./Allwmake', or set SOLIDS4FOAM_INST_DIR to a built install. Using lightweight for now."
    else
        echo "NOTE: solids4foam not compiled or not initialized."
        echo "To use solids4foam, set SOLIDS4FOAM_INST_DIR or initialise submodules:"
        echo "  git submodule update --init --recursive"
        useLightweightPhysicsModel "Falling back to lightweight electrophysiology compilation."
    fi
fi

echo "Using SOLIDS4FOAM_INST_DIR=$SOLIDS4FOAM_INST_DIR"
echo

export SOLIDS4FOAM_INST_DIR
