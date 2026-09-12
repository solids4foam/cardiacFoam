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
#
# ${BASH_SOURCE[0]} is empty when this script is sourced directly into a
# zsh shell (BASH_SOURCE is a bash-only array; zsh never populates it) -
# silently falling back to $0 in that case (zsh does set $0 to the sourced
# file's own path, unlike bash). Without this fallback, dirname on an empty
# string resolves _thisDir to the caller's cwd, which pushes _repoRoot one
# directory too high (e.g. to $HOME instead of the repo) - the bundled
# modules/solids4foam path then silently fails its existence check and this
# script falls through to picking a DIFFERENT solids4foam tree than the one
# other already-built libraries were compiled against. That mismatch is
# the "two solids4foam trees" ABI hazard: everything compiles and links
# without a single warning, then crashes at runtime.
_thisDir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
_repoRoot="$(cd "$_thisDir/.." && pwd)"
_bundledSolids4Foam="$_repoRoot/modules/solids4foam"
_bundledPhysicsModel="$_repoRoot/modules/physicsModel"

# Fail loudly, not silently, if the above still didn't land on the repo
# root - e.g. some other shell/invocation pattern makes both
# BASH_SOURCE[0] and $0 unreliable. Every path derived from _repoRoot
# below depends on this being correct; a wrong value has no other
# symptom until a build silently links against the wrong solids4foam tree.
if [ ! -d "$_repoRoot/src/ionicModels" ]; then
    echo "FATAL: resolveSolids4Foam.sh could not locate the cardiacFoam repo root." >&2
    echo "  Computed _repoRoot='$_repoRoot' (missing src/ionicModels)." >&2
    echo "  Source this script as 'etc/resolveSolids4Foam.sh' from the repo root." >&2
    return 1 2>/dev/null || exit 1
fi
# solidModel.H only exists in the full solids4foam, never in the lightweight
# physicsModel replacement — use it as a reliable discriminator.
_solids4FoamHeader="src/solids4FoamModels/solidModels/solidModel/solidModel.H"
# A solids4foam tree is "built" once its lnInclude has been generated; test for
# a representative header there to avoid selecting an un-built source tree.
_s4fLnHeader="src/solids4FoamModels/lnInclude/solidModel.H"

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
    if [ -f "$SOLIDS4FOAM_INST_DIR/$_s4fLnHeader" ]
    then
        echo "Using full solids4foam installation (built lnInclude found)."
        echo "SOLIDS4FOAM_INST_DIR=$SOLIDS4FOAM_INST_DIR"
        echo
        export USE_LIGHTWEIGHT_PHYSICSMODEL=0
    elif [ -f "$SOLIDS4FOAM_INST_DIR/$_solids4FoamHeader" ]
    then
        useLightweightPhysicsModel \
            "SOLIDS4FOAM_INST_DIR is set to a solids4foam source tree, but it is not built (missing $_s4fLnHeader)."
    else
        useLightweightPhysicsModel \
            "SOLIDS4FOAM_INST_DIR is set, but it is not a valid solids4foam tree (missing $_solids4FoamHeader)."
    fi
    echo
else
    # SOLIDS4FOAM_INST_DIR not set: auto-discover a *built* solids4foam so it
    # doesn't have to be exported every time. A candidate counts only when its
    # lnInclude is populated (i.e. actually compiled) — this avoids silently
    # selecting an un-built source tree (e.g. a freshly checked-out submodule),
    # which otherwise fails the build with missing headers (solidModel.H, ...).
    #
    # The bundled submodule is tried FIRST, ahead of any out-of-tree checkout.
    # When several built solids4foam trees exist on one machine they are not
    # interchangeable: this repo's submodule carries local changes, so its
    # solidModel has a different memory layout from an external checkout's.
    # Picking external headers while $FOAM_USER_LIBBIN holds a
    # libsolids4FoamModels.dylib built from the submodule (or vice versa)
    # compiles and links without a single warning, then crashes at run time
    # inside a constructor — e.g. sequentialElectroMechanical registering its
    # Ta field via solid().mesh(), where solid() resolves to the wrong offset
    # and the objectRegistry reference is garbage. Defaulting to the repo's own
    # submodule makes the common case self-consistent; set SOLIDS4FOAM_INST_DIR
    # explicitly to override.
    _s4fFound=""
    for _cand in "$_bundledSolids4Foam" "$HOME/solids4foam" "$WM_PROJECT_USER_DIR/solids4foam"
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
