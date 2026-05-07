#!/bin/bash

# Resolve and export SOLIDS4FOAM_INST_DIR, preferring a system solids4foam
# install and falling back to the bundled submodule if necessary

# Guard against repeated sourcing (e.g. via AllwmakeParseArguments)
if [ -n "$_SOLIDS4FOAM_RESOLVED" ]; then
    return 0 2>/dev/null || exit 0
fi
export _SOLIDS4FOAM_RESOLVED=1

# Directory containing this script
_thisDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_bundledSolids4Foam="$_thisDir/../modules/solids4foam"
_bundledPhysicsModel="$_thisDir/../modules/physicsModel"
_solids4FoamHeader="src/solids4FoamModels/physicsModel/physicsModel.H"

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
elif [ -f "$_bundledSolids4Foam/$_solids4FoamHeader" ]
then
    echo
    echo "Using bundled compiled solids4foam."
    SOLIDS4FOAM_INST_DIR="$_bundledSolids4Foam"
    export USE_LIGHTWEIGHT_PHYSICSMODEL=0
else
    echo "NOTE: solids4foam not compiled or not initialized."
    echo "To use solids4foam, set SOLIDS4FOAM_INST_DIR or initialise submodules:"
    echo "  git submodule update --init --recursive"
    useLightweightPhysicsModel "Falling back to lightweight electrophysiology compilation."
fi

echo "Using SOLIDS4FOAM_INST_DIR=$SOLIDS4FOAM_INST_DIR"
echo

export SOLIDS4FOAM_INST_DIR
