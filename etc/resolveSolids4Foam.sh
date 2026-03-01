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

# Force lightweight mode if requested
if [ "${FORCE_LIGHTWEIGHT_PHYSICSMODEL:-0}" = "1" ]; then
    SOLIDS4FOAM_INST_DIR="$_thisDir/../modules/physicsModel"
    export USE_LIGHTWEIGHT_PHYSICSMODEL=1

    wmakeLnInclude $SOLIDS4FOAM_INST_DIR/src/solids4FoamModels

    echo "Using lightweight physicsModel (forced)"
    echo "SOLIDS4FOAM_INST_DIR=$SOLIDS4FOAM_INST_DIR"
    echo

    export SOLIDS4FOAM_INST_DIR
    return 0 2>/dev/null || exit 0
fi

if [ -n "$SOLIDS4FOAM_INST_DIR" ] && [ -d "$SOLIDS4FOAM_INST_DIR" ]
then
    # Do nothing
    export USE_LIGHTWEIGHT_PHYSICSMODEL=0
    echo
elif [ -d "$_thisDir/../modules/solids4foam" ]
then
    echo
    SOLIDS4FOAM_INST_DIR="$_thisDir/../modules/solids4foam"
    export USE_LIGHTWEIGHT_PHYSICSMODEL=0
else
    echo "NOTE: solids4foam not found."
    echo "To us solids4foam, set SOLIDS4FOAM_INST_DIR or initialise submodules:"
    echo "  git submodule update --init --recursive"

    # Use physicsModel from modules
    SOLIDS4FOAM_INST_DIR="$_thisDir/../modules/physicsModel"

    # Build light-weight physics model replacement
    export USE_LIGHTWEIGHT_PHYSICSMODEL=1
    wmake libso $SOLIDS4FOAM_INST_DIR/src/solids4FoamModels
fi

echo "Using SOLIDS4FOAM_INST_DIR=$SOLIDS4FOAM_INST_DIR"
echo

export SOLIDS4FOAM_INST_DIR
