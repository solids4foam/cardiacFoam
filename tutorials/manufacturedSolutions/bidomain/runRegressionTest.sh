#!/bin/bash

# Source OpenFOAM run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

findFirstMatch()
{
    local pattern="$1"

    for candidate in ${pattern}
    do
        if [ -s "${candidate}" ]; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

findManufacturedErrorFile()
{
    local candidate

    candidate="$(findFirstMatch 'postProcessing/*.dat')" || true
    if [ -n "${candidate}" ] && grep -Eq '(^| )([Bb]idomain )?[Mm]anufactured-solution error summary' "${candidate}"; then
        echo "${candidate}"
        return 0
    fi

    for candidate in postProcessing/*.dat processor*/postProcessing/*.dat
    do
        if [ -s "${candidate}" ] && grep -Eq '(^| )([Bb]idomain )?[Mm]anufactured-solution error summary' "${candidate}"; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

checkManufacturedErrorSummary()
{
    local errorFile
    errorFile="$(findManufacturedErrorFile)" || {
        echo "FAIL: manufactured error summary file not found."
        return 1
    }

    if ! grep -q 'Field     L1-error' "${errorFile}"; then
        echo "FAIL: manufactured error table missing in ${errorFile}"
        return 1
    fi

    echo "PASS: manufactured error summary detected in ${errorFile}"
    return 0
}

parallelRun=false
for arg in "$@"; do
    case "$arg" in
        parallel)
            parallelRun=true
            ;;
    esac
done

if $parallelRun; then
    runApplication decomposePar
    runParallel cardiacFoam
    runApplication reconstructPar
else
    runApplication cardiacFoam
fi

checkManufacturedErrorSummary || exit 1
