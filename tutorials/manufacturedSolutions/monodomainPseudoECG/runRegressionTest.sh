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
    if [ -n "${candidate}" ] && grep -q 'Manufactured-solution error summary' "${candidate}"; then
        echo "${candidate}"
        return 0
    fi

    for candidate in postProcessing/*.dat processor*/postProcessing/*.dat
    do
        if [ -s "${candidate}" ] && grep -q 'Manufactured-solution error summary' "${candidate}"; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

findPseudoECGFile()
{
    findFirstMatch 'postProcessing/pseudoECG.dat' \
        || findFirstMatch 'processor*/postProcessing/pseudoECG.dat'
}

findManufacturedPseudoECGSummary()
{
    findFirstMatch 'postProcessing/manufacturedPseudoECGSummary.dat' \
        || findFirstMatch 'processor*/postProcessing/manufacturedPseudoECGSummary.dat'
}

checkPseudoECGResults()
{
    local dataFile
    dataFile="$(findPseudoECGFile)" || {
        echo "FAIL: pseudoECG.dat not found in postProcessing/ or processor*/postProcessing/."
        return 1
    }

    if ! grep -q '^#.*time' "${dataFile}"; then
        echo "FAIL: pseudoECG header missing in ${dataFile}"
        return 1
    fi

    if ! awk '
        $1 !~ /^#/ && NF >= 2 { found = 1; exit 0 }
        END { exit !found }
    ' "${dataFile}"; then
        echo "FAIL: no numeric pseudoECG samples found in ${dataFile}"
        return 1
    fi

    echo "PASS: pseudoECG output detected in ${dataFile}"
    return 0
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

checkManufacturedPseudoECGSummary()
{
    local summaryFile
    summaryFile="$(findManufacturedPseudoECGSummary)" || {
        echo "FAIL: manufacturedPseudoECGSummary.dat not found."
        return 1
    }

    if ! grep -q '^Manufactured pseudo-ECG summary' "${summaryFile}"; then
        echo "FAIL: manufactured pseudo-ECG summary header missing in ${summaryFile}"
        return 1
    fi

    echo "PASS: manufactured pseudo-ECG summary detected in ${summaryFile}"
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
checkPseudoECGResults || exit 1
checkManufacturedPseudoECGSummary || exit 1
