#!/bin/bash

scriptDir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
regressionRoot="${scriptDir}"

# Output results and temporary cases inside the regression folder
regressionCaseDir="${regressionRoot}/case"

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

prepareRegressionCase()
{
    rRoot="$1"
    cDir="$2"

    rm -rf "${cDir}"
    mkdir -p "${cDir}"
    
    # Copy all folders from regressionRoot into the caseDir
    for item in constant system
    do
        if [ -d "${rRoot}/${item}" ]; then
            cp -r "${rRoot}/${item}" "${cDir}/"
        fi
    done

    # Ensure log folders exist
    mkdir -p "${cDir}/postProcessing"
}

checkSingleCellResults()
{
    cDir="$1"
    rRoot="$2"
    refFile="${rRoot}/singleCell.reference"
    outDir="${cDir}/postProcessing"

    if [ ! -f "${refFile}" ]; then
        echo "Reference file not found: ${refFile}"
        return 1
    fi

    if [ ! -d "${outDir}" ]; then
        echo "Post-processing directory not found: ${outDir}"
        return 1
    fi

    nFail=0
    nCheck=0

    while read -r fileName time expected tolerance; do
        if [ -z "${fileName}" ] || [[ "${fileName}" == \#* ]]; then
            continue
        fi

        dataFile="${outDir}/${fileName}"
        if [ ! -f "${dataFile}" ]; then
            echo "FAIL: missing output file ${dataFile}"
            nFail=$((nFail + 1))
            nCheck=$((nCheck + 1))
            continue
        fi

        actual="$(
            awk -v target="${time}" '
                BEGIN { bestDiff = 1e99; found = 0; actual = 0.0; }
                NR > 1 && NF >= 2 {
                    d = $1 - target;
                    if (d < 0) d = -d;
                    if (d < bestDiff) {
                        bestDiff = d;
                        actual = $2;
                        found = 1;
                    }
                }
                END {
                    if (found && bestDiff <= 1e-9) {
                        print actual;
                        exit 0;
                    }
                    exit 1;
                }
            ' "${dataFile}"
        )" || true

        nCheck=$((nCheck + 1))

        if [ -z "${actual}" ]; then
            echo "FAIL: ${fileName} at t=${time} not found in ${dataFile}"
            nFail=$((nFail + 1))
            continue
        fi

        if awk -v a="${actual}" -v e="${expected}" -v t="${tolerance}" '
            BEGIN {
                d = a - e;
                if (d < 0) d = -d;
                exit !(d <= t);
            }
        '; then
            echo "PASS: ${fileName} t=${time} actual=${actual} expected=${expected} tol=${tolerance}"
        else
            echo "FAIL: ${fileName} t=${time} actual=${actual} expected=${expected} tol=${tolerance}"
            nFail=$((nFail + 1))
        fi
    done < "${refFile}"

    echo "singleCell comparison: ${nCheck} checks, ${nFail} failures"

    if [ "${nFail}" -ne 0 ]; then
        return 1
    fi

    return 0
}

prepareRegressionCase "${regressionRoot}" "${regressionCaseDir}"

(
    cd "${regressionCaseDir}" || exit 1
    runApplication cardiacFoam
) || exit 1

# Compare single-cell trace against regression checkpoints
checkSingleCellResults "${regressionCaseDir}" "${regressionRoot}" || exit 1
