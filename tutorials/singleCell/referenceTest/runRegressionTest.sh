#!/bin/bash

# Run from tutorial root regardless of where this script is called from
cd "${0%/*}/.." || exit 1

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

prepareRegressionCase()
{
    tutorialRoot="$1"
    caseDir="$2"

    rm -rf "${caseDir}"
    mkdir -p "${caseDir}"

    for item in constant system
    do
        if [ -d "${tutorialRoot}/${item}" ]; then
            cp -r "${tutorialRoot}/${item}" "${caseDir}/"
        fi
    done

    for item in Allrun Allclean README.md plotSweep plotVoltage
    do
        if [ -e "${tutorialRoot}/${item}" ]; then
            cp -r "${tutorialRoot}/${item}" "${caseDir}/"
        fi
    done

    for dir in constant system
    do
        [ -d "${tutorialRoot}/referenceTest/${dir}" ] || continue
        find "${tutorialRoot}/referenceTest/${dir}" -type f | while IFS= read -r refFile; do
            target="${caseDir}/${refFile#${tutorialRoot}/referenceTest/}"
            mkdir -p "$(dirname "${target}")"
            echo "Restoring ${target#${caseDir}/} from ${refFile}"
            cp "${refFile}" "${target}"
        done
    done
}

checkSingleCellResults()
{
    caseDir="$1"
    tutorialRoot="$2"
    refFile="${tutorialRoot}/referenceTest/singleCell.reference"
    outDir="${caseDir}/postProcessing"

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

tutorialRoot="$(pwd)"
regressionCaseDir="${tutorialRoot}/regressionTests/case"

prepareRegressionCase "${tutorialRoot}" "${regressionCaseDir}"

(
    cd "${regressionCaseDir}" || exit 1
    runApplication cardiacFoam
) || exit 1

# Compare single-cell trace against regression checkpoints
checkSingleCellResults "${regressionCaseDir}" "${tutorialRoot}" || exit 1

rm -rf "${regressionCaseDir}"
