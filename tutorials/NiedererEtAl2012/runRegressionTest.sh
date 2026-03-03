#!/bin/bash

# Source OpenFOAM run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

checkSmokeCheckResults()
{
    caseDir="$(pwd)"
    refFile="${caseDir}/system/smokeCheck.reference"
    postRoot="${caseDir}/postProcessing/smokeCheck"
    postDir="${postRoot}/0"

    if [ ! -f "${refFile}" ]; then
        echo "Reference file not found: ${refFile}"
        return 1
    fi

    if [ ! -d "${postDir}" ]; then
        latestDir=""
        latestVal=""

        if [ -d "${postRoot}" ]; then
            for dir in "${postRoot}"/*; do
                if [ ! -d "${dir}" ]; then
                    continue
                fi

                base="$(basename "${dir}")"
                if [[ ! "${base}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
                    continue
                fi

                if [ -z "${latestDir}" ] || awk -v a="${base}" -v b="${latestVal}" 'BEGIN { exit !(a > b) }'; then
                    latestDir="${dir}"
                    latestVal="${base}"
                fi
            done
        fi

        if [ -n "${latestDir}" ]; then
            postDir="${latestDir}"
            echo "Using smokeCheck output directory: ${postDir}"
        else
            echo "Post-processing directory not found: ${postDir}"
            echo "No numeric smokeCheck output folders found under: ${postRoot}"
            return 1
        fi
    fi

    nFail=0
    nCheck=0

    while read -r field time expected tolerance; do
        if [ -z "${field}" ] || [[ "${field}" == \#* ]]; then
            continue
        fi

        dataFile="${postDir}/${field}"
        if [ ! -f "${dataFile}" ]; then
            echo "FAIL: missing output file ${dataFile}"
            nFail=$((nFail + 1))
            nCheck=$((nCheck + 1))
            continue
        fi

        actual="$(
            awk -v target="${time}" '
                BEGIN { bestDiff = 1e99; found = 0; actual = 0.0; }
                $1 !~ /^#/ && NF >= 2 {
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
            echo "FAIL: ${field} at t=${time} not found in ${dataFile}"
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
            echo "PASS: ${field} t=${time} actual=${actual} expected=${expected} tol=${tolerance}"
        else
            echo "FAIL: ${field} t=${time} actual=${actual} expected=${expected} tol=${tolerance}"
            nFail=$((nFail + 1))
        fi
    done < "${refFile}"

    echo "smokeCheck comparison: ${nCheck} checks, ${nFail} failures"

    if [ "${nFail}" -ne 0 ]; then
        return 1
    fi

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

runApplication blockMesh

if $parallelRun; then
    runApplication decomposePar
    runParallel cardiacFoam
    runApplication reconstructPar
else
    runApplication cardiacFoam
fi

# Regression protocol for this case: smokeCheck only.
runApplication -o postProcess -func smokeCheck
checkSmokeCheckResults || exit 1
