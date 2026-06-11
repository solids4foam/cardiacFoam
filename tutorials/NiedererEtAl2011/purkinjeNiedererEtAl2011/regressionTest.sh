#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ============================================================
# Niederer Purkinje slab regression test
#
# Phase 1 — monodomain 1-D Purkinje + 3-D reaction-diffusion
#   Graph utility run + coupled run with reference-value checks.
#
# Phase 2 — eikonal 1-D Purkinje + 3-D steady-state eikonal
#   Run to completion and check quantitative activationTime values
#   at PVJ nodes against eikonalSlab.reference.
# ============================================================

REF_FILE="purkinjeSlab.reference"
END_TIME=0.02
DT=1e-5
GRAPH_STEPS=2000
BLOCKMESH_LOGFILE="log.blockMesh"
GRAPH_LOGFILE="log.runPurkinjeGraph"
ALLRUN_LOGFILE="log.Allrun"
EIKONAL_LOGFILE="log.Allrun.eikonal"

if [[ "$(uname -s)" == "Darwin" && -n "${WM_PROJECT_DIR:-}" && -n "${WM_OPTIONS:-}" ]]; then
    openfoamLibDir="${WM_PROJECT_DIR}/platforms/${WM_OPTIONS}/lib"
    if [[ -d "${openfoamLibDir}" ]]; then
        export DYLD_LIBRARY_PATH="${openfoamLibDir}:${DYLD_LIBRARY_PATH:-}"
    fi
fi

# ----------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------

absDiff()
{
    awk -v a="$1" -v e="$2" '
        BEGIN {
            d = a - e;
            if (d < 0) d = -d;
            print d;
        }
    '
}

checkWithinTolerance()
{
    local label="$1"
    local actual="$2"
    local expected="$3"
    local tolerance="$4"
    local diffAbs

    diffAbs="$(absDiff "${actual}" "${expected}")"

    if awk -v d="${diffAbs}" -v t="${tolerance}" 'BEGIN {exit !(d <= t)}'; then
        printf "PASS: %s actual=%.8g expected=%.8g difference=%.3g tolerance=%.3g\n" \
            "${label}" "${actual}" "${expected}" "${diffAbs}" "${tolerance}"
        return 0
    fi

    printf "FAIL: %s actual=%.8g expected=%.8g difference=%.3g tolerance=%.3g\n" \
        "${label}" "${actual}" "${expected}" "${diffAbs}" "${tolerance}"
    return 1
}

graphDictValue()
{
    local key="$1"

    awk -v key="${key}" '
        $1 == key {
            while (getline line) {
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", line);
                if (line == "") {
                    continue;
                }
                sub(/;.*/, "", line);
                sub(/\(.*/, "", line);
                print line;
                exit;
            }
        }
    ' constant/purkinjeGraph
}

latestVtkFile()
{
    find postProcessing/purkinjeNetworkVTK -name 'purkinjeNetwork_*.vtk' | sort | tail -n 1
}

extractVtkNodeScalar()
{
    local vtkFile="$1"
    local scalarName="$2"
    local nodeIndex="$3"

    awk -v scalarName="${scalarName}" -v nodeIndex="${nodeIndex}" '
        $1 == "SCALARS" && $2 == scalarName {
            inScalar = 1;
            getline;
            count = 0;
            next;
        }
        inScalar && $1 == "SCALARS" {
            inScalar = 0;
        }
        inScalar && NF > 0 {
            if (count == nodeIndex) {
                print $1;
                exit;
            }
            count++;
        }
    ' "${vtkFile}"
}

extractFinalPvjValue()
{
    local pvjName="$1"
    local column

    case "${pvjName}" in
        pvj0) column=2 ;;
        pvj1) column=3 ;;
        pvj2) column=4 ;;
        pvj3) column=5 ;;
        pvj4) column=6 ;;
        pvj5) column=7 ;;
        pvj6) column=8 ;;
        pvj7) column=9 ;;
        *) return 1 ;;
    esac

    awk -v column="${column}" '$1 !~ /^#/ {value = $column} END {print value}' \
        postProcessing/purkinjeNetwork.dat
}

# Eikonal case: Purkinje activation times written to purkinjeNetwork.dat
# columns: time, node0_activationTime, node1_activationTime, ...
# nodeN is at column N+2.
extractEikonalPurkinjeAT()
{
    local nodeKey="$1"
    local nodeNum="${nodeKey#node}"
    local column=$((nodeNum + 2))

    awk -v col="${column}" '$1 !~ /^#/ {print $col; exit}' \
        postProcessing/purkinjeNetwork.dat
}

extractReferenceValue()
{
    local kind="$1"
    local key="$2"
    local metric="$3"
    local graphUtilityVtk="$4"
    local coupledVtk="$5"

    case "${kind}" in
        topology)
            graphDictValue "${key}"
            ;;
        graphUtilityVm)
            extractVtkNodeScalar "${graphUtilityVtk}" Vm_V "${key#node}"
            ;;
        coupledVm)
            extractVtkNodeScalar "${coupledVtk}" Vm_V "${key#node}"
            ;;
        coupledPVJ)
            extractFinalPvjValue "${key}"
            ;;
        eikonalPurkinjeAT)
            extractEikonalPurkinjeAT "${key}"
            ;;
        *)
            return 1
            ;;
    esac
}

checkReferenceValues()
{
    local graphUtilityVtk="$1"
    local coupledVtk="$2"
    local failures=0
    local checks=0
    local kind key metric expected tolerance actual

    if [[ ! -f "${REF_FILE}" ]]; then
        echo "FAIL: reference file not found: ${REF_FILE}"
        return 1
    fi

    while IFS=' ' read -r kind key metric expected tolerance; do
        if [[ -z "${kind}" || "${kind}" == \#* ]]; then
            continue
        fi

        actual="$(extractReferenceValue "${kind}" "${key}" "${metric}" "${graphUtilityVtk}" "${coupledVtk}")"
        checks=$((checks + 1))

        if [[ -z "${actual}" ]]; then
            echo "FAIL: could not extract ${kind} ${key} ${metric}"
            failures=$((failures + 1))
            continue
        fi

        checkWithinTolerance "${kind} ${key} ${metric}" "${actual}" "${expected}" "${tolerance}" \
            || failures=$((failures + 1))
    done < "${REF_FILE}"

    echo "Purkinje slab reference comparison: ${checks} checks, ${failures} failures"
    return "${failures}"
}

# ================================================================
# Phase 1 — monodomain 1-D Purkinje + 3-D reaction-diffusion
# ================================================================

echo "============================================================"
echo "Phase 1: monodomain 1-D Purkinje + 3-D reaction-diffusion"
echo "Graph utility and coupled PVJ slab run to t=${END_TIME}"
echo "============================================================"
echo

./Allclean > /dev/null 2>&1 || true
foamDictionary system/controlDict -entry endTime -set "${END_TIME}" > /dev/null 2>&1
cp constant/electroProperties.monodomain constant/electroProperties

blockMesh > "${BLOCKMESH_LOGFILE}" 2>&1
runPurkinjeGraph -case . -conductionDomain purkinjeNetwork -nSteps "${GRAPH_STEPS}" -deltaT "${DT}" \
    > "${GRAPH_LOGFILE}" 2>&1
graphUtilityVtk="$(latestVtkFile)"
cp "${graphUtilityVtk}" graphUtilityFinal.vtk
graphUtilityVtk="graphUtilityFinal.vtk"

echo "PASS: graph-only Purkinje utility wrote ${graphUtilityVtk}"

./Allclean > /dev/null 2>&1 || true
foamDictionary system/controlDict -entry endTime -set "${END_TIME}" > /dev/null 2>&1

./Allrun parallel > "${ALLRUN_LOGFILE}" 2>&1
coupledVtk="$(latestVtkFile)"

if [[ ! -s postProcessing/purkinjeNetwork.dat ]]; then
    echo "FAIL: coupled run did not write postProcessing/purkinjeNetwork.dat"
    exit 1
fi

echo "PASS: coupled Purkinje-slab run wrote ${coupledVtk}"
checkReferenceValues "${graphUtilityVtk}" "${coupledVtk}"

echo
echo "Phase 1 PASSED"

# ================================================================
# Phase 2 — eikonal 1-D Purkinje + 3-D steady-state eikonal
# ================================================================

echo
echo "============================================================"
echo "Phase 2: eikonal 1-D Purkinje + 3-D steady-state eikonal"
echo "Quantitative checks: PVJ activationTime vs eikonalSlab.reference"
echo "============================================================"
echo

./Allclean > /dev/null 2>&1 || true

./Allrun solver=eikonal > "${EIKONAL_LOGFILE}" 2>&1

# runApplication inside Allrun redirects solver output to log.cardiacFoam,
# not to the wrapper log — check the solver log for FatalError and End.
EIKONAL_SOLVER_LOG="log.cardiacFoam"
if grep -q "FatalError" "${EIKONAL_SOLVER_LOG:-/dev/null}"; then
    echo "FAIL: eikonal run produced a FatalError"
    echo "--- last 20 lines of ${EIKONAL_SOLVER_LOG} ---"
    tail -20 "${EIKONAL_SOLVER_LOG}"
    exit 1
fi

if ! grep -q "^End" "${EIKONAL_SOLVER_LOG:-/dev/null}"; then
    echo "FAIL: eikonal run did not reach normal End"
    echo "--- last 20 lines of ${EIKONAL_LOGFILE} ---"
    tail -20 "${EIKONAL_LOGFILE}"
    exit 1
fi

# Eikonal solver uses applyModelTimeControls → endTime = deltaT = 1.0
# so results land in 1/
if [[ ! -s 1/activationTime ]]; then
    echo "FAIL: 1/activationTime not written by eikonal run"
    exit 1
fi

echo "PASS: eikonal run completed and wrote 1/activationTime"

# Quantitative checks: Purkinje Dijkstra activation times at key nodes.
# These are deterministic (Dijkstra) and bitwise-reproducible.
EIKONAL_REF_FILE="eikonalSlab.reference"
if [[ ! -f "${EIKONAL_REF_FILE}" ]]; then
    echo "FAIL: eikonal reference file not found: ${EIKONAL_REF_FILE}"
    exit 1
fi

eikonalFailures=0
eikonalChecks=0
while IFS=' ' read -r kind key metric expected tolerance; do
    if [[ -z "${kind}" || "${kind}" == \#* ]]; then
        continue
    fi
    actual="$(extractReferenceValue "${kind}" "${key}" "${metric}" "" "")"
    eikonalChecks=$((eikonalChecks + 1))
    if [[ -z "${actual}" ]]; then
        echo "FAIL: could not extract ${kind} ${key} ${metric}"
        eikonalFailures=$((eikonalFailures + 1))
        continue
    fi
    checkWithinTolerance "${kind} ${key} ${metric}" "${actual}" "${expected}" "${tolerance}" \
        || eikonalFailures=$((eikonalFailures + 1))
done < "${EIKONAL_REF_FILE}"

echo "Eikonal slab reference comparison: ${eikonalChecks} checks, ${eikonalFailures} failures"
if (( eikonalFailures > 0 )); then
    exit 1
fi

echo
echo "Phase 2 PASSED"

echo
echo "============================================================"
echo "Regression test PASSED"
echo "============================================================"
exit 0
