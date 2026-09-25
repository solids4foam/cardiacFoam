#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# Shared regression helpers (tutorials/regressionFunctions)
helperDir="$(cd "${BASH_SOURCE[0]%/*}" && pwd)"
until [[ -f "${helperDir}/regressionFunctions" || "${helperDir}" == / ]]
do
    helperDir="$(dirname "${helperDir}")"
done
. "${helperDir}/regressionFunctions"

# OpenFOAM run functions; also restore DYLD_LIBRARY_PATH on macOS, where
# System Integrity Protection clears it for child shells. RunFunctions
# reads unset variables, so relax 'set -u' while sourcing it.
set +u
. "${WM_PROJECT_DIR:?}/bin/tools/RunFunctions"
set -u

# ============================================================
# Configuration checks
#
# Runs every check listed in ./checks: a short cardiacFoam run of one
# dictionary configuration on a small base case. A "pass" check must run to
# completion and write exactly the files listed in
# regression/manifests/<name>; a "fail:<text>" check must stop with a FOAM
# FATAL error whose message contains <text>. Then ./comparisons checks that
# options change (or leave unchanged) the results they should.
#
# Usage: regression/regressionTest.sh [--update-manifests] [check ...]
#
#   check               run only the named checks (default: all)
#   --update-manifests  rewrite the manifests of the passing checks that ran
# ============================================================

caseDir="$(cd "${BASH_SOURCE[0]%/*}/.." && pwd)"
manifestDir="${caseDir}/regression/manifests"
runRoot="${caseDir}/runs"

updateManifests=0
selected=()
for arg in "$@"; do
    case "${arg}" in
        --update-manifests) updateManifests=1 ;;
        -h|--help)
            sed -n '/^# Usage:/,/^# ====/p' "${BASH_SOURCE[0]}" | sed '$d; s/^# \{0,1\}//'
            exit 0
            ;;
        -*) echo "Unknown option: ${arg}"; exit 2 ;;
        *) selected+=("${arg}") ;;
    esac
done

# Files a run wrote, relative to the run directory: every time directory
# except the initial one, postProcessing/, and constant/*.withDefaultValues
writtenFiles()
{
    (
        cd "$1"
        {
            find . -mindepth 1 -maxdepth 1 -type d -name '[0-9]*' ! -name 0 \
                -exec find {} -type f \; 2>/dev/null
            [[ -d postProcessing ]] && find postProcessing -type f
            find constant -maxdepth 2 -name '*.withDefaultValues' -type f
        } | sed 's|^\./||' | LC_ALL=C sort
    )
}

# True when no checks were named on the command line, or <name> was
isSelected()
{
    local s
    (( ${#selected[@]} == 0 )) && return 0
    for s in "${selected[@]}"; do
        [[ "${s}" == "$1" ]] && return 0
    done
    return 1
}

# Apply "path=value" (set) and "path!" (remove) edits to electroProperties
applyEdits()
{
    local edit
    for edit in "$@"; do
        if [[ "${edit}" == *'!' ]]; then
            foamDictionary constant/electroProperties \
                -entry "${edit%!}" -remove > /dev/null
        else
            foamDictionary constant/electroProperties \
                -entry "${edit%%=*}" -set "${edit#*=}" > /dev/null
        fi
    done
}

# runCheck <name> <base> <variant> <expect> [edit ...]
runCheck()
{
    local name="$1" base="$2" variant="$3" expect="$4"
    shift 4
    local runDir="${runRoot}/${name}"
    local rc=0

    rm -rf "${runDir}"
    mkdir -p "${runRoot}"
    cp -a "${caseDir}/base/${base}" "${runDir}"
    cp -a "${caseDir}/variants/${variant}/." "${runDir}/"

    (
        cd "${runDir}"
        applyEdits "$@"
        ./Allprepare > log.Allprepare 2>&1
    ) || {
        echo "FAIL: ${name}: case preparation failed"
        regressionDumpLog "${runDir}/log.Allprepare"
        return 1
    }

    (cd "${runDir}" && cardiacFoam > log.cardiacFoam 2>&1) || rc=$?

    if [[ "${expect}" == pass ]]; then
        if (( rc != 0 )); then
            echo "FAIL: ${name}: cardiacFoam exited with status ${rc}"
            regressionDumpLog "${runDir}/log.cardiacFoam"
            return 1
        fi
        (cd "${runDir}" && checkSolverLogs > /dev/null) || {
            echo "FAIL: ${name}: run did not complete cleanly"
            regressionDumpLog "${runDir}/log.cardiacFoam"
            return 1
        }

        local manifest="${manifestDir}/${name}"
        if (( updateManifests )); then
            mkdir -p "${manifestDir}"
            writtenFiles "${runDir}" > "${manifest}"
            echo "PASS: ${name} (manifest updated)"
            return 0
        fi
        if [[ ! -f "${manifest}" ]]; then
            echo "FAIL: ${name}: no manifest ${manifest#"${caseDir}"/}"
            return 1
        fi
        if ! diff -u "${manifest}" <(writtenFiles "${runDir}") \
            > "${runDir}/manifest.diff"; then
            echo "FAIL: ${name}: written files differ from the manifest"
            cat "${runDir}/manifest.diff"
            return 1
        fi
        echo "PASS: ${name}"
        return 0
    fi

    local text="${expect#fail:}"
    text="${text//\~/ }"
    if (( rc == 0 )); then
        echo "FAIL: ${name}: cardiacFoam succeeded, expected a fatal error"
        return 1
    fi
    if ! grep -q 'FOAM FATAL' "${runDir}/log.cardiacFoam"; then
        echo "FAIL: ${name}: cardiacFoam failed without a FOAM FATAL error"
        regressionDumpLog "${runDir}/log.cardiacFoam"
        return 1
    fi
    if ! sed -n '/FOAM FATAL/,$p' "${runDir}/log.cardiacFoam" \
        | grep -qF -- "${text}"; then
        echo "FAIL: ${name}: fatal error does not mention '${text}'"
        regressionDumpLog "${runDir}/log.cardiacFoam"
        return 1
    fi
    echo "PASS: ${name} (fails as expected, mentions '${text}')"
    return 0
}

echo "============================================================"
echo "Configuration checks"
echo "============================================================"

checks=0
failures=0
failed=()

while IFS=$' \t' read -r name base variant expect edits; do
    [[ -z "${name}" || "${name}" == \#* ]] && continue

    isSelected "${name}" || continue

    checks=$((checks + 1))
    IFS=' ' read -r -a editList <<< "${edits:-}"
    if ! runCheck "${name}" "${base}" "${variant}" "${expect}" \
        "${editList[@]+"${editList[@]}"}"; then
        failures=$((failures + 1))
        failed+=("${name}")
    fi
done < "${caseDir}/checks"

# Comparisons between the outputs of two checks that both ran (see
# ./comparisons): "same" requires identical files, "differ" different ones
while IFS=$' \t' read -r relation nameA nameB file; do
    [[ -z "${relation}" || "${relation}" == \#* ]] && continue
    [[ -d "${runRoot}/${nameA}" && -d "${runRoot}/${nameB}" ]] || continue
    isSelected "${nameA}" && isSelected "${nameB}" || continue

    fileA="${runRoot}/${nameA}/${file}"
    fileB="${runRoot}/${nameB}/${file}"
    checks=$((checks + 1))

    if [[ ! -f "${fileA}" || ! -f "${fileB}" ]]; then
        echo "FAIL: ${relation} ${nameA} ${nameB}: ${file} missing"
        failures=$((failures + 1))
        failed+=("${relation}:${nameA}:${nameB}")
    elif cmp -s "${fileA}" "${fileB}"; then
        if [[ "${relation}" == same ]]; then
            echo "PASS: ${nameA} and ${nameB} give the same ${file}"
        else
            echo "FAIL: ${nameA} and ${nameB} give the same ${file}, expected a difference"
            failures=$((failures + 1))
            failed+=("${relation}:${nameA}:${nameB}")
        fi
    else
        if [[ "${relation}" == differ ]]; then
            echo "PASS: ${nameA} and ${nameB} give different ${file}"
        else
            echo "FAIL: ${nameA} and ${nameB} give different ${file}, expected the same"
            failures=$((failures + 1))
            failed+=("${relation}:${nameA}:${nameB}")
        fi
    fi
done < "${caseDir}/comparisons"

echo
if (( checks == 0 )); then
    echo "FAIL: no check matched: $(printf '%s ' "${selected[@]}")"
    exit 1
fi
if (( failures > 0 )); then
    echo "Configuration checks FAILED (${failures}/${checks}): $(printf '%s ' "${failed[@]}")"
    exit 1
fi
echo "Configuration checks PASSED (${checks} checks)"
exit 0
