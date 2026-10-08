# shellcheck shell=bash
# ============================================================
# cardiacFoam tutorial regression library
# ============================================================
#
# Sourced by every tutorial's regression/regressionTest.sh. A case script
# declares its reference file and how to run the case; this file owns the
# option parsing, the build-mode gate, the run step, the reference grammar,
# the comparator and the JSON report.
#
# Case script skeleton:
#
#     #!/usr/bin/env bash
#     set -euo pipefail
#     . "$(dirname "${BASH_SOURCE[0]}")/../../../regression/lib.sh"
#     regression_init "Case title" regression/case.reference "$@"
#     regression_set system/controlDict endTime 0.05    # regression-only
#     regression_run parallel
#     regression_compare
#     regression_finish
#
# Options accepted by regression_init (from the script's "$@"):
#     --check-only     compare existing outputs; regression_run does nothing
#     --report PATH    write the comparison evidence as JSON to PATH
#     --help
#
# Environment:
#     CARDIAC_REGRESSION_BUILD_MODE   with-solids4foam | lightweight
#         Set by Alltest-regression. regression_require_solids4foam exits
#         with REGRESSION_SKIP_CODE (77) under lightweight.
#     CARDIAC_REGRESSION_SCOPE        standard | full   (default standard)
#         Read by cases that offer a wider check than the suite gate, see
#         regression_scope.
#
# Reference grammar. One check per line, whitespace separated, '#' starts
# a comment. Paths are relative to the case root; a path that does not
# exist is also looked up under processor*/ (parallel runs that were not
# reconstructed). The last two fields of every row are the expected value
# and the absolute tolerance; a check passes when
# |actual - expected| <= tolerance.
#
#     probe   <file> <time> <column> <expected> <tolerance>
#         Whitespace-separated time series (OpenFOAM probes, surfaceFieldValue,
#         ECG traces, cardiacFoam .txt traces). Row = the one whose first
#         field is nearest to <time>, within REGRESSION_TIME_WINDOW.
#         <column> is a 1-based index after parentheses are stripped (so a
#         vector probe spans three columns) or a header name; a name is
#         matched against the first header line, with or without a leading
#         '#' and with or without a 'numeric_' prefix.
#     csv     <file> <key> <column> <expected> <tolerance>
#         Comma-separated table with a header line. Row = nearest first-field
#         value to <key> within REGRESSION_TIME_WINDOW; <column> by name.
#     final   <file> <column> <expected> <tolerance>
#         Last data row of a whitespace-separated series; <column> by name
#         or index.
#     summary <key> <expected> <tolerance>
#         Scalar from the manufactured-solution summary file
#         (REGRESSION_SUMMARY_FILE). Keys: cells, cellsPerDirection, finalTime.
#     error   <field> <L1|L2|Linf> <expected> <tolerance>
#         Error norm of <field> from the summary file's table.
#     <other> <fields...>
#         Passed to the case's regression_case_check function, if defined.
#
# A case may set before regression_compare:
#     REGRESSION_TIME_WINDOW    nearest-row window for probe/csv (default 2.5e-3)
#     REGRESSION_SUMMARY_FILE   file read by the summary and error kinds

IFS=$' \t\n'

REGRESSION_SKIP_CODE=77
REGRESSION_EXIT_HOOKS=()

REGRESSION_TITLE=""
REGRESSION_REF_FILE=""
REGRESSION_CHECK_ONLY=0
REGRESSION_REPORT_PATH=""
REGRESSION_REPORT_ROWS=""
REGRESSION_TIME_WINDOW="2.5e-3"
REGRESSION_SUMMARY_FILE=""
REGRESSION_SCOPE=standard
REGRESSION_ALLRUN_LOG="log.Allrun"
REGRESSION_CHECKS=0
REGRESSION_FAILURES=0
REGRESSION_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[1]:-$0}")" && pwd)"
REGRESSION_CASE_DIR="$(cd "${REGRESSION_SCRIPT_DIR}/.." && pwd)"
REGRESSION_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

regression_usage()
{
    cat <<EOF
Usage: regressionTest.sh [--check-only] [--report PATH]

Without options, clean the case, run Allrun, and compare its outputs to the
reference data ${REGRESSION_REF_FILE}. --check-only performs the same
comparison on existing outputs without cleaning or running the case.
--report writes the comparison evidence as JSON.
EOF
}

# regression_init <title> <reference file> [script arguments...]
regression_init()
{
    REGRESSION_TITLE="$1"
    REGRESSION_REF_FILE="$2"
    shift 2

    while (( $# > 0 )); do
        case "$1" in
            --check-only)
                REGRESSION_CHECK_ONLY=1
                ;;
            --report)
                if (( $# < 2 )); then
                    echo "FAIL: --report requires a path" >&2
                    exit 2
                fi
                REGRESSION_REPORT_PATH="$2"
                shift
                ;;
            --help|-h)
                regression_usage
                exit 0
                ;;
            *)
                echo "FAIL: unknown option: $1" >&2
                regression_usage >&2
                exit 2
                ;;
        esac
        shift
    done

    cd "${REGRESSION_CASE_DIR}" || exit 1

    trap regression_run_exit_hooks EXIT
    if [[ -n "${REGRESSION_REPORT_PATH}" ]]; then
        REGRESSION_REPORT_ROWS="$(mktemp)"
        regression_add_exit_hook 'rm -f "${REGRESSION_REPORT_ROWS}"'
    fi

    echo "============================================================"
    echo "${REGRESSION_TITLE}"
    echo "============================================================"
    echo
}

# regression_add_exit_hook <command>: run <command> when the script exits,
# whatever the exit path. Hooks run in reverse order of registration.
regression_add_exit_hook()
{
    REGRESSION_EXIT_HOOKS+=("$1")
}

regression_run_exit_hooks()
{
    local i
    for (( i = ${#REGRESSION_EXIT_HOOKS[@]} - 1; i >= 0; i-- )); do
        eval "${REGRESSION_EXIT_HOOKS[i]}"
    done
}

# regression_fail <message>: one failed check that is not a reference row.
regression_fail()
{
    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1))
    REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
    echo "FAIL: $1"
    regression_report_row check "" "" "" "" "" "" "" failed "$(regression_json_string "$1")"
}

# regression_scope: sets REGRESSION_SCOPE to standard or full from
# CARDIAC_REGRESSION_SCOPE; exits on any other value. --check-only always
# selects standard, since only one set of outputs exists to compare.
regression_scope()
{
    REGRESSION_SCOPE="${CARDIAC_REGRESSION_SCOPE:-standard}"
    case "${REGRESSION_SCOPE}" in
        standard|full) ;;
        *)
            echo "FAIL: CARDIAC_REGRESSION_SCOPE must be 'standard' or 'full', got '${REGRESSION_SCOPE}'." >&2
            exit 2
            ;;
    esac
    if (( REGRESSION_CHECK_ONLY )); then
        REGRESSION_SCOPE=standard
    fi
    echo "Scope: ${REGRESSION_SCOPE}"
}

# regression_require_run_mode: for cases whose variants cannot be compared
# against one set of existing outputs.
regression_require_run_mode()
{
    if (( REGRESSION_CHECK_ONLY )); then
        echo "FAIL: --check-only is not supported by this case: each variant needs its own run." >&2
        exit 2
    fi
}

# ------------------------------------------------------------
# Build-mode gate
# ------------------------------------------------------------

regression_solids4foam_tree()
{
    local dir="${SOLIDS4FOAM_INST_DIR:-}"
    local repoRoot
    local candidate
    local marker="src/solids4FoamModels/lnInclude/solidModel.H"

    if [[ -n "${dir}" && -f "${dir}/src/solids4FoamModels/solidModels/solidModel/solidModel.H" ]]; then
        return 0
    fi

    repoRoot="$(cd "${REGRESSION_LIB_DIR}/../.." && pwd)"
    for candidate in \
        "${repoRoot}/modules/solids4foam" \
        "${HOME}/solids4foam" \
        "${WM_PROJECT_USER_DIR:-}/solids4foam"
    do
        if [[ -n "${candidate}" && -f "${candidate}/${marker}" ]]; then
            export SOLIDS4FOAM_INST_DIR="${candidate}"
            return 0
        fi
    done

    return 1
}

regression_electromechanics_lib_built()
{
    local lib
    [[ -n "${FOAM_USER_LIBBIN:-}" ]] || return 1
    for lib in "${FOAM_USER_LIBBIN}"/libelectroMechanicalModels.*; do
        [[ -e "${lib}" ]] && return 0
    done
    return 1
}

# regression_require_solids4foam: skip (77) under lightweight mode or without
# a solids4foam tree; fail when the tree exists but cardiacFoam was built
# without libelectroMechanicalModels.
regression_require_solids4foam()
{
    if [[ "${CARDIAC_REGRESSION_BUILD_MODE:-}" == "lightweight" ]]; then
        echo "SKIP: this regression requires a full solids4foam build, but lightweight mode was specified."
        exit "${REGRESSION_SKIP_CODE}"
    fi

    if ! regression_solids4foam_tree; then
        echo "SKIP: this regression requires a full solids4foam build."
        echo "      SOLIDS4FOAM_INST_DIR does not point to a compiled solids4foam tree."
        exit "${REGRESSION_SKIP_CODE}"
    fi

    if ! regression_electromechanics_lib_built; then
        echo "FAIL: solids4foam is available, but libelectroMechanicalModels is not compiled."
        echo "      Rebuild cardiacFoam in full mode before running this regression."
        exit 1
    fi
}

# ------------------------------------------------------------
# Regression configuration
# ------------------------------------------------------------

# Regression-only settings. A regression checks that the case still runs and
# that its outputs have not changed, so it may run a smaller configuration
# than the tutorial: a shorter endTime, a coarser mesh, fewer electrodes.
# These helpers edit a case file for the script's invocation only: the file
# is restored when the script exits, so the tutorial keeps its own settings.
# They do nothing under --check-only, and an edit that changes nothing fails
# the script, so a setting cannot silently go missing.

regression_backup_file()
{
    local file="$1" backup="$1.regressionTest.bak"
    if [[ ! -e "${backup}" ]]; then
        cp -p "${file}" "${backup}"
        regression_add_exit_hook "mv -f $(printf '%q' "${backup}") $(printf '%q' "${file}")"
    fi
}

# regression_edit <file> <sed -E expression>...
regression_edit()
{
    local file="$1" expr
    shift

    (( REGRESSION_CHECK_ONLY )) && return 0
    regression_backup_file "${file}"

    for expr in "$@"; do
        sed -E -e "${expr}" "${file}" > "${file}.regressionTest.tmp"
        if cmp -s "${file}" "${file}.regressionTest.tmp"; then
            rm -f "${file}.regressionTest.tmp"
            regression_fail "regression setting '${expr}' changes nothing in ${file}"
            regression_finish
        fi
        mv -f "${file}.regressionTest.tmp" "${file}"
        echo "Regression setting: ${file}: ${expr}"
    done
}

# regression_set <file> <keyword> <value>: sets the one entry of the file
# whose line starts with <keyword> (at any indentation) to <value>.
regression_set()
{
    local file="$1" keyword="$2" value="$3" count

    (( REGRESSION_CHECK_ONLY )) && return 0

    count="$(grep -cE "^[[:space:]]*${keyword}[[:space:]]+[^;]*;" "${file}" || true)"
    if (( count != 1 )); then
        regression_fail "regression setting ${keyword} ${value}: ${file} has ${count} '${keyword}' entries, expected one"
        regression_finish
    fi
    regression_edit "${file}" "s/^([[:space:]]*${keyword}[[:space:]]+)[^;]*;/\\1${value};/"
}

# ------------------------------------------------------------
# Running the case
# ------------------------------------------------------------

regression_log_tail()
{
    local label="$1" logFile="$2" maxLines="${3:-80}"

    if [[ -s "${logFile}" ]]; then
        echo "----- last ${maxLines} lines of ${label} (${logFile}) -----"
        tail -n "${maxLines}" "${logFile}"
        echo "----- end of ${label} -----"
    else
        echo "(no log file at ${logFile})"
    fi
}

regression_dump_logs()
{
    local logFile
    regression_log_tail "Allrun" "${REGRESSION_ALLRUN_LOG}"
    for logFile in log.*; do
        [[ "${logFile}" == "${REGRESSION_ALLRUN_LOG}" ]] && continue
        regression_log_tail "${logFile#log.}" "${logFile}" 40
    done
}

# regression_run [Allrun arguments...]
# Cleans the case and runs Allrun. Returns 1 (after printing the logs) when
# Allrun exits non-zero, log.cardiacFoam does not end with "End", or any
# log.* reports a fatal error.
# Does nothing under --check-only.
regression_run()
{
    if (( REGRESSION_CHECK_ONLY )); then
        echo "Checking existing outputs (no cleanup or run)"
        return 0
    fi

    ./Allclean > /dev/null 2>&1 || true

    if ! ./Allrun "$@" > "${REGRESSION_ALLRUN_LOG}" 2>&1; then
        echo "FAIL: Allrun $* exited non-zero. Surfacing logs:"
        regression_dump_logs
        return 1
    fi

    if ! regression_check_solver_logs; then
        echo "FAIL: Allrun $* finished but a solver log is incomplete or reports a fatal error. Surfacing logs:"
        regression_dump_logs
        return 1
    fi

    return 0
}

# regression_check_solver_logs [-e <regex>] [log ...]
# Every named log (default log.cardiacFoam, when present) must contain a line
# matching <regex> (default the OpenFOAM '^End' line), and no log.* file in
# the case may report "FOAM FATAL" or "FOAM aborting". Returns 1 otherwise.
regression_check_solver_logs()
{
    local endPattern='^End' logFile failures=0
    if [[ "${1:-}" == "-e" ]]; then
        endPattern="$2"
        shift 2
    fi
    local -a logs=("$@")
    if (( ${#logs[@]} == 0 )) && [[ -f log.cardiacFoam ]]; then
        logs=(log.cardiacFoam)
    fi
    for logFile in "${logs[@]+"${logs[@]}"}"; do
        if [[ ! -f "${logFile}" ]]; then
            echo "FAIL: expected log file not found: ${logFile}"
            failures=$((failures + 1))
        elif ! grep -Eq "${endPattern}" "${logFile}"; then
            echo "FAIL: ${logFile} has no line matching '${endPattern}' (run did not complete)"
            failures=$((failures + 1))
        fi
    done
    for logFile in log.*; do
        [[ -f "${logFile}" ]] || continue
        if grep -Eq 'FOAM FATAL|FOAM aborting' "${logFile}"; then
            echo "FAIL: ${logFile} reports a fatal error"
            failures=$((failures + 1))
        fi
    done
    (( failures == 0 ))
}

# regression_run_or_fail [Allrun arguments...]: regression_run, exiting 1 on failure.
regression_run_or_fail()
{
    regression_run "$@" || { regression_write_report; exit 1; }
}

# ------------------------------------------------------------
# Output discovery
# ------------------------------------------------------------

# regression_resolve_file <path>: prints <path>, or processor*/<path> when the
# serial location does not exist. Returns 1 when neither exists.
regression_resolve_file()
{
    local path="$1" candidate

    if [[ -s "${path}" ]]; then
        echo "${path}"
        return 0
    fi

    for candidate in processor*/"${path}"; do
        if [[ -s "${candidate}" ]]; then
            echo "${candidate}"
            return 0
        fi
    done

    return 1
}

# regression_find_output <text> <glob>...: first non-empty file matching one
# of the globs (serial location, then processor*/) that contains <text>.
regression_find_output()
{
    local text="$1" pattern candidate
    shift

    for pattern in "$@"; do
        for candidate in ${pattern} processor*/${pattern}; do
            if [[ -s "${candidate}" ]] && grep -q -- "${text}" "${candidate}"; then
                echo "${candidate}"
                return 0
            fi
        done
    done

    return 1
}

# regression_require_output <label> <file> [text]: fails the script when the
# file is missing (or does not contain <text>), otherwise prints a PASS line.
regression_require_output()
{
    local label="$1" file="$2" text="${3:-}"

    if [[ -z "${file}" || ! -s "${file}" ]]; then
        echo "FAIL: ${label} not found${file:+ at ${file}}"
        [[ -f log.cardiacFoam ]] && regression_log_tail cardiacFoam log.cardiacFoam 40
        regression_write_report
        exit 1
    fi
    if [[ -n "${text}" ]] && ! grep -q -- "${text}" "${file}"; then
        echo "FAIL: ${label} at ${file} does not contain '${text}'"
        regression_write_report
        exit 1
    fi
    echo "PASS: ${label} detected in ${file}"
}

# ------------------------------------------------------------
# Value extraction
# ------------------------------------------------------------

# Column lookup shared by the awk programs below. header_column() resolves a
# name against the first header line: "# time E1 E2" maps E1 to column 2
# (the lone '#' is dropped), "#time E1" and "time E1" map E1 to column 2.
regression_awk_common='
function strip_hash(s) { sub(/^#/, "", s); return s }
function header_column(line, name,    n, i, parts, offset, token) {
    n = split(line, parts, FS)
    offset = 0
    if (parts[1] == "#") offset = 1
    for (i = 1; i <= n; i++) {
        token = strip_hash(parts[i])
        sub(/^numeric_/, "", token)
        if (token == name) return i - offset
    }
    return 0
}
function is_number(s) { return s ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/ }
'

# regression_probe_value <file> <time> <column|name>
regression_probe_value()
{
    local file="$1" time="$2" selector="$3"

    awk -v target="${time}" -v selector="${selector}" -v window="${REGRESSION_TIME_WINDOW}" \
        "${regression_awk_common}"'
        BEGIN { bestDiff = 1e99; found = 0; col = 0; headerDone = 0 }
        {
            line = $0
            if (!headerDone && !is_number($1)) {
                if (col == 0 && !is_number(selector)) {
                    c = header_column(line, selector)
                    if (c > 0) col = c
                }
                next
            }
            headerDone = 1
            if (col == 0) {
                if (is_number(selector)) col = selector + 0
                else exit 1
            }
            gsub(/[()]/, " ", line)
            n = split(line, f, " ")
            if (n < col) next
            d = f[1] - target
            if (d < 0) d = -d
            if (d < bestDiff) { bestDiff = d; actual = f[col]; found = 1 }
        }
        END {
            if (found && bestDiff <= window) { print actual; exit 0 }
            exit 1
        }
    ' "${file}"
}

# regression_csv_value <file> <key> <column name>
regression_csv_value()
{
    local file="$1" key="$2" name="$3"

    awk -F, -v target="${key}" -v name="${name}" -v window="${REGRESSION_TIME_WINDOW}" \
        "${regression_awk_common}"'
        BEGIN { bestDiff = 1e99; found = 0; col = 0 }
        NR == 1 { col = header_column($0, name); next }
        col > 0 && NF >= col {
            d = $1 - target
            if (d < 0) d = -d
            if (d < bestDiff) { bestDiff = d; actual = $col; found = 1 }
        }
        END {
            if (found && bestDiff <= window) { print actual; exit 0 }
            exit 1
        }
    ' "${file}"
}

# regression_final_value <file> <column|name>
regression_final_value()
{
    local file="$1" selector="$2"

    awk -v selector="${selector}" "${regression_awk_common}"'
        BEGIN { col = 0; found = 0 }
        {
            if (!is_number($1)) {
                if (col == 0 && !is_number(selector)) {
                    c = header_column($0, selector)
                    if (c > 0) col = c
                }
                next
            }
            if (col == 0) {
                if (is_number(selector)) col = selector + 0
                else exit 1
            }
            if (NF >= col) { value = $col; found = 1 }
        }
        END {
            if (found) { print value; exit 0 }
            exit 1
        }
    ' "${file}"
}

# regression_summary_value <key>: reads REGRESSION_SUMMARY_FILE. The three
# verifier summary formats are recognised:
#   "Number of cells = N" / "Final simulation time = t"      (field verifiers)
#   "# cellsPerDirection N" / "# time t"                     (bath verifier)
#   "... error summary (t = t):"                             (electromechanics)
regression_summary_value()
{
    local key="$1" file="${REGRESSION_SUMMARY_FILE}"
    [[ -s "${file}" ]] || return 1

    case "${key}" in
        cells)
            awk -F= '/Number of cells/ { gsub(/[[:space:]]/, "", $2); print $2; exit }' "${file}"
            ;;
        cellsPerDirection)
            awk '/^# cellsPerDirection/ { print $3; exit }' "${file}"
            ;;
        finalTime)
            awk -F= '/Final simulation time/ { gsub(/[[:space:]]/, "", $2); print $2; exit }' "${file}" \
                | grep . \
                || awk '/^# time/ { print $3; exit }' "${file}" | grep . \
                || sed -nE 's/.*error summary \(t = ([^)]+)\):.*/\1/p' "${file}" | head -n 1
            ;;
        *)
            return 1
            ;;
    esac
}

# regression_error_metric <field> <L1|L2|Linf>: reads REGRESSION_SUMMARY_FILE.
regression_error_metric()
{
    local field="$1" metric="$2" column file="${REGRESSION_SUMMARY_FILE}"
    [[ -s "${file}" ]] || return 1

    case "${metric}" in
        L1)   column=2 ;;
        L2)   column=3 ;;
        Linf) column=4 ;;
        *) return 1 ;;
    esac

    awk -v field="${field}" -v column="${column}" '$1 == field { print $column; exit }' "${file}"
}

# ------------------------------------------------------------
# Comparison and report
# ------------------------------------------------------------

regression_json_string()
{
    local s="$1"
    s=${s//\\/\\\\}
    s=${s//\"/\\\"}
    s=${s//$'\n'/\\n}
    printf '%s' "$s"
}

regression_json_number()
{
    if [[ -n "$1" ]]; then printf '%s' "$1"; else printf 'null'; fi
}

# regression_report_row <kind> <file> <time> <column> <expected> <tolerance> <actual> <difference> <status> <reason>
regression_report_row()
{
    [[ -n "${REGRESSION_REPORT_ROWS}" ]] || return 0
    printf '{"kind":"%s","file":"%s","time":%s,"column":"%s","expected":%s,"tolerance":%s,"actual":%s,"difference":%s,"status":"%s","reason":"%s"}\n' \
        "$1" "$(regression_json_string "$2")" "$(regression_json_number "$3")" \
        "$(regression_json_string "$4")" "$(regression_json_number "$5")" \
        "$(regression_json_number "$6")" "$(regression_json_number "$7")" \
        "$(regression_json_number "$8")" "$9" "${10}" >> "${REGRESSION_REPORT_ROWS}"
}

regression_write_report()
{
    [[ -n "${REGRESSION_REPORT_PATH}" ]] || return 0

    mkdir -p "$(dirname "${REGRESSION_REPORT_PATH}")"
    {
        printf '{\n'
        printf '  "schema_version": 1,\n'
        printf '  "mode": "%s",\n' "$( (( REGRESSION_CHECK_ONLY )) && printf 'check-only' || printf 'run-and-check')"
        printf '  "reference_file": "%s",\n' "$(regression_json_string "${REGRESSION_REF_FILE}")"
        printf '  "status": "%s",\n' "$( (( REGRESSION_FAILURES == 0 )) && printf 'passed' || printf 'failed')"
        printf '  "checks": %s,\n' "${REGRESSION_CHECKS}"
        printf '  "failures": %s,\n' "${REGRESSION_FAILURES}"
        printf '  "results": ['
        if [[ -n "${REGRESSION_REPORT_ROWS}" && -s "${REGRESSION_REPORT_ROWS}" ]]; then
            paste -sd, "${REGRESSION_REPORT_ROWS}"
        fi
        printf ']\n}\n'
    } > "${REGRESSION_REPORT_PATH}"
}

regression_abs_diff()
{
    awk -v a="$1" -v e="$2" 'BEGIN { d = a - e; if (d < 0) d = -d; print d }'
}

# regression_check <label> <actual> <expected> <tolerance> [kind file time column]
# One comparison; counts it, prints PASS/FAIL and records the report row.
regression_check()
{
    local label="$1" actual="$2" expected="$3" tolerance="$4"
    local kind="${5:-check}" file="${6:-}" time="${7:-}" column="${8:-}"
    local diffAbs status reason

    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1))

    if [[ -z "${actual}" ]]; then
        echo "FAIL: ${label}: value not found"
        REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
        regression_report_row "${kind}" "${file}" "${time}" "${column}" "${expected}" "${tolerance}" "" "" failed value-not-found
        return 1
    fi

    diffAbs="$(regression_abs_diff "${actual}" "${expected}")"

    if awk -v d="${diffAbs}" -v t="${tolerance}" 'BEGIN { exit !(d <= t) }'; then
        status=PASS; reason=within-tolerance
    else
        status=FAIL; reason=outside-tolerance
        REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
    fi

    printf "%s: %s actual=%.9g expected=%.9g difference=%.3g tolerance=%.3g\n" \
        "${status}" "${label}" "${actual}" "${expected}" "${diffAbs}" "${tolerance}"
    regression_report_row "${kind}" "${file}" "${time}" "${column}" "${expected}" "${tolerance}" \
        "${actual}" "${diffAbs}" "$(tr 'A-Z' 'a-z' <<< "${status}")ed" "${reason}"

    [[ "${status}" == PASS ]]
}

regression_missing()
{
    local label="$1" kind="$2" file="$3" expected="$4" tolerance="$5"
    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1))
    REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
    echo "FAIL: ${label}: missing output file ${file}"
    regression_report_row "${kind}" "${file}" "" "" "${expected}" "${tolerance}" "" "" failed missing-output
}

# regression_compare [reference file]: one check per reference row.
regression_compare()
{
    local refFile="${1:-${REGRESSION_REF_FILE}}"
    local kind rest
    local -a f
    local file resolved actual expected tolerance n

    if [[ ! -f "${refFile}" ]]; then
        echo "FAIL: reference file not found: ${refFile}"
        REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1))
        REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
        return 1
    fi

    while IFS= read -r line || [[ -n "${line}" ]]; do
        line="${line%%#*}"
        read -r -a f <<< "${line}"
        n=${#f[@]}
        (( n == 0 )) && continue
        kind="${f[0]}"

        case "${kind}" in
            probe|csv)
                if (( n != 6 )); then
                    echo "FAIL: ${refFile}: '${kind}' rows need 5 fields: ${line}"
                    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1)); REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
                    continue
                fi
                file="${f[1]}"; expected="${f[4]}"; tolerance="${f[5]}"
                if ! resolved="$(regression_resolve_file "${file}")"; then
                    regression_missing "${kind} ${file} ${f[2]} ${f[3]}" "${kind}" "${file}" "${expected}" "${tolerance}"
                    continue
                fi
                if [[ "${kind}" == probe ]]; then
                    actual="$(regression_probe_value "${resolved}" "${f[2]}" "${f[3]}")" || actual=""
                else
                    actual="$(regression_csv_value "${resolved}" "${f[2]}" "${f[3]}")" || actual=""
                fi
                regression_check "${kind} ${file} t=${f[2]} ${f[3]}" "${actual}" "${expected}" "${tolerance}" \
                    "${kind}" "${file}" "${f[2]}" "${f[3]}" || true
                ;;
            final)
                if (( n != 5 )); then
                    echo "FAIL: ${refFile}: 'final' rows need 4 fields: ${line}"
                    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1)); REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
                    continue
                fi
                file="${f[1]}"; expected="${f[3]}"; tolerance="${f[4]}"
                if ! resolved="$(regression_resolve_file "${file}")"; then
                    regression_missing "final ${file} ${f[2]}" final "${file}" "${expected}" "${tolerance}"
                    continue
                fi
                actual="$(regression_final_value "${resolved}" "${f[2]}")" || actual=""
                regression_check "final ${file} ${f[2]}" "${actual}" "${expected}" "${tolerance}" \
                    final "${file}" "" "${f[2]}" || true
                ;;
            summary)
                if (( n != 4 )); then
                    echo "FAIL: ${refFile}: 'summary' rows need 3 fields: ${line}"
                    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1)); REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
                    continue
                fi
                actual="$(regression_summary_value "${f[1]}")" || actual=""
                regression_check "summary ${f[1]}" "${actual}" "${f[2]}" "${f[3]}" \
                    summary "${REGRESSION_SUMMARY_FILE}" "" "${f[1]}" || true
                ;;
            error)
                if (( n != 5 )); then
                    echo "FAIL: ${refFile}: 'error' rows need 4 fields: ${line}"
                    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1)); REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
                    continue
                fi
                actual="$(regression_error_metric "${f[1]}" "${f[2]}")" || actual=""
                regression_check "error ${f[1]} ${f[2]}" "${actual}" "${f[3]}" "${f[4]}" \
                    error "${REGRESSION_SUMMARY_FILE}" "" "${f[1]} ${f[2]}" || true
                ;;
            *)
                if declare -F regression_case_check > /dev/null; then
                    regression_case_check "${f[@]}" || true
                else
                    echo "FAIL: ${refFile}: unknown reference kind '${kind}'"
                    REGRESSION_CHECKS=$((REGRESSION_CHECKS + 1)); REGRESSION_FAILURES=$((REGRESSION_FAILURES + 1))
                fi
                ;;
        esac
    done < "${refFile}"

    (( REGRESSION_FAILURES == 0 ))
}

# regression_finish: summary line, report, exit status.
regression_finish()
{
    echo
    echo "============================================================"
    if (( REGRESSION_FAILURES == 0 )); then
        echo "Regression test PASSED (${REGRESSION_CHECKS} checks)"
        echo "============================================================"
        regression_write_report
        exit 0
    fi
    echo "Regression test FAILED (${REGRESSION_FAILURES}/${REGRESSION_CHECKS} checks)"
    echo "============================================================"
    regression_write_report
    exit 1
}
