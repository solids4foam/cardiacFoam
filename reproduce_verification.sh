#!/usr/bin/env bash
# Run or re-check the normalized verification experiments.
# Usage: ./reproduce_verification.sh [--dry-run] [--skip-run] [experiment_id ...]
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="$REPO_ROOT/applications/scripts/paperI_results"
DRIVER="$REPO_ROOT/applications/scripts/driverFoam/bin/driverFoam"

DRY_RUN=0; SKIP_RUN=0; SELECT=()
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --skip-run) SKIP_RUN=1 ;;
        -*) echo "unknown flag: $arg" >&2; exit 64 ;;
        *) SELECT+=("$arg") ;;
    esac
done

if [[ "$DRY_RUN" -eq 0 && "$SKIP_RUN" -eq 0 && -z "${WM_PROJECT_VERSION:-}" ]]; then
    echo "ERROR: OpenFOAM environment not sourced (WM_PROJECT_VERSION unset)." >&2
    exit 3
fi

selected() {
    [[ "${#SELECT[@]}" -eq 0 ]] && return 0
    local candidate
    for candidate in "${SELECT[@]}"; do
        [[ "$candidate" == "$1" ]] && return 0
    done
    return 1
}

rc=0
while IFS=$'\t' read -r experiment_id case_dir runner aggregator result reference; do
    [[ -z "$experiment_id" || "$experiment_id" == "experiment_id" ]] && continue
    selected "$experiment_id" || continue
    case_root="$REPO_ROOT/$case_dir"
    result_path="$case_root/$result"
    echo "== $experiment_id =="

    if [[ "$DRY_RUN" -eq 1 ]]; then
        echo "  run: $case_dir/$runner"
        echo "  result: $case_dir/$result"
        [[ "$aggregator" != "-" ]] && echo "  aggregate: $aggregator"
        [[ "$reference" != "-" ]] && echo "  reference: $case_dir/$reference"
        continue
    fi

    if [[ "$SKIP_RUN" -eq 0 ]]; then
        (cd "$case_root" && bash "$runner") < /dev/null || {
            echo "  RUN FAILED"; rc=1; continue;
        }
        provenance_tmp="$case_root/provenance.json.tmp"
        PROVENANCE_CMD="./reproduce_verification.sh $experiment_id" \
          PROVENANCE_EXIT=0 \
          bash "$TOOLS/capture_provenance.sh" "$case_root" > "$provenance_tmp" || {
            echo "  PROVENANCE CAPTURE FAILED"; rc=1; continue;
          }
        mv "$provenance_tmp" "$case_root/provenance.json"
    fi

    # --skip-run checks the existing normalized result. Re-aggregation belongs
    # to an actual run because archived raw sweep folders are not required to
    # remain in a compact verification checkout.
    if [[ "$SKIP_RUN" -eq 0 && "$aggregator" != "-" ]]; then
        python3 "$TOOLS/aggregate.py" "$aggregator" --repo-root "$REPO_ROOT" || {
            echo "  AGGREGATE FAILED"; rc=1; continue;
        }
    fi
    if [[ ! -s "$result_path" ]]; then
        echo "  FAIL (missing or empty result: $result)"; rc=1; continue
    fi
    if [[ "$reference" == "-" ]]; then
        echo "  PASS (result produced; no numerical reference registered)"
        continue
    fi

    reference_path="$case_root/$reference"
    python3 "$TOOLS/keyset_gate.py" "$result_path" "$reference_path" || {
        echo "  FAIL (canonical row set differs from reference)"; rc=1; continue;
    }
    python3 "$TOOLS/check_against_reference.py" "$result_path" "$reference_path" || {
        echo "  FAIL"; rc=1; continue;
    }
    echo "  PASS"
done < <("$DRIVER" experiment-plan --format tsv)
exit "$rc"
