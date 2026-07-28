#!/usr/bin/env bash
# Reproduce Paper I convergence tables: provenance -> run -> canonicalize -> diff.
# Usage: ./reproduce_paperI.sh [--dry-run] [--skip-run] [key ...]
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$REPO_ROOT/applications/scripts/paperI_results"
REGISTRY="$PKG/paperI_cases.tsv"

DRY_RUN=0; SKIP_RUN=0; SELECT=()
for a in "$@"; do
  case "$a" in
    --dry-run) DRY_RUN=1 ;;
    --skip-run) SKIP_RUN=1 ;;
    -*) echo "unknown flag: $a" >&2; exit 64 ;;
    *) SELECT+=("$a") ;;
  esac
done

# OF env is only needed to actually run solvers (not for --dry-run or --skip-run).
if [ "$DRY_RUN" -eq 0 ] && [ "$SKIP_RUN" -eq 0 ] && [ -z "${WM_PROJECT_VERSION:-}" ]; then
  echo "ERROR: OpenFOAM environment not sourced (WM_PROJECT_VERSION unset)." >&2
  exit 3
fi

selected() {
  [ "${#SELECT[@]}" -eq 0 ] && return 0
  local s; for s in "${SELECT[@]}"; do [ "$s" = "$1" ] && return 0; done
  return 1
}

rc=0
# skip the header line, then read tab-separated columns (uniform bash: no run_kind)
while IFS=$'\t' read -r key case_dir run_entry agg_key fresh_rel ref_rel; do
  [ -z "$key" ] && continue
  case "$key" in \#*) continue ;; esac
  selected "$key" || continue
  cdir="$REPO_ROOT/$case_dir"
  fresh="$cdir/$fresh_rel"
  ref="$cdir/$ref_rel"
  echo "== $key =="
  if [ "$DRY_RUN" -eq 1 ]; then
    echo "  run: bash $run_entry"
    echo "  canonicalize: aggregate.py $agg_key -> $fresh_rel"
    echo "  diff vs: $ref_rel"
    continue
  fi
  PROVENANCE_CMD="reproduce_paperI.sh $key" \
    bash "$PKG/capture_provenance.sh" "$cdir" < /dev/null > "$cdir/provenance.json" || \
    echo "  (provenance capture warning)"
  # uniform bash: run_entry is a script under the case dir; "-" means no run step
  # < /dev/null is required: without it, run_entry inherits the same stdin fd
  # as this while-loop's `read` (fed by the process substitution below), and
  # if anything the run entry invokes so much as peeks at stdin, registry
  # lines vanish silently -- cases get skipped with no error at all.
  if [ "$SKIP_RUN" -eq 0 ] && [ "$run_entry" != "-" ]; then
    ( cd "$cdir" && bash "$run_entry" ) < /dev/null || { echo "  RUN FAILED"; rc=1; continue; }
  fi
  if ! agg_err=$(python3 "$PKG/aggregate.py" "$agg_key" 2>&1 >/dev/null); then
    if [ ! -f "$ref" ]; then echo "  SKIP (no data + no reference yet)"; continue; fi
    echo "  AGGREGATE FAILED"; echo "$agg_err" | sed 's/^/    /'; rc=1; continue
  fi
  python3 "$PKG/check_against_reference.py" "$fresh" "$ref"
  case $? in
    0) echo "  PASS" ;;
    2) echo "  SKIP (no committed reference yet)" ;;
    *) echo "  FAIL"; rc=1 ;;
  esac
done < <(tail -n +2 "$REGISTRY")
exit $rc
