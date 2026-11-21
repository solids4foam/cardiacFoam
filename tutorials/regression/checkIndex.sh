#!/usr/bin/env bash
set -euo pipefail

# Checks that the canonical case table in tutorials/README.md agrees with
# the regression suite. The Regression column of every case row names its
# tier: `standard`, `full` or `none`. A row marked `standard` must be a
# discovered case without regression/full-only, a row marked `full` a
# discovered case with it, and every discovered case must have such a row.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

indexed()
{
    awk -F'|' -v tier="\`$1\`" '
        /^\| \[`/ {
            path = $2; sub(/^ *\[`/, "", path); sub(/`\].*$/, "", path)
            if (index($5, tier)) print path
        }' "${ROOT_DIR}/README.md" | sort
}

listed="$("${ROOT_DIR}/Alltest-regression" --list)"
strip() { sed -E 's/  \((full-only|needs-solids4foam|needs-tools|skips on:[^)]*)\)//g' | sort; }
discoveredStandard="$(grep -v '(full-only)' <<< "${listed}" | strip || true)"
discoveredFull="$(grep '(full-only)' <<< "${listed}" | strip || true)"

status=0
compareTier()
{
    local tier="$1" indexedCases="$2" discoveredCases="$3" missing
    missing="$(comm -13 <(echo "${indexedCases}") <(echo "${discoveredCases}") | grep . || true)"
    if [[ -n "${missing}" ]]; then
        echo "FAIL: in the ${tier} tier but not marked \`${tier}\` in tutorials/README.md:"
        printf '  %s\n' ${missing}
        status=1
    fi
    missing="$(comm -23 <(echo "${indexedCases}") <(echo "${discoveredCases}") | grep . || true)"
    if [[ -n "${missing}" ]]; then
        echo "FAIL: marked \`${tier}\` in tutorials/README.md but not in the ${tier} tier:"
        printf '  %s\n' ${missing}
        status=1
    fi
}
compareTier standard "$(indexed standard)" "${discoveredStandard}"
compareTier full "$(indexed full)" "${discoveredFull}"

if (( status == 0 )); then
    echo "PASS: tutorials/README.md index and Alltest-regression agree" \
         "($(grep -c . <<< "${discoveredStandard}") standard, $(grep -c . <<< "${discoveredFull}") full)"
fi
exit "${status}"
