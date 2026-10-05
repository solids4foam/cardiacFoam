#!/usr/bin/env bash
set -euo pipefail

# Checks that the canonical case table in tutorials/README.md agrees with
# the regression suite: every row whose Regression column names
# Alltest-regression must be a discovered case, and every discovered case
# must have such a row.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

indexed="$(awk -F'|' '
    /^\| `/ {
        path = $2; gsub(/[` ]/, "", path)
        if ($6 ~ /Alltest-regression/) print path
    }' "${ROOT_DIR}/README.md" | sort)"

discovered="$("${ROOT_DIR}/Alltest-regression" --list | sed -E 's/  \((needs-solids4foam|skips on:[^)]*)\)//g' | sort)"

missingFromIndex="$(comm -13 <(echo "${indexed}") <(echo "${discovered}"))"
missingFromSuite="$(comm -23 <(echo "${indexed}") <(echo "${discovered}"))"

status=0
if [[ -n "${missingFromIndex}" ]]; then
    echo "FAIL: discovered by Alltest-regression but not marked in tutorials/README.md:"
    printf '  %s\n' ${missingFromIndex}
    status=1
fi
if [[ -n "${missingFromSuite}" ]]; then
    echo "FAIL: marked Alltest-regression in tutorials/README.md but not discovered:"
    printf '  %s\n' ${missingFromSuite}
    status=1
fi
if (( status == 0 )); then
    echo "PASS: tutorials/README.md index and Alltest-regression agree ($(echo "${discovered}" | grep -c .) cases)"
fi
exit "${status}"
