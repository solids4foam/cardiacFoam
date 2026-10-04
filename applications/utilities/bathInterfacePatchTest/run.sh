#!/usr/bin/env bash

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    # Local default; set WM_PROJECT_DIR beforehand on other installations.
    source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null
fi
set -euo pipefail

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${test_dir}/../../.." && pwd)"
template="${repo_dir}/tutorials/manufacturedSolutions/bathBidomainConormal/system"
scratch="$(mktemp -d "${TMPDIR:-/tmp}/bath-interface-patch.XXXXXX")"
if [[ "${KEEP_CASES:-0}" != 1 ]]; then
    trap 'rm -rf "${scratch}"' EXIT
else
    echo "Keeping generated cases in ${scratch}"
fi

wmake "${test_dir}" >/dev/null

resolutions=("$@")
if (( ${#resolutions[@]} == 0 )); then
    resolutions=(10 20 40)
fi

printf 'N traceRms traceMax fluxRms fluxMax stripResidualRms stripResidualMax\n'
results="${scratch}/results.txt"
for n in "${resolutions[@]}"; do
    case_dir="${scratch}/N${n}"
    mkdir -p "${case_dir}/system" "${case_dir}/constant"
    cp "${template}/controlDict" "${template}/fvSchemes" \
        "${template}/fvSolution" \
        "${template}/blockMeshDict" "${template}/topoSetDict" \
        "${case_dir}/system/"
    python3 - "${case_dir}/system/blockMeshDict" "${n}" <<'PY'
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
n = int(sys.argv[2])
text = path.read_text()
assert text.count("(20 20 1)") == 3
path.write_text(text.replace("(20 20 1)", f"({n} {n} 1)"))
PY
    blockMesh -case "${case_dir}" >"${case_dir}/blockMesh.log" 2>&1
    topoSet -case "${case_dir}" >"${case_dir}/topoSet.log" 2>&1
    "${FOAM_USER_APPBIN}/bathInterfacePatchTest" -case "${case_dir}" \
        >"${case_dir}/patchTest.log" 2>&1
    python3 - "${case_dir}/patchTest.log" "${n}" >>"${results}" <<'PY'
import pathlib
import re
import sys

text = pathlib.Path(sys.argv[1]).read_text()
match = re.search(r"nExposed=.*stripResidualMax=[^\s]+", text)
if not match:
    raise SystemExit(text)
items = dict(part.split("=", 1) for part in match.group().split())
names = (
    "traceRms", "traceMax", "fluxRms", "fluxMax",
    "stripResidualRms", "stripResidualMax",
)
print(sys.argv[2], *(items[name] for name in names))
PY
    tail -n 1 "${results}"
done

control_log="${scratch}/N${resolutions[0]}/normalControl.log"
"${FOAM_USER_APPBIN}/bathInterfacePatchTest" \
    -case "${scratch}/N${resolutions[0]}" -tangentSlope 0 \
    >"${control_log}" 2>&1
python3 - "${control_log}" <<'PY'
import pathlib
import re
import sys

text = pathlib.Path(sys.argv[1]).read_text()
match = re.search(r"nExposed=.*stripResidualMax=[^\s]+", text)
if not match:
    raise SystemExit(text)
items = dict(part.split("=", 1) for part in match.group().split())
names = ("traceMax", "fluxMax", "stripResidualMax")
values = [float(items[name]) for name in names]
if max(values) > 1e-10:
    raise SystemExit("FAIL: normal-only control " + str(dict(zip(names, values))))
print("PASS normal-only control: " + " ".join(f"{name}={value:g}" for name, value in zip(names, values)))
PY

if [[ "${DIAGNOSTIC_ONLY:-0}" != 1 ]]; then
    python3 - "${results}" <<'PY'
import pathlib
import sys

rows = [list(map(float, line.split())) for line in pathlib.Path(sys.argv[1]).read_text().splitlines()]
tol = 1e-10
failed = [
    (int(row[0]), row[2], row[4], row[6])
    for row in rows
    if max(row[2], row[4], row[6]) > tol
]
if failed:
    for n, trace, flux, residual in failed:
        print(
            f"FAIL N={n}: affine interface inconsistency "
            f"(traceMax={trace:g}, fluxMax={flux:g}, "
            f"stripResidualMax={residual:g}; tolerance={tol:g})",
            file=sys.stderr,
        )
    raise SystemExit(1)
print("PASS: affine interface trace, flux, and local residual are consistent")
PY
fi
