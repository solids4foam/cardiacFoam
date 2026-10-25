#!/usr/bin/env bash
# Usage: run.sh [hex|tet] [N ...]
# Meshes the monodomain sphericalShell case at each N, runs wallGradientProbe
# serial and on 4 scotch ranks, and checks the exact properties.

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    # Local default; set WM_PROJECT_DIR beforehand on other installations.
    source /Volumes/OpenFOAM-v2412/etc/bashrc >/dev/null
fi
. "${WM_PROJECT_DIR}/bin/tools/RunFunctions"
set -euo pipefail

mesh="${1:-hex}"
shift || true
resolutions=("$@")
if (( ${#resolutions[@]} == 0 )); then
    if [[ "${mesh}" == hex ]]; then resolutions=(8 16 32); else resolutions=(6 10 14); fi
fi

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${test_dir}/../../.." && pwd)"
template="${repo_dir}/tutorials/manufacturedSolutions/monodomainPseudoECG/sphericalShell"
scratch="$(mktemp -d "${TMPDIR:-/tmp}/wall-gradient-probe.XXXXXX")"
if [[ "${KEEP_CASES:-0}" != 1 ]]; then
    trap 'rm -rf "${scratch}"' EXIT
else
    echo "Keeping generated cases in ${scratch}"
fi

wmake "${test_dir}" >/dev/null

results="${scratch}/results.txt"
for n in "${resolutions[@]}"; do
    case_dir="${scratch}/${mesh}_N${n}"
    mkdir -p "${case_dir}/system" "${case_dir}/constant" "${case_dir}/local"
    cp "${template}/system/"{blockMeshDict,controlDict,fvSolution,decomposeParDict} \
        "${case_dir}/system/"
    cat > "${case_dir}/system/fvSchemes" <<'EOF'
FoamFile { version 2; format ascii; class dictionary; object fvSchemes; }
ddtSchemes { default none; }
gradSchemes
{
    default                 none;
    grad(psi)               sealedWallLeastSquares;
    grad(psi.component(0))  sealedWallLeastSquares;
    grad(stock)             sealedWallLeastSquares;
    grad(reference)         leastSquares;
}
divSchemes { default none; }
laplacianSchemes { default none; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
EOF
    echo "N ${n};" > "${case_dir}/local/mesh"
    if [[ "${mesh}" == hex ]]; then
        blockMesh -case "${case_dir}" >"${case_dir}/log.blockMesh" 2>&1
    else
        lc=$(awk -v N="${n}" 'BEGIN { printf "%.17g", 1/N }')
        gmsh -3 "${template}/shell.geo.template" -setnumber lc "${lc}" \
            -setnumber seed 1 -o "${case_dir}/local/shell.msh" \
            >"${case_dir}/log.gmsh" 2>&1
        gmshToFoam -case "${case_dir}" "${case_dir}/local/shell.msh" \
            >"${case_dir}/log.gmshToFoam" 2>&1
    fi
    decomposePar -case "${case_dir}" >"${case_dir}/log.decomposePar" 2>&1

    for f in affine shell; do
        wallGradientProbe -case "${case_dir}" -field "${f}" \
            >"${case_dir}/log.${f}.serial" 2>&1
        mpirun -np 4 wallGradientProbe -case "${case_dir}" -field "${f}" \
            -parallel >"${case_dir}/log.${f}.parallel" 2>&1
        for run in serial parallel; do
            echo "${n} ${run} $(grep '^field ' "${case_dir}/log.${f}.${run}")" \
                >>"${results}"
        done
    done
done

python3 - "${results}" <<'PY'
import math
import sys

rows = []
for line in open(sys.argv[1]):
    words = line.split()
    item = dict(zip(words[2::2], words[3::2]))
    item["N"], item["run"] = int(words[0]), words[1]
    rows.append(item)

failed = []
for r in rows:
    tag = f"N={r['N']} {r['field']} {r['run']}"
    if float(r["stockDiff"]) != 0:
        failed.append(f"{tag}: stockDiff {r['stockDiff']} (must be 0)")
    if float(r["nameDiff"]) != 0:
        failed.append(f"{tag}: nameDiff {r['nameDiff']} (must be 0)")
    if r["field"] == "affine":
        for key in ("wallMax", "interiorMax"):
            if float(r[key]) > 1e-10:
                failed.append(f"{tag}: affine {key} {r[key]} (> 1e-10)")

shell = {(r["N"], r["run"]): r for r in rows if r["field"] == "shell"}
for (n, run), r in sorted(shell.items()):
    if run != "parallel":
        continue
    s = shell[(n, "serial")]
    for key in ("wallRMS", "wallMax", "interiorRMS", "interiorMax"):
        a, b = float(s[key]), float(r[key])
        if abs(a - b) > 1e-10*max(abs(a), 1e-300):
            failed.append(f"N={n} shell: serial {key} {a} != parallel {b}")

print("N h wallRMS wallMax interiorRMS interiorMax  p(wallRMS) p(interiorRMS)")
previous = None
for n in sorted({n for n, _ in shell}):
    r = shell[(n, "serial")]
    h = float(r["h"])
    line = f"{n} {h:.4g} {r['wallRMS']} {r['wallMax']} {r['interiorRMS']} {r['interiorMax']}"
    if previous:
        hp, wp, ip = previous
        ratio = math.log(hp/h)
        line += f"  {math.log(wp/float(r['wallRMS']))/ratio:.2f}"
        line += f" {math.log(ip/float(r['interiorRMS']))/ratio:.2f}"
    previous = (h, float(r["wallRMS"]), float(r["interiorRMS"]))
    print(line)

for message in failed:
    print("FAIL " + message, file=sys.stderr)
if failed:
    raise SystemExit(1)
print("PASS: stock reproduction, gradient names, affine exactness, serial = parallel")
PY
