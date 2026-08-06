#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
TET_TEMPLATE="$CASE_DIR/setup/mesh/tet/box.geo.template"
TET_SCHEMES="$CASE_DIR/setup/mesh/tet/fvSchemes"
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
PYTHON="${PYTHON:-python3}"
# tet   -- gmsh Delaunay meshes + this case's own setup/mesh/tet overlay
#          (leastSquares gradScheme; unit-cube geometry shared verbatim with
#          monodomainPseudoECG/eikonalECG's tet overlays).
# ortho -- the case's own blockMesh ladder and fvSchemes, i.e. the configuration
#          behind the reported Cartesian bidomain convergence table.
MESH_MODE="${MESH_MODE:-tet}"
RESOLUTIONS_STR="${RESOLUTIONS:-10 20 40}"
VARIANTS_STR="${VARIANTS:-baseline outer2 nonorth1 combined}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
WORK_ROOT="${WORK_ROOT:-$(mktemp -d /tmp/cardiacfoam-bidomain-corrector-study-XXXXXX)}"
KEEP_WORK="${KEEP_WORK:-0}"

read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"
read -r -a VARIANTS <<< "$VARIANTS_STR"

set +eu
source "$OPENFOAM_BASHRC" >/dev/null
# Darwin/SIP strips DYLD_LIBRARY_PATH across a fresh bash exec; RunFunctions
# restores it from FOAM_LD_LIBRARY_PATH (see its own "Darwin workaround" block).
# Without this, cardiacFoam aborts with "Library not loaded: @rpath/libOpenFOAM.dylib".
source "$WM_PROJECT_DIR/bin/tools/RunFunctions" >/dev/null 2>&1
set -eu

required=(cardiacFoam checkMesh foamDictionary)
case "$MESH_MODE" in
    tet)   required+=(gmsh gmshToFoam) ;;
    ortho) required+=(blockMesh) ;;
    *) echo "Unknown MESH_MODE '$MESH_MODE' (expected tet or ortho)" >&2; exit 2 ;;
esac

for executable in "${required[@]}"; do
    if ! command -v "$executable" >/dev/null 2>&1; then
        echo "Required executable not found: $executable" >&2
        exit 2
    fi
done

cleanup()
{
    if [[ "$KEEP_WORK" == "1" ]]; then
        echo "Keeping work directory: $WORK_ROOT"
    else
        rm -rf "$WORK_ROOT"
    fi
}
trap cleanup EXIT

rm -rf "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR"

dt_for_n()
{
    case "$1" in
        10) echo 0.00892857 ;;
        20) echo 0.00224215 ;;
        40) echo 0.000560538 ;;
        80) echo 0.000140174 ;;
        *) echo "No timestep configured for N=$1" >&2; exit 2 ;;
    esac
}

steps_for_n()
{
    case "$1" in
        10) echo 2 ;;
        20) echo 9 ;;
        40) echo 36 ;;
        *) echo "No step count configured for N=$1" >&2; exit 2 ;;
    esac
}

variant_controls()
{
    case "$1" in
        baseline) echo "1 0" ;;
        outer2)   echo "2 0" ;;
        outer3)   echo "3 0" ;;
        outer4)   echo "4 0" ;;
        outer8)   echo "8 0" ;;
        outer16)  echo "16 0" ;;
        nonorth1) echo "1 1" ;;
        combined) echo "2 1" ;;
        *) echo "Unknown variant '$1'" >&2; exit 2 ;;
    esac
}

clean_run_state()
{
    local run_case="$1"
    find "$run_case" -maxdepth 1 -type d -name '[0-9]*' -exec rm -rf {} +
    rm -rf "$run_case/postProcessing"
}

for n in "${RESOLUTIONS[@]}"; do
    run_case="$WORK_ROOT/N$n"
    mkdir -p "$run_case"
    rsync -a \
        --exclude '/[0-9]*/' \
        --exclude '/constant/polyMesh/' \
        --exclude '/logs/' \
        --exclude '/postProcessing/' \
        --exclude '/setup/mesh/tet/setup/corrector/results/' \
        "$CASE_DIR/" "$run_case/"

    if [[ "$MESH_MODE" == "tet" ]]; then
        lc="$($PYTHON -c "print(1.0/float('$n'))")"
        sed "s/__LC__/$lc/" "$TET_TEMPLATE" > "$run_case/setup/box.geo"
        cp "$TET_SCHEMES" "$run_case/system/fvSchemes"
        (
            cd "$run_case"
            gmsh -3 setup/box.geo -o box.msh -format msh2 > log.gmsh 2>&1
            gmshToFoam box.msh > log.gmshToFoam 2>&1
            rm box.msh
            checkMesh > log.checkMesh 2>&1
        )
    else
        # Keep the case's own fvSchemes: the reported Cartesian ladder runs
        # Gauss linear, and substituting leastSquares would not test it.
        sed -E "s/^(hex \(0 1 2 3 4 5 6 7\)) \([0-9]+ [0-9]+ [0-9]+\)/\1 ($n $n $n)/" \
            "$CASE_DIR/system/blockMeshDict.3D" > "$run_case/system/blockMeshDict.3D"
        if ! grep -q "($n $n $n)" "$run_case/system/blockMeshDict.3D"; then
            echo "Failed to set cell count $n in blockMeshDict.3D" >&2
            exit 2
        fi
        (
            cd "$run_case"
            blockMesh -dict system/blockMeshDict.3D > log.blockMesh 2>&1
            checkMesh > log.checkMesh 2>&1
        )
    fi
    cp "$run_case/log.checkMesh" "$RESULTS_DIR/N${n}_mesh.log"

    dt="$(dt_for_n "$n")"
    if [[ -n "${ENDTIME:-}" ]]; then
        # Reproduce a reported run: hold the physical window fixed and let the
        # step count follow dt, rather than using the short screening window.
        steps="$($PYTHON -c "import math;print(math.ceil(float('$ENDTIME')/float('$dt')))")"
    else
        steps="$(steps_for_n "$n")"
    fi
    end_time="$($PYTHON -c "print(float('$dt')*int('$steps'))")"
    foamDictionary "$run_case/system/controlDict" -entry deltaT -set "$dt" >/dev/null
    foamDictionary "$run_case/system/controlDict" -entry endTime -set "$end_time" >/dev/null
    foamDictionary "$run_case/system/controlDict" -entry writeControl -set timeStep >/dev/null
    foamDictionary "$run_case/system/controlDict" -entry writeInterval -set "$steps" >/dev/null
    foamDictionary "$run_case/system/controlDict" -entry writeFormat -set ascii >/dev/null

    for variant in "${VARIANTS[@]}"; do
        read -r outer nonorth <<< "$(variant_controls "$variant")"
        clean_run_state "$run_case"
        foamDictionary "$run_case/system/fvSolution" \
            -entry PIMPLE.nOuterCorrectors -set "$outer" >/dev/null
        if foamDictionary "$run_case/system/fvSolution" \
            -entry PIMPLE.nNonOrthogonalCorrectors >/dev/null 2>&1; then
            foamDictionary "$run_case/system/fvSolution" \
                -entry PIMPLE.nNonOrthogonalCorrectors -set "$nonorth" >/dev/null
        else
            foamDictionary "$run_case/system/fvSolution" \
                -entry PIMPLE.nNonOrthogonalCorrectors -add "$nonorth" >/dev/null
        fi

        run_id="N${n}_${variant}"
        out_dir="$RESULTS_DIR/$run_id"
        mkdir -p "$out_dir"
        printf 'run_id=%s\nresolution=%s\nvariant=%s\nn_outer=%s\nn_nonorth=%s\ndelta_t=%s\nn_steps=%s\nend_time=%s\n' \
            "$run_id" "$n" "$variant" "$outer" "$nonorth" \
            "$dt" "$steps" "$end_time" > "$out_dir/metadata.env"

        echo "=== $run_id ==="
        (cd "$run_case" && cardiacFoam > "$out_dir/log.cardiacFoam" 2>&1)
        dat="$(find "$run_case/postProcessing" -maxdepth 1 \
            -name '3D_*_cells_implicit.dat' -print | head -1)"
        if [[ -z "$dat" ]]; then
            echo "No bidomain MMS summary produced for $run_id" >&2
            exit 3
        fi
        cp "$dat" "$out_dir/summary.dat"
        cp "$run_case/system/fvSolution" "$out_dir/fvSolution"

        # Final coupled fields, needed to measure the iteration error against a
        # converged outer reference; the MMS norm alone cannot select a count.
        latest="$(find "$run_case" -maxdepth 1 -type d -name '[0-9]*' \
            | sed "s|.*/||" | sort -g | tail -1)"
        if [[ -z "$latest" || "$latest" == "0" ]]; then
            echo "No written time directory for $run_id" >&2
            exit 3
        fi
        for field in Vm phiE; do
            cp "$run_case/$latest/$field" "$out_dir/$field"
        done
        echo "latest_time=$latest" >> "$out_dir/metadata.env"
    done
done

"$PYTHON" "$SCRIPT_DIR/summarize_corrector_study.py" "$RESULTS_DIR"
echo "Bidomain corrector study complete: $RESULTS_DIR"

