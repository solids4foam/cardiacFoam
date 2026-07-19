#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
OPENFOAM_BASHRC="${OPENFOAM_BASHRC:-/Volumes/OpenFOAM-v2412/etc/bashrc}"
PYTHON="${PYTHON:-python3}"
RESOLUTIONS_STR="${RESOLUTIONS:-10 20}"
VARIANTS_STR="${VARIANTS:-baseline predictor phi8}"
RESULTS_DIR="${RESULTS_DIR:-$SCRIPT_DIR/results}"
WORK_ROOT="${WORK_ROOT:-$(mktemp -d /tmp/cardiacfoam-bath-coupling-study-XXXXXX)}"
KEEP_WORK="${KEEP_WORK:-0}"
NPROCS="${NPROCS:-6}"

read -r -a RESOLUTIONS <<< "$RESOLUTIONS_STR"
read -r -a VARIANTS <<< "$VARIANTS_STR"

if [[ ! -f "$OPENFOAM_BASHRC" ]]; then
    echo "OpenFOAM bashrc not found: $OPENFOAM_BASHRC" >&2
    exit 2
fi

set +eu
source "$OPENFOAM_BASHRC" >/dev/null
set -eu

required_executables=(cardiacFoam bathBidomainInterfaceMetrics \
    setTorsoOrganConductivityField gmsh gmshToFoam foamDictionary)
if [[ "$NPROCS" -gt 1 ]]; then
    required_executables+=(decomposePar mpirun reconstructPar)
fi
for executable in "${required_executables[@]}"
do
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
        baseline)  echo "onePass 0" ;;
        predictor) echo "predictorCorrector 0" ;;
        phi8)      echo "onePass 8" ;;
        *) echo "Unknown variant '$1'" >&2; exit 2 ;;
    esac
}

clean_time_state()
{
    local run_case="$1"
    find "$run_case" -maxdepth 1 -type d \
        \( -name '0' -o -name '[0-9]*' -o -name 'postProcessing' \
           -o -name 'processor*' \) \
        -exec rm -rf {} +
}

for n in "${RESOLUTIONS[@]}"
do
    run_case="$WORK_ROOT/N$n"
    mkdir -p "$run_case"
    rsync -a \
        --exclude '/[0-9]*/' \
        --exclude '/postProcessing/' \
        --exclude '/processor*/' \
        --exclude '/constant/polyMesh/' \
        --exclude '/setup/couplingStudy/results/' \
        --exclude '/setup/interfaceStudy/' \
        --exclude '/setup/interfaceMeshBank/' \
        --exclude '/setup/matchedSubmeshStudy/' \
        --exclude '/setup/results/' \
        "$CASE_DIR/" "$run_case/"

    lc="$($PYTHON -c "print(1.0/float('$n'))")"
    sed "s/__LC__/$lc/" "$CASE_DIR/setup/three_domain_box.geo.template" \
        > "$run_case/setup/three_domain_box.geo"
    (
        cd "$run_case"
        gmsh -3 setup/three_domain_box.geo -o three_domain_box.msh -format msh2 \
            > log.gmsh 2>&1
        gmshToFoam three_domain_box.msh > log.gmshToFoam 2>&1
        rm -f three_domain_box.msh
        checkMesh > log.checkMesh 2>&1
        "$PYTHON" setup/verify_mesh.py "$run_case" > log.verifyMesh 2>&1
    )
    cp "$run_case/log.checkMesh" "$RESULTS_DIR/N${n}_mesh.log"

    dt="$(dt_for_n "$n")"
    steps="$(steps_for_n "$n")"
    foamDictionary "$run_case/system/controlDict" -entry deltaT -set "$dt" >/dev/null
    foamDictionary "$run_case/system/controlDict" -entry endTime -set 0.02 >/dev/null
    foamDictionary "$run_case/system/controlDict" -entry writeControl -set timeStep >/dev/null
    foamDictionary "$run_case/system/controlDict" -entry writeInterval -set "$steps" >/dev/null
    foamDictionary "$run_case/system/fvSolution" \
        -entry PIMPLE.nOuterCorrectors -set 1 >/dev/null
    foamDictionary "$run_case/system/fvSolution" \
        -entry PIMPLE.nNonOrthogonalCorrectors -set 1 >/dev/null
    foamDictionary "$run_case/constant/electroProperties" \
        -entry bidomainSolverCoeffs.bathPotentialDomain.interfaceConductivityInterpolation \
        -set distanceWeightedHarmonic >/dev/null
    foamDictionary "$run_case/constant/electroProperties" \
        -entry bidomainSolverCoeffs.bathPotentialDomain.intracellularAssembly \
        -set matchedSubmesh >/dev/null
    if [[ "$NPROCS" -gt 1 ]]; then
        foamDictionary "$run_case/system/decomposeParDict" \
            -entry numberOfSubdomains -set "$NPROCS" >/dev/null
    fi

    for variant in "${VARIANTS[@]}"
    do
        read -r method phi_nonorth <<< "$(variant_controls "$variant")"
        run_id="N${n}_${variant}"
        out_dir="$RESULTS_DIR/$run_id"
        mkdir -p "$out_dir"
        clean_time_state "$run_case"

        foamDictionary "$run_case/constant/electroProperties" \
            -entry bidomainSolverCoeffs.bathPdeCouplingMethod \
            -set "$method" >/dev/null
        foamDictionary "$run_case/constant/electroProperties" \
            -entry bidomainSolverCoeffs.bathPotentialDomain.phiENonOrthogonalCorrectors \
            -set "$phi_nonorth" >/dev/null

        printf 'run_id=%s\nresolution=%s\nvariant=%s\nmethod=%s\nphi_nonorth=%s\ndelta_t=%s\n' \
            "$run_id" "$n" "$variant" "$method" \
            "$phi_nonorth" "$dt" > "$out_dir/metadata.env"

        echo "=== $run_id (nprocs=$NPROCS) ==="
        (
            cd "$run_case"
            setTorsoOrganConductivityField > "$out_dir/log.setConductivity" 2>&1
            if [[ "$NPROCS" -gt 1 ]]; then
                decomposePar -force > "$out_dir/log.decomposePar" 2>&1
                mpirun -np "$NPROCS" cardiacFoam -parallel \
                    > "$out_dir/log.cardiacFoam" 2>&1
                reconstructPar -latestTime > "$out_dir/log.reconstructPar" 2>&1
            else
                cardiacFoam > "$out_dir/log.cardiacFoam" 2>&1
            fi
            bathBidomainInterfaceMetrics -latestTime \
                > "$out_dir/log.interfaceMetrics" 2>&1
        )
        cp "$run_case/postProcessing/bathBidomainInterfaceMetrics.csv" \
            "$out_dir/metrics.csv"
        cp "$run_case/constant/electroProperties" "$out_dir/electroProperties"
        cp "$run_case/system/fvSolution" "$out_dir/fvSolution"
    done
done

"$PYTHON" "$SCRIPT_DIR/summarize_coupling_study.py" "$RESULTS_DIR"
echo "Bath coupling study complete: $RESULTS_DIR"
