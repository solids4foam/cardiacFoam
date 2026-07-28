#!/bin/bash

set -e

# Scenarios defined based on gaur_xenotransplant_test_scenarios.md mapping:
# GCaL -> AC_PCa
# Jrel -> AC_grelbarjsrol
# INaCa -> AC_Gncx
# INaK -> AC_Pnak

scenarios=(
    "0_Baseline:::"
    "1_PrimaryCaReleaseFailure:AC_PCa:0.95:AC_grelbarjsrol:0.60"
    "2_SERCARelaxationDysfunction:AC_grelbarjsrol:0.90:AC_Gncx:1.10"
    "3_MyofilamentContractilityLoss:::" # We can't scale myofilaments without contractility variable, so leaving empty or scaling something else if user specified. The doc says GCaL=1, Jrel=1.
    "4_RejectionCaInjury:AC_PCa:0.80:AC_grelbarjsrol:0.65:AC_Pnak:0.80:AC_Gncx:1.20:AC_GKr:0.85:AC_GK1:0.90"
)

# Template for electroProperties
TEMPLATE=$(cat << 'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      electroProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
myocardiumSolver singleCellSolver;
singleCellSolverCoeffs
{
    // Cell model selection
    ionicModel    Gaur;
    tissue        myocyte;

    ionicConstantOverrides
    {
        myocyte
        {
            scale
            {
__OVERRIDES__
            }
        }
    }

    // Time integration settings
    solutionAlgorithm explicit;
    solver          RKF45;
    initialODEStep  1e-5;
    maxSteps        1000000000;
    absTol          1e-6;
    relTol          1e-4;

    // Stimulus protocol
    singleCellStimulus
    {
        stim_start      20;
        stim_duration   1;
        stim_amplitude  60;
        stim_period_S1  1000;
        nstim1          10;
        stim_period_S2  0;
        nstim2          0;
    }

    // optional I/O control
    writeAfterTime  8;
    writeFrequency  0.001;

    outputVariables
    {
        ionic
        {
            export (Vm cai cansr cajsr AV_ICaL_ICaL AV_INaCa AV_Jrel AV_Jup AV_Jleak);
            debug (Vm);
        }
        activeTension
        {
            export (AV_Ta Ca_TRPN TmBlocked XW XS);
        }
    }

    activeTensionModel LandNiederer;
}

EOF
)

for scenario in "${scenarios[@]}"; do
    # parse scenario
    IFS=':' read -r name v1 s1 v2 s2 v3 s3 v4 s4 v5 s5 v6 s6 <<< "$scenario"
    echo "Running Scenario: $name"

    overrides=""
    if [ -n "$v1" ]; then overrides="${overrides}${v1} ${s1};"; fi
    if [ -n "$v2" ]; then overrides="${overrides}
                ${v2} ${s2};"; fi
    if [ -n "$v3" ]; then overrides="${overrides}
                ${v3} ${s3};"; fi
    if [ -n "$v4" ]; then overrides="${overrides}
                ${v4} ${s4};"; fi
    if [ -n "$v5" ]; then overrides="${overrides}
                ${v5} ${s5};"; fi
    if [ -n "$v6" ]; then overrides="${overrides}
                ${v6} ${s6};"; fi

    # generate dict
    echo "${TEMPLATE//__OVERRIDES__/$overrides}" > constant/electroProperties

    # run
    rm -f log.cardiacFoam
    export CF_SKIP_PLOTS=1
    ./Allrun
    if [ -f log.cardiacFoam ]; then
        mv log.cardiacFoam "log.$name"
    fi
done

echo "All simulations completed."
