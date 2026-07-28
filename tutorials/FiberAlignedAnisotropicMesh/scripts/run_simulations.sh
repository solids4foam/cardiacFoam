#!/bin/bash
source /Volumes/OpenFOAM-v2412/etc/bashrc

cases="monodomain_coarse monodomain_coarse_tet monodomain_medium monodomain_medium_tet monodomain_fine monodomain_fine_tet"
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

for case in $cases; do
    echo "=========================================================="
    echo "Running cardiacFoam for $case..."
    echo "=========================================================="
    
    cd $script_dir/../simulations/$case
    
    # Run cardiacFoam and output to log.cardiacFoam
    cardiacFoam > log.cardiacFoam 2>&1 &
    PID=$!
    wait $PID
    
    # Check if run was successful
    if [ $? -eq 0 ]; then
        echo "Simulation finished successfully for $case."
    else
        echo "Simulation failed for $case! Check log.cardiacFoam."
    fi
done

echo "All 6 simulations completed!"
