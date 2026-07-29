#!/bin/bash
source /Volumes/OpenFOAM-v2412/etc/bashrc

cases="monodomain_coarse monodomain_coarse_tet"
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

for case in $cases; do
    echo "=========================================================="
    echo "Running cardiacFoam in parallel for $case..."
    echo "=========================================================="
    
    cd $script_dir/../simulations/$case
    
    # Clean previous run
    echo "Cleaning old parallel and time step data..."
    rm -f log.cardiacFoam log.decomposePar log.reconstructPar
    find . -maxdepth 1 -name "0.*" -type d -exec rm -rf {} +
    find . -maxdepth 1 -name "processor*" -type d -exec rm -rf {} +
    
    # Decompose
    echo "Decomposing domain into 6 subdomains..."
    decomposePar > log.decomposePar 2>&1
    if [ $? -ne 0 ]; then
        echo "decomposePar failed for $case! Check log.decomposePar"
        continue
    fi
    
    # Run simulation
    echo "Starting parallel simulation (endTime=0.4s)..."
    mpirun -np 6 cardiacFoam -parallel > log.cardiacFoam 2>&1
    if [ $? -ne 0 ]; then
        echo "cardiacFoam failed for $case! Check log.cardiacFoam"
        # We continue anyway to reconstruct whatever it managed to output
    else
        echo "Simulation finished successfully for $case."
    fi
    
    # Reconstruct
    echo "Reconstructing parallel results..."
    reconstructPar > log.reconstructPar 2>&1
    if [ $? -eq 0 ]; then
        echo "Removing processor directories..."
        find . -maxdepth 1 -name "processor*" -type d -exec rm -rf {} +
    else
        echo "reconstructPar failed! Check log.reconstructPar"
    fi
    
done

echo "All simulations have been processed."
