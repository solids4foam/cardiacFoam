#!/bin/bash


# Source OpenFOAM directly
source /Volumes/OpenFOAM-v2412/etc/bashrc

CWD=$(pwd)

echo "Starting CV Convergence Test for dx = 0.5, 0.2, 0.1, 0.05 mm"
echo "--------------------------------------------------------"

for cells in 40 100 200 400; do
    dx=$(echo "scale=3; 20/$cells" | bc)
    echo "Running dx = $dx mm ($cells cells)"
    
    # Update blockMeshDict
    sed -i '' "s/hex (0 1 2 3 4 5 6 7) ([0-9]* 1 1)/hex (0 1 2 3 4 5 6 7) ($cells 1 1)/" system/blockMeshDict
    
    # Run
    ./Allclean >/dev/null 2>&1 || true
    blockMesh > log.blockMesh 2>&1 || true
    cardiacFoam > log.cardiacFoam 2>&1 || true
    python3 setup/extract_cv.py > log.extract 2>&1 || true
    
    if [ -f postProcessing/cv_summary.txt ]; then
        cv=$(cat postProcessing/cv_summary.txt | grep 'Central calibration CV:' -A 1 | tail -n 1 | grep -Eo 'CV = [0-9]+\.[0-9]+' | grep -Eo '[0-9]+\.[0-9]+')
        echo "  Result: CV = $cv m/s"
    else
        echo "  Result: FAILED to compute CV"
    fi
done

echo "--------------------------------------------------------"
echo "Done."
