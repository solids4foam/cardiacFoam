#!/bin/bash
# setup_tet_cases.sh
# Generates OpenFOAM mesh from gmsh, scales it, maps MATLAB Cobiveco fields, and runs the solver

source /Volumes/OpenFOAM-v2412/etc/bashrc
export PATH=/Users/simaocastro/cobivecco-OpenFOam/platforms/darwin64ClangDPInt32Opt/bin:$PATH

# Canonical physical scale for BOTH pipelines (must match setup_all_cases.sh).
SCALE="(0.02 0.02 0.02)"

cd tutorials/FiberAlignedAnisotropicMesh

# Ensure VTU output is generated for MATLAB
echo "Generating VTU inputs for MATLAB Cobiveco..."
cd meshes
python3 ../scripts/meshio_to_cobiveco.py biv_ellipsoid.msh monodomain_coarse_tet
python3 ../scripts/meshio_to_cobiveco.py biv_ellipsoid.msh monodomain_medium_tet
python3 ../scripts/meshio_to_cobiveco.py biv_ellipsoid_fine.msh monodomain_fine_tet
cd ..

# Run MATLAB Cobiveco
echo "Running MATLAB Cobiveco for all tet meshes..."
cd ../..
/Applications/MATLAB_R2026a.app/bin/matlab -batch "run_batch_cobiveco" > tutorials/FiberAlignedAnisotropicMesh/log.matlab 2>&1
cd tutorials/FiberAlignedAnisotropicMesh

# Now setup the OpenFOAM cases
for res in coarse medium fine; do
    echo "Setting up monodomain_${res}_tet..."
    
    rm -rf simulations/monodomain_${res}_tet
    mkdir -p simulations/monodomain_${res}_tet
    
    cp -r templates/monodomainHeartTissue/* simulations/monodomain_${res}_tet/
    cp templates/monodomainHeartTissue/system/controlDict simulations/monodomain_${res}_tet/system/
    
    cd simulations/monodomain_${res}_tet
    
    if [ "$res" == "fine" ]; then
        cp ../../meshes/biv_ellipsoid_fine.msh .
        gmshToFoam biv_ellipsoid_fine.msh > log.gmshToFoam 2>&1
    else
        cp ../../meshes/biv_ellipsoid.msh .
        gmshToFoam biv_ellipsoid.msh > log.gmshToFoam 2>&1
    fi
    
    # No early scaling: generate content on the native mesh, scale last.
    # Write cell centres for VTU mapping (native coordinates).
    writeCellCentres > log.cellCentres 2>&1
    
    # Run the VTU to OpenFOAM converter
    # Note: MATLAB Cobiveco UVC output usually doesn't have fibers directly.
    # We first run our Python script to calculate the fibers on the MATLAB VTU
    echo "Calculating fibers on MATLAB VTU..."
    python3 ../../../scripts/calc_gold_fibers.py ../../meshes/result_${res}/monodomain_${res}_tet_cobiveco_result.vtu ../../meshes/result_${res}/monodomain_${res}_tet_with_fibers.vtu
    
    echo "Mapping MATLAB fields into OpenFOAM..."
    # VTU (cobiveco output) is at native units, matching the un-scaled OpenFOAM mesh.
    python3 ../../../scripts/vtu_to_openfoam.py ../../meshes/result_${res}/monodomain_${res}_tet_with_fibers.vtu 0 1.0
    
    # Set conductivity
    setCardiacConductivity > log.conductivity 2>&1

    # Scale-last: single transformPoints as the final geometry step.
    transformPoints -scale "$SCALE" > log.transform 2>&1

    touch results.foam
    
    cd ../..
done

echo "Starting cardiacFoam on coarse tet mesh (for quick validation)..."
cd simulations/monodomain_coarse_tet
cardiacFoam > log.cardiacFoam 2>&1 &
cd ../..

echo "Done!"
