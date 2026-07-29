# Automation and Processing Scripts

Contains scripts for generating geometries, automating the simulation setup, and launching cases.

- `generate_mesh.py`: Gmsh python API script to generate anisotropic and isotropic tetrahedral meshes.
- `setup_all_cases.sh` / `setup_tet_cases.sh`: Automates mesh conversion (`gmshToFoam`, `polyDualMesh`) and `fiberFoam` coordinate generation.
- `run_simulations.sh`: Sequential execution script for all configured `cardiacFoam` simulations.
