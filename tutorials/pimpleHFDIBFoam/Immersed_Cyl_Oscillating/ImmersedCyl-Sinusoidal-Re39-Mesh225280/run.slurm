#!/bin/bash -l
#SBATCH --job-name=Immersed-Osc-Sin-Re39-Mesh225280
#SBATCH --output=log.%x_%j.out
#SBATCH --error=log.%x_%j.err
#SBATCH --partition=main
#SBATCH --time=48:00:00
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1

# Good practice to navigate to the submission directory
cd $SLURM_SUBMIT_DIR

# Load OpenFOAM
#source /usr/lib/openfoam/openfoam8/etc/bashrc
#source /usr/lib/openfoam/openfoam2312/etc/bashrc
#source /usr/lib/openfoam/openfoam2212/etc/bashrc
#source /usr/lib/openfoam/openfoam2112/etc/bashrc

# Execute script
./Allrun-parallel
