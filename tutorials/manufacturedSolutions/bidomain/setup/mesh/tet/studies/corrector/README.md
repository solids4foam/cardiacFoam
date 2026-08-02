# corrector - bidomain

## Purpose
This study validates the predictor-corrector inner loop convergence and stability for the bidomain solver.

## Execution
Run `./run_corrector_study.sh` and parse metrics with `./summarize_corrector_study.py`.

## Tracking & Outputs
All generated outputs, mesh files, and metric archives are saved to the local `results/` folder, which is explicitly ignored by git. Do not commit generated OpenFOAM data.
