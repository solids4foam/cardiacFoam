# Small 2D ionic coupling slab

This is a 20 mm by 3 mm TNNP monodomain slab with one empty mesh layer in
the third direction (1,500 cells). It is a simple propagation case, not a
reproduction of the published Niederer benchmark. The stimulus occupies the
first 1.5 mm in x and y and lasts 2 ms. The supplied baseline uses
`TNNPcompactBatched`, `rushLarsen`, `batchedSubsteps 5`, and a tissue
`deltaT` of 2 microseconds. Its effective ionic substep is 0.4 microseconds.

With OpenFOAM v2412 loaded and `cardiacFoam` built, run:

```bash
cd tutorials/ioniGPUKernel/slab2D
./Allrun
```

The case writes `Vm`, `ionicCurrent`, and `activationTime` at 5, 10, and
15 ms. To compare with the scalar reference, set `ionicModel TNNP;` and
remove the two `batched*` settings in `constant/electroProperties`, then run
in a separate copy of this case. Keep mesh, stimulus, tissue, `deltaT`, and
output times identical.

The parent directory has a case generator for all 12 models in the
single-cell inventory. It copies this geometry and stimulus to new
directories, and selects the tissue type and time step used in that
inventory. Generate an individual case from the repository root:

```bash
python3 tutorials/ioniGPUKernel/prepare_slab_case.py TWorld scalar /tmp/tworld-scalar
python3 tutorials/ioniGPUKernel/prepare_slab_case.py TWorld batched /tmp/tworld-batched
(cd /tmp/tworld-scalar && ./Allrun)
(cd /tmp/tworld-batched && ./Allrun)
python3 tutorials/ioniGPUKernel/slab2D/setup/post_processing_slab_field_parity.py \
    /tmp/tworld-scalar /tmp/tworld-batched
```

Replace `TWorld` with any model name listed by the script. The short paired
matrix can also create and run all 12 cases:

```bash
python3 tutorials/ioniGPUKernel/run_slab_matrix.py \
    /tmp/cardiac-slab-matrix --end-time 0.005 --workers 2 --require-gpu
```

It writes `summary.csv`, plus logs and fields under the output directory.
Use a fresh output path for every invocation. Default steps are 20
microseconds for BuenoOrovio and 2 microseconds for the other models, with
five ionic substeps for batched cases. The 5 ms default is an immediate
failure and coupling screen; increase `--end-time` for a propagation run.
These cases use a common external stimulus and are not a substitute for
model-specific pacing and validation. They do not reproduce the published
Niederer benchmark.

## Integrator and time-step study

The parent tutorial provides a sweep for Euler and Rush–Larsen plus Euler. It
compares each batched GPU run with scalar RKF45 at the same tissue `deltaT`, a
fine scalar case, and a fine Rush–Larsen GPU baseline. For a quick TNNP study:

```bash
python3 tutorials/ioniGPUKernel/run_slab_integrator_study.py \
    /tmp/tnnp-integrator-sweep \
    --models TNNP --integrators euler rushLarsen \
    --dt-multipliers 1 2 5 --steps 1 5 10 25 --require-gpu
python3 tutorials/ioniGPUKernel/slab2D/setup/post_processing_integrator_sweep_summary.py \
    /tmp/tnnp-integrator-sweep
```

The ionic substep is `deltaT / batchedSubsteps`. Hold `deltaT` fixed while
increasing substeps to test the ionic integrator at a fixed PDE step. Increase
both proportionally to keep the ionic substep approximately fixed while
testing a coarser PDE step. Increasing `deltaT` alone changes both. Failed
cases remain in the sweep output and should inform the stable range, rather
than being omitted from the report.

The sweep summary is an early screen based on activated-cell mask agreement
and activation-time differences; it is not a full APD/state/current
validation. The post-processing scripts and their metric definitions are in
[`setup`](setup/):

```bash
python3 tutorials/ioniGPUKernel/slab2D/setup/post_processing_slab_field_parity.py \
    /path/to/reference-case /path/to/candidate-case
```

The scripts in `setup` are copies of the comparison and sweep-summary tools
used by the parent study. The runners remain in `tutorials/ioniGPUKernel` so
they can also be used with the model matrix.

For a CUDA build, use the architecture of the target GPU in `NVARCH`:

```bash
source /usr/lib/openfoam/openfoam2412/etc/bashrc
export CARDIAC_ENABLE_CUDA=1
export CUDA_HOME=/usr
export PATH=/usr/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
export NVARCH=75
export CUDA_HOST_CXX=/usr/bin/g++-10
(cd src/ionicModels && wclean libso)
(cd src/ionicModels && wmake libso)
```

This command reflects the tested CUDA 11.5 build on xenosim; its RTX 4000 Ada
run used `NVARCH=75` and GCC 10 for nvcc. Keep `FOAM_SIGFPE=false` in the GPU
job environment. The login node has no visible GPU, so run the case in a
Slurm GPU allocation. Add `--require-gpu` to the matrix runner and confirm the
case log reports a CUDA device; a completed CPU fallback is not a GPU pass.

To run a longer GPU propagation comparison, use a fresh output directory and
extend the end time:

```bash
python3 tutorials/ioniGPUKernel/run_slab_matrix.py \
    /tmp/cardiac-slab-15ms --models TNNP TWorld BuenoOrovio \
    --end-time 0.015 --workers 1 --require-gpu
```

This remains a simple slab test. It does not reproduce the published
Niederer benchmark conditions.
