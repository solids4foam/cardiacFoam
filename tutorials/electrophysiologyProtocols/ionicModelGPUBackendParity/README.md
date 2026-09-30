# Ionic-model GPU backend parity

This is a small 2-D monodomain case for comparing a physiological ionic model
in its scalar, batched-host, and batched-CUDA forms. The tissue PDE remains on
the CPU; CUDA applies only to the batched ionic update.

The checked-in base case is a TNNP epicardial-cell slab. The study definition
in [`setup/studies/modelMatrix/sweep_scalar_batched.json`](setup/studies/modelMatrix/sweep_scalar_batched.json)
lists the 11 myocardial models, their tissue choices, and the matched scalar,
host-batched, and CUDA-batched execution modes. Fabbri is included in the
single-cell study only because it is an AV-node model, not a myocardial slab
model.

The JSON is deliberately declarative: it contains case inputs and comparison
settings, but no scheduler, machine, or orchestration-tool commands. Generated
cases, fields, logs, and result tables remain local and are not committed.

## Base case

The slab is 20 mm by 3 mm with 1,500 cells and a stimulus at one end. Its
baseline uses `deltaT = 2e-6 s`, five ionic substeps, and the batched
Rush--Larsen/Euler integrator. For a scalar reference, select the scalar
ionic-model name and remove `batchedIntegrator` and `batchedSubsteps` from
`constant/electroProperties`; keep the mesh, stimulus, tissue timestep, and
written times unchanged.

Run the checked-in base case directly with:

```bash
cd tutorials/electrophysiologyProtocols/ionicModelGPUBackendParity
./Allrun
```

With CUDA compiled into `libionicModels`, a batched model selects CUDA when a
device is visible to the solver process and otherwise uses its host path.
Confirm the device-selection line in `log.cardiacFoam` before treating a run
as a CUDA result.

## Post-processing

Compare two matching case directories after their outputs are written:

```bash
python3 tutorials/electrophysiologyProtocols/ionicModelGPUBackendParity/setup/postprocess_backend_parity.py \
    /path/to/reference-case /path/to/candidate-case \
    --times 0.005 0.01 0.015
```

It reports voltage error, activated-cell counts and masks, and activation-time
differences. It is appropriate for scalar-versus-batched or
host-batched-versus-CUDA comparisons. A moving wavefront can make direct
voltage error large for a small activation-time shift, so interpret the field
metrics together.
