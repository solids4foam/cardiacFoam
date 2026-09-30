# Single-cell electrophysiology

This is the zero-dimensional `singleCellSolver` protocol. It evolves one ionic
model without a spatial PDE and writes traces below `postProcessing/`.

## Base case

`constant/electroProperties` selects the ionic model, tissue type, stimulus,
and ODE controls. The checked-in case can be run directly:

```bash
cd tutorials/electrophysiologyProtocols/singleCell
./Allrun
./regression/regressionTest.sh
```

## Studies

Each multi-case protocol lives under `setup/studies/` and contains declarative
JSON inputs plus only the post-processing specific to that study. Generated
traces and reports are local, ignored results.

- [`tworldVsGaur`](setup/studies/tworldVsGaur/README.md) compares focused pig
  and human ventricular pacing responses.
- [`ionicModelGPUBackendParity`](setup/studies/ionicModelGPUBackendParity/README.md)
  is the scalar, host-batched, and CUDA-batched model/tissue matrix. It is the
  0-D companion to the 2-D
  [`ionicModelGPUBackendParity`](../ionicModelGPUBackendParity/README.md)
  case.

The study JSON selects the model/tissue rows. A batched model uses its host
implementation without a visible CUDA device and its CUDA backend when CUDA is
built and a device is visible; this is runtime selection, not a different
model name.
