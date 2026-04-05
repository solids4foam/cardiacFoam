# domains

This directory holds domain implementations used by the electro core.

Current contents:

- `conductionSystemDomain/` owns pre-primary auxiliary conduction-network
  domain models and graph solvers.
- `ecgDomain/` owns agnostic post-primary ECG domain contracts and solver
  interfaces (`ECGDomain`, `ecgSolver`).
- `reactionDiffusionMyocardiumDomain/` owns the shared tissue-domain
  lifecycle, state, activation tracking, stimulus handling, and post-processing.
  The monodomain and bidomain PDE kernels now live under `solvers/`.

Concrete ECG domain model implementations live under
`src/electroModels/ecgModels/`.
