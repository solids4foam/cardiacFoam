# electroModels

The electrophysiology part of the core: the core itself, the domains it builds (myocardium, Purkinje network, ECG, bath), their solvers, and the couplers between them. The top-level model, `electroModel`, is what the solver advances.

## What's available

| Piece | Names |
|---|---|
| Tissue solvers | `monodomainSolver`, `bidomainSolver`, `eikonalSolver` |
| Single cell | `singleCellSolver` |
| Purkinje solvers | `monodomain1DSolver`, `eikonalSolver1D`, `restitutionEikonalSolver1D` |
| ECG models | `pseudoECG` (class `pseudoECGSolver`), `torsoECG`, `eikonalECG` |
| Bath | `extracellularPotentialDomain`: a global `phiE` over heart and bath |
| Purkinje–muscle junction couplers | `eikonalPvjCoupler`, `eikonalMonodomainPvjCoupler`, `reactionDiffusionPvjCoupler` |
| Advance scheme | `staggeredElectrophysicsAdvanceScheme` |

## Folders

| Folder | What it holds |
|---|---|
| [core/](core/README.md) | builds the system of domains and advances it |
| [electroDomains/](electroDomains/README.md) | the domains: myocardium, Purkinje, ECG, bath |
| [myocardiumModels/](myocardiumModels/README.md) | the tissue solvers and the single-cell solver |
| [conductionSystemModels/](conductionSystemModels/README.md) | the Purkinje solvers |
| [ecgModels/](ecgModels/README.md) | the ECG models |
| [electroCouplers/](electroCouplers/README.md) | the couplers between domains |

**Deep dive:** [ARCHITECTURE.md](./ARCHITECTURE.md) and [core/ARCHITECTURE.md](./core/ARCHITECTURE.md) explain how this library is built inside.

## What this does not own

- The cell models: [ionicModels](../ionicModels/README.md).
- The concrete verifiers: [verificationModels](../verificationModels/README.md). Their abstract bases are here, in `core/verificationModels/`.
- Mechanics: [electroMechanicalModels](../electroMechanicalModels/README.md).
- The `Vm`/`Cai` contract: [couplingModels](../couplingModels/README.md).
