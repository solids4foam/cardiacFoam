# electroDomains

The domains the core builds. Each one owns the state of one physical region and hands the numerics to a solver.

## What's available

| Domain | What it holds |
|---|---|
| `myocardiumDomain` | the heart tissue: `Vm`, its ionic model and a tissue solver (monodomain or bidomain). `eikonalMyocardiumDomain` is its eikonal form. |
| `conductionSystemDomain` | the Purkinje network as a graph, with its ionic model and a Purkinje solver |
| `ecgDomain` | electrodes and an ECG model; it reads the tissue state and never changes it |
| `extracellularPotentialDomain` | the global extracellular potential `phiE` over heart and bath, shared with a bidomain myocardium. Its heart+torso operator, `passiveConductor`, is also the conductor of the lead fields of `torsoECG`. |

The domains exchange data only through interfaces: coupling endpoints, used by the Purkinje–muscle junction couplers, and state providers, used by the ECG.

## Time stepping of `myocardiumDomain`

A step over `[t0, t0+dt]` is a reaction step, the ionic ODE, followed by the diffusion solve of `Vm`. Without a bath `advance` does both. With a bath the advance scheme calls them separately (`solveReactionStep`, the `phiE`/`Vm` solves, `finalizeDiffusionStep`), and both paths share the same reaction and end-of-step code.

`timeCouplingScheme` sets how the ionic current enters the diffusion solve:

- `godunov`: the current after the ODE step, held fixed during the diffusion solve. First order in time.
- `sbdf2`: the implicit diffusion solve takes the extrapolated current `I* = I^n + r (I^n - I^{n-1})`, with `r = dt/dt_old`. After the solve the current is evaluated again at the new `Vm`, and the ODE step receives the rate of change of `Vm` from the previous step. It requires `solutionAlgorithm implicit` and an ionic model that supports both couplings; any other combination stops at input. With `ddt(Vm) backward` it is second order.

`sbdf2` needs two steps of history. The first step of a run uses no rate hint and the first two use no extrapolation. These histories are not written. A restarted run therefore rebuilds them exactly as a run from `t = 0` does, and the log records the restart.

In the implicit solve, PIMPLE `nOuterCorrectors` and `nNonOrthogonalCorrectors` repeat the same operation: they reassemble the `Vm` equation with updated explicit terms and solve it again. Neither evaluates the ionic current again, so only their product sets the number of solves. With a bath, `nOuterCorrectors` sets the number of `phiE`/`Vm` passes instead.

## Folders

```text
src/electroModels/electroDomains/
├── myocardiumDomain/
├── conductionSystemDomain/
├── ecgDomain/
└── extracellularPotentialDomain/
```

The solvers these domains use are in [myocardiumModels](../myocardiumModels/README.md), [conductionSystemModels](../conductionSystemModels/README.md) and [ecgModels](../ecgModels/README.md); the couplers are in [electroCouplers](../electroCouplers/README.md).
