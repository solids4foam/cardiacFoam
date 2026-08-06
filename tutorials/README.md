# tutorials architecture

This folder contains maintained reference, protocol, and verification cases for
`cardiacFoam`. The table below is the canonical index of committed runnable
cases; local research cases and generated sweep directories are not part of
this documented contract.

## Canonical tutorial cases

`Driver entry` is the strict `foamctl` / `openfoam_driver` registry name. A dash
means that the folder can still be addressed as a generic case but has no
dedicated registered spec. `Regression` means coverage by the cross-case
`Alltest-regression` runner, not merely the presence of an `Allrun` script.

| Canonical path | Purpose | Solver or workflow | Build mode | Driver entry | Regression | Principal outputs |
| --- | --- | --- | --- | --- | --- | --- |
| `electrophysiologyProtocols/singleCell` | Single-point action-potential runs and ionic-model sweeps | `singleCellSolver` | lightweight or full | `singleCell` | `Alltest-regression` | voltage traces under `postProcessing/` |
| `electrophysiologyProtocols/ionicHeterogeneityProbe` | Probe transmural ionic heterogeneity without a tissue PDE | Bueno-Orovio probe and plotting workflow | lightweight or full | — | not covered | trace, AP-metric, and smoothness CSV files under `postProcessing/ionicHeterogeneityProbe/` |
| `electrophysiologyProtocols/restitutionCurves_s1s2Protocol` | Generate S1-S2 action-potential-duration restitution curves | `singleCellSolver` pacing sweep | lightweight or full | `restitutionCurves` | not covered | per-run traces and restitution tables/plots |
| `electrophysiologyProtocols/cableProtocol/monodomain1DCableCV` | Measure and refine 1D conduction velocity | `monodomainSolver` cable convergence workflow | lightweight or full | `monodomainAndEikonal1DCableCVConvergence` | not covered | activation probes, CV summaries, convergence CSV files and plots |
| `electrophysiologyProtocols/cableProtocol/eikonal1DCableCV` | Compare eikonal 1D conduction velocity across resolutions | `eikonalSolver` cable workflow | lightweight or full | — | not covered | activation probes and CV summaries |
| `electrophysiologyProtocols/rotorInstability` | Exercise sustained re-entry and activation-time behavior | monodomain rotor protocol | lightweight or full | — | `Alltest-regression` | probe traces and activation-time metrics |
| `NiedererEtAl2011/NiedererEtAl2011verification` | Run the Niederer slab electrophysiology benchmark | `monodomainSolver` verification workflow | lightweight or full | `niederer2012` | `Alltest-regression` | activation probes, smoke-check fields, summaries and plots |
| `NiedererEtAl2011/purkinjeNiedererEtAl2011` | Couple a small 1D Purkinje graph to the Niederer slab | `monodomainSolver`, `monodomain1DSolver`, and PVJ coupling | lightweight or full | — | `Alltest-regression` | Purkinje graph data/VTK and PVJ activation-time checks |
| `NiedererEtAl2011/electroMechanicalNiedererEtAl2011` | Demonstrate sequential electrophysiology-solid coupling | `electroMechanicalModel` with monodomain electrophysiology | full only | — | `Alltest-regression` (expected skip in lightweight mode) | active-tension probes and coupled solid/electro fields |
| `manufacturedSolutions/monodomainPseudoECG` | Verify monodomain fields and pseudo-ECG against manufactured solutions | `monodomainFDAManufactured` plus pseudo-ECG verification | lightweight or full | `manufacturedFDA` | `Alltest-regression` | manufactured error summaries and pseudo-ECG series under `postProcessing/` |
| `manufacturedSolutions/bidomain` | Verify spatial bidomain convergence | `bidomainFDAManufactured` | lightweight or full | `manufacturedFDABidomain` | `Alltest-regression` | manufactured field-error summaries under `postProcessing/` |
| `manufacturedSolutions/bathBidomain` | Verify bidomain-with-bath fields and ECG ownership | `bathBidomainFDAManufactured` and optional `torsoECG` | lightweight or full | `manufacturedFDABathBidomain` | `Alltest-regression` | global bath fields and manufactured error/ECG summaries |
| `manufacturedSolutions/eikonalECG` | Verify activation time and template/quadrature ECG calculations | `eikonalSolver` with manufactured eikonal verification | lightweight or full | `manufacturedEikonalECG` | `Alltest-regression` | activation-time and ECG reference/error series plus summary CSV files |
| `manufacturedSolutions/monodomainTotalLagrangianEM` | Verify coupled monodomain and nonlinear solid mechanics | manufactured total-Lagrangian electromechanics workflow | full only | `manufacturedMonodomainTotalLagrangianEM` | not covered | `Vm`, `D`, `lambda`, and `Ta` error/convergence tables and plots |

## Common script pattern

Most runnable cases provide:

- `Allrun`: run the simulation and, where configured, post-process it;
- `Allclean`: remove generated output;
- optional `regressionTest.sh`: perform case-local quantitative checks.

Read the case-local `README.md` before running a case. It defines its inputs,
parallel options, expected cost, and output details.

## Cross-case regression entrypoint

The aggregate runner copies each covered case outside the source tutorial tree
before executing its local regression script. Select the build mode explicitly:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight ./tutorials/Alltest-regression
CARDIAC_REGRESSION_BUILD_MODE=with-solids4foam ./tutorials/Alltest-regression
```

The electromechanical Niederer regression is the only expected skip in
lightweight mode. Any other exit-77 skip fails the aggregate run.

## Python automation integration

Registered driver specs reuse the same case roots, refinement patterns,
dictionary mutations, output collection, and post-processing entrypoints shown
above. Generic-case discovery can address other runnable folders, but discovery
does not make a folder a dedicated registered tutorial or add it to the
cross-case regression suite.
