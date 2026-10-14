# tutorials architecture

This folder contains maintained reference, cardiac protocols, and verification cases for
`cardiacFoam`. The table below is the canonical index of committed runnable
cases; local research cases and generated sweep directories are not part of
this documented contract.

## Folder layout

```text
tutorials/
├── electrophysiologyProtocols/    small EP checks: cells, restitution, CV, re-entry, Purkinje
├── electromechanicsProtocols/     small electromechanics checks: supports, active tension
├── NiedererEtAl2011verification/  Niederer et al. (2011) monodomain slab benchmark
├── manufacturedSolutions/         verification against known (manufactured) solutions
└── idealizedHeart/                anatomical biventricular cases (EP, EM, pathologies)
```

## Canonical tutorial cases

`Regression` means coverage by the cross-case `Alltest-regression` runner,
not merely the presence of an `Allrun` script. The runner discovers every
case with a `regression/regressionTest.sh`, and `regression/checkIndex.sh`
fails when this column and that discovery disagree.

| Canonical path | Purpose | Solver or workflow | Build mode | Regression | Principal outputs |
| --- | --- | --- | --- | --- | --- |
| `electrophysiologyProtocols/singleCell` | Single-point action-potential runs and ionic-model sweeps | `singleCellSolver` | lightweight or full | `Alltest-regression` | voltage traces under `postProcessing/` |
| `electrophysiologyProtocols/ionicHeterogeneity` | Probe transmural ionic heterogeneity without a tissue PDE | Bueno-Orovio probe and plotting workflow | lightweight or full | `Alltest-regression` | trace, AP-metric, and smoothness CSV files under `postProcessing/ionicHeterogeneityProbe/` |
| `electrophysiologyProtocols/restitutionCurves_s1s2Protocol` | Generate S1-S2 action-potential-duration restitution curves | `singleCellSolver` pacing sweep | lightweight or full | not covered | per-run traces and restitution tables/plots |
| `electrophysiologyProtocols/cableProtocol/monodomain1DCableCV` | Measure and refine 1D conduction velocity | `monodomainSolver` cable convergence workflow | lightweight or full | not covered | activation probes, CV summaries, convergence CSV files and plots |
| `electrophysiologyProtocols/cableProtocol/eikonal1DCableCV` | Compare eikonal 1D conduction velocity across resolutions | `eikonalSolver` cable workflow | lightweight or full | not covered | activation probes and CV summaries |
| `electrophysiologyProtocols/rotorInstability` | Exercise sustained re-entry and activation-time behavior | monodomain rotor protocol | lightweight or full | `Alltest-regression` | probe traces and activation-time metrics |
| `electrophysiologyProtocols/purkinjeRestitution2D` | Pace a Purkinje network coupled to a 2D slab over several beats, antegrade and retrograde | `monodomainSolver` with `restitutionEikonalSolver1D` or `monodomain1DSolver`, variant selected at run time | lightweight or full | `Alltest-regression` | per-node network activation time, `Vm` and junction current under `postProcessing/` |
| `NiedererEtAl2011verification` | Run the Niederer slab electrophysiology benchmark | `monodomainSolver` verification workflow | lightweight or full | `Alltest-regression` | activation probes, smoke-check fields, summaries and plots |
| `electromechanicsProtocols/springSupportedSlab` | Show spring (Robin) supports in electromechanics: slab shortening from isometric to free as the end-spring stiffness drops | `electroMechanicalModel` with monodomain electrophysiology and `solidSpringDashpot` end supports | full only | `Alltest-regression` (expected skip in lightweight mode) | end displacement and end force (`solidForces`), `Ta` and `Vm` probes |
| `manufacturedSolutions/monodomainPseudoECG` | Verify monodomain fields and pseudo-ECG against manufactured solutions | `monodomainFDAManufactured` plus pseudo-ECG verification | lightweight or full | `Alltest-regression` | manufactured error summaries and pseudo-ECG series under `postProcessing/` |
| `manufacturedSolutions/monodomainPseudoECG/sphericalShell` | Verify the monodomain on a spherical heart shell with sealed curved walls, a rotated full conductivity tensor and a cubed-sphere or gmsh-tetrahedra mesh | `monodomainFDAManufactured` with `sealedHeartBoundary`, `setManufacturedFields` and `manufacturedShellMonodomainVerifier`; gmsh and `gmshToFoam` for tetrahedra | lightweight or full | `Alltest-regression` | manufactured field-error summaries under `postProcessing/` |
| `manufacturedSolutions/monodomainLeadFieldECG` | Monodomain lead-field functional verification | FDA manufactured voltage | lightweight or full | `Alltest-regression` | lead-field trace and analytical functional errors |
| `manufacturedSolutions/bidomainLeadFieldECG` | Bidomain lead-field functional verification | FDA manufactured voltage | lightweight or full | `Alltest-regression` | lead-field trace and analytical functional errors |
| `manufacturedSolutions/bidomain` | Verify spatial bidomain convergence | `bidomainFDAManufactured` | lightweight or full | `Alltest-regression` | manufactured field-error summaries under `postProcessing/` |
| `manufacturedSolutions/bidomain/sphericalShell` | Verify the heart-only bidomain on a spherical heart shell with sealed curved walls, unequal rotated full conductivity tensors and a cubed-sphere or gmsh-tetrahedra mesh | `monodomainFDAManufactured` with `sealedHeartBoundary`, `setManufacturedFields` and `manufacturedShellBidomainVerifier`; gmsh and `gmshToFoam` for tetrahedra | lightweight or full | `Alltest-regression` | manufactured field-error summaries under `postProcessing/` |
| `manufacturedSolutions/bathBidomain` | Verify the bidomain with a bath: heart and bath fields and the heart-bath interface | `bathBidomainFDAManufactured` | lightweight or full | `Alltest-regression` | manufactured error summaries and interface metrics |
| `manufacturedSolutions/eikonalECG` | Verify activation time and template/quadrature ECG calculations | `eikonalSolver` with manufactured eikonal verification | lightweight or full | `Alltest-regression` | activation-time and ECG reference/error series plus summary CSV files |
| `manufacturedSolutions/eikonalLeadFieldECG` | Verify template-based reciprocal and pseudo ECG together | Existing manufactured eikonal cube, uniform template and prescribed sensitivities | lightweight or full | `Alltest-regression` | independent surface-integral errors, template/source checks and dual traces |
| `manufacturedSolutions/eikonalECG/sphericalShell` | Verify the eikonal solver on a spherical heart shell with sealed curved walls, a rotated full conductivity tensor and a cubed-sphere or gmsh-tetrahedra mesh | `eikonalSolver` with `sealedHeartBoundary`, `setManufacturedFields` and `manufacturedShellEikonalVerifier`; gmsh and `gmshToFoam` for tetrahedra | lightweight or full | `Alltest-regression` | manufactured activation-time error summaries under `postProcessing/` |
| `manufacturedSolutions/monodomain1D3D` | Verify 1D-3D monodomain coupling against manufactured solutions | `monodomainSolver` with 1D graph coupling and manufactured verifier | lightweight or full | `Alltest-regression` | 3-D and graph error summaries and coupling diagnostics under `postProcessing/` and `verification/` |
| `manufacturedSolutions/monodomainTotalLagrangianEM` | Verify coupled monodomain and nonlinear solid mechanics | manufactured total-Lagrangian electromechanics workflow | full only | `Alltest-regression` (expected skip in lightweight mode) | `Vm`, `D`, `lambda`, and `Ta` error/convergence tables and plots |
| `idealizedHeart/electroHeart` | Run whole-ventricle electrophysiology with a Purkinje network and pseudo-ECG | `monodomainSolver`, `eikonalSolver` or the hybrid combination, selected at run time | lightweight or full | `Alltest-regression` (`monodomain` on the human tree; all variants and trees with `CARDIAC_REGRESSION_SCOPE=full`) | activation-time probes and pseudo-ECG series under `postProcessing/` |
| `idealizedHeart/electroMechHeart` | Couple whole-ventricle electrophysiology to the solid on the same anatomy | `electroMechanicalModel` with monodomain electrophysiology | full only | `Alltest-regression` (expected skip in lightweight mode) | `Ta` and `D` probes on the solid, activation-time and `Vm` probes on the electro region |
| `idealizedHeart/pathos/conductionBlock` | Model left or right bundle-branch block structurally | `monodomainSolver` with one Purkinje bundle severed, variant selected at run time | lightweight or full | `Alltest-regression` | junction activation-time probes and pseudo-ECG series |
| `idealizedHeart/pathos/ionicPathology` | Model acute ischemia or Brugada syndrome type-1 | `monodomainSolver` with `ionicConstantOverrides` on `TNNP`, variant selected at run time | lightweight or full | not covered | pseudo-ECG series and monodomain fields |
| `configurationChecks` | Check which `electroProperties` configurations are accepted, what they write, and which are rejected | short runs of every myocardium, conduction, coupler and ECG selection, plus invalid ones | lightweight or full | `Alltest-regression` | per-check runs under `runs/`, compared with `regression/manifests/` |

## Common script pattern

Most runnable cases provide:

- `Allrun`: run the simulation and, where configured, post-process it;
- `Allclean`: remove generated output;
- optional `regression/regressionTest.sh`: compare the run against the case's
  `regression/*.reference` rows. Every script sources the shared library
  `regression/lib.sh` at this level, which owns the options (`--check-only`,
  `--report PATH`), the build-mode gate, the reference grammar and the
  comparator; a case script only declares its reference, its Allrun
  arguments and any case-local check.

Read the case-local `README.md` before running a case. It defines its inputs,
parallel options, expected cost, and output details.

## Cross-case regression entrypoint

The aggregate runner discovers every case with a `regression/regressionTest.sh`,
copies it outside the source tutorial tree without its generated output, and
runs the script there. `--list` prints the discovered cases and `--only PATTERN`
restricts a run. Select the build mode explicitly:

```bash
CARDIAC_REGRESSION_BUILD_MODE=lightweight ./tutorials/Alltest-regression
CARDIAC_REGRESSION_BUILD_MODE=with-solids4foam ./tutorials/Alltest-regression
```

A case that needs cardiacFoam built with solids4foam carries an empty marker
file `regression/needs-solids4foam`; today those are the spring-supported
electromechanics slab (`electromechanicsProtocols/springSupportedSlab`), the
electromechanics manufactured solution
(`manufacturedSolutions/monodomainTotalLagrangianEM`) and the idealized-heart
electromechanics case (`idealizedHeart/electroMechHeart`). In lightweight
mode their regression scripts exit 77 and the runner reports an expected
skip. A case that cannot run on some OpenFOAM versions lists them in
`regression/skips-on-openfoam`; today `manufacturedSolutions/bathBidomain`
lists `v2312`, where its `laplacian(conductivityIntracellular,Vm)` lookup
fails, and `manufacturedSolutions/bidomain/sphericalShell` lists `v2312`,
where its Vm and ionic-state error rows differ from the v2412 and v2512
values. Any other exit-77 skip fails the aggregate run.

`CARDIAC_REGRESSION_SCOPE` selects how much a case checks: `standard` (the
default) is the suite gate, `full` widens cases that offer a larger matrix.
Today only `idealizedHeart/electroHeart` reads it: `standard` runs the
`monodomain` variant on the human Purkinje tree, `full` runs all three
solver variants on both trees, six whole-heart runs. In `with-solids4foam`
mode no skip is expected, so a missing solids4foam build fails these cases
instead of skipping them. In CI only the `with-solids4foam` job that builds
solids4foam from the submodule runs them; every other job runs the
regressions in lightweight mode.

### What a regression runs

A regression checks that a case still runs, through every feature it exists
to exercise, and that its results have not changed. It is not the tutorial:
**every regression finishes in a few minutes on a 4-vCPU CI runner**, about
five on the slowest. When the tutorial takes longer, the case's
`regression/regressionTest.sh` runs a smaller configuration of it:

- a shorter `endTime`, ending soon after the last event the reference checks;
- a coarser mesh;
- only the electrodes or variants the reference checks.

The script sets it with `regression_set` or `regression_edit` from
`regression/lib.sh`, for its own run only: the case's `Allrun` and
dictionaries keep the tutorial's settings, so `./Allrun` still runs the full
tutorial. The reference is generated from the regression configuration, and
every reference checks values that move when the result changes; a probe
that reads `-1` (not yet activated) is checked only alongside probes that
have activated.

One case does not meet the time bound: `idealizedHeart/electroMechHeart`.
Its ionic model is integrated over the 123k-cell heart and its solid solved
every step, and its stimulus starts at 2 ms, so even the 10 ms the regression
runs (the tutorial runs 20 ms) take 9 to 15 minutes on a 14-core workstation.

### Idealized heart/torso ECG application control

[idealizedHeart/electroHeartECG](idealizedHeart/electroHeartECG/README.md)
reuses the supplied heart–lung–torso geometry for monodomain/bidomain source
runs, Wilson lead fields and frozen-source forward recording. Run it explicitly;
it is outside the default lightweight regression sweep. Its contract checker
retains the measured difference between the two recording discretisations.
