# Tutorials

Maintained reference, protocol and verification cases for `cardiacFoam`. The
tables below are the canonical index of the committed runnable cases, grouped by
the equation each case solves. Local research cases and generated sweep
directories are not part of this documented contract.

Every case has its own `README.md`. It opens with the same two sections:
**What it does** (what the case runs, what it exercises, what to use it for)
and **Regression** (its tier, and what the regression runs and checks, or
which case covers it). Inputs, options, cost and outputs follow. The tables
here only say what a case is for and where it is; read the case README before
running it.

## Folder layout

```text
tutorials/
├── electrophysiologyProtocols/    small EP checks: cells, restitution, CV, re-entry, Purkinje
├── electromechanicsProtocols/     small electromechanics checks: supports
├── benchmarks/                    published benchmarks (Niederer et al. 2011 slab, monodomain and bidomain)
├── manufacturedSolutions/         verification against known solutions, by equation
│   ├── generators/                the definition of each manufactured solution
│   ├── monodomain/  eikonal/  bidomain/  bathBidomain/  leadFieldECG/  electromechanics/
├── idealizedHeart/                anatomical biventricular cases (EP, EM, pathologies)
├── configurationChecks/           which electroProperties configurations are accepted
└── regression/                    shared regression library and index check
```

## Run a case

```bash
cd tutorials/manufacturedSolutions/bidomain/cube
./Allrun -h        # options of this case
./Allrun           # mesh, run, post-process
./Allclean         # remove the generated output
```

Most runnable cases provide:

- `Allrun`: run the simulation and, where configured, post-process it;
- `Allclean`: remove generated output;
- optional `regression/regressionTest.sh`: compare the run against the case's
  `regression/*.reference` rows. Every script sources the shared library
  `regression/lib.sh` at this level, which owns the options (`--check-only`,
  `--report PATH`), the build-mode gate, the reference grammar and the
  comparator; a case script only declares its reference, its `Allrun`
  arguments and any case-local check.

`Build mode` below is how `cardiacFoam` has to be built: `lightweight` is
without solids4foam, `full` with it.

## Cases by equation

`Regression` is the case's tier in the cross-case `Alltest-regression` runner
(see [Regression tiers](#regression-tiers)): `standard`, `full`, or `none` with
the cases that cover it. `regression/checkIndex.sh` fails when this column and
the runner disagree.

### Ionic models, no tissue equation

Single-cell and meshless cases: an ionic model integrated without a spatial PDE.

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`electrophysiologyProtocols/singleCell`](electrophysiologyProtocols/singleCell/README.md) | Single-point action-potential runs and ionic-model sweeps | lightweight or full | `standard` |
| [`electrophysiologyProtocols/ionicHeterogeneity`](electrophysiologyProtocols/ionicHeterogeneity/README.md) | Probe transmural ionic heterogeneity without a tissue PDE | lightweight or full | `standard` |
| [`electrophysiologyProtocols/restitutionCurves_s1s2Protocol`](electrophysiologyProtocols/restitutionCurves_s1s2Protocol/README.md) | Generate S1-S2 action-potential-duration restitution curves | lightweight or full | `none` (`singleCellSolver`: `singleCell`; Bueno-Orovio: `ionicHeterogeneity`) |

### Monodomain

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`electrophysiologyProtocols/cableProtocol/monodomain1DCableCV`](electrophysiologyProtocols/cableProtocol/monodomain1DCableCV/README.md) | Measure and refine 1D conduction velocity | lightweight or full | `none` (monodomain: `monodomain/cube`, Niederer) |
| [`electrophysiologyProtocols/rotorInstability`](electrophysiologyProtocols/rotorInstability/README.md) | Exercise sustained re-entry and activation-time behavior | lightweight or full | `full` |
| [`electrophysiologyProtocols/purkinjeRestitution2D`](electrophysiologyProtocols/purkinjeRestitution2D/README.md) | Pace a Purkinje network coupled to a 2D slab over several beats, antegrade and retrograde | lightweight or full | `full` |
| [`benchmarks/NiedererEtAl2011verification`](benchmarks/NiedererEtAl2011verification/README.md) | Run the Niederer et al. (2011) slab benchmark | lightweight or full | `full` |
| [`manufacturedSolutions/monodomain/cube`](manufacturedSolutions/monodomain/cube/README.md) | Verify the monodomain and the pseudo-ECG on the unit cube | lightweight or full | `standard` |
| [`manufacturedSolutions/monodomain/sphericalShell`](manufacturedSolutions/monodomain/sphericalShell/README.md) | Verify the monodomain on a spherical shell with sealed curved walls | lightweight or full | `standard` |
| [`manufacturedSolutions/monodomain/coupled1D3DCube`](manufacturedSolutions/monodomain/coupled1D3DCube/README.md) | Verify the 1D-3D Purkinje-myocardial coupling on the unit cube | lightweight or full | `standard` |
| [`manufacturedSolutions/monodomain/coupled1D3DShell`](manufacturedSolutions/monodomain/coupled1D3DShell/README.md) | Verify the 1D-3D junction against its continuum definition on the shell (needs gmsh) | lightweight or full | `standard` (expected skip without gmsh) |
| [`idealizedHeart/pathos/conductionBlock`](idealizedHeart/pathos/conductionBlock/README.md) | Model left or right bundle-branch block structurally | lightweight or full | `full` |
| [`idealizedHeart/pathos/ionicPathology`](idealizedHeart/pathos/ionicPathology/README.md) | Model acute ischemia or Brugada syndrome type-1 | lightweight or full | `none` (solver: `electroHeart`; overrides: `ionicHeterogeneity`) |

### Eikonal

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`electrophysiologyProtocols/cableProtocol/eikonal1DCableCV`](electrophysiologyProtocols/cableProtocol/eikonal1DCableCV/README.md) | Compare eikonal 1D conduction velocity across resolutions | lightweight or full | `none` (eikonal: `eikonal/cube`) |
| [`manufacturedSolutions/eikonal/cube`](manufacturedSolutions/eikonal/cube/README.md) | Verify activation time and the eikonal ECG on the unit cube | lightweight or full | `standard` |
| [`manufacturedSolutions/eikonal/sphericalShell`](manufacturedSolutions/eikonal/sphericalShell/README.md) | Verify the eikonal solver on a spherical shell with sealed curved walls | lightweight or full | `standard` |
| [`manufacturedSolutions/eikonal/coupled1D3DSphere`](manufacturedSolutions/eikonal/coupled1D3DSphere/README.md) | Verify the path from a graph into tissue against an exact solution (needs gmsh) | lightweight or full | `standard` (expected skip without gmsh) |

### Bidomain

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`benchmarks/NiedererEtAl2011verification/bidomain`](benchmarks/NiedererEtAl2011verification/bidomain/README.md) | Run the Niederer et al. (2011) slab benchmark with the bidomain model, with the published conductivities or an equal-anisotropy set that reproduces the monodomain case | lightweight or full | `full` |
| [`manufacturedSolutions/bidomain/cube`](manufacturedSolutions/bidomain/cube/README.md) | Verify spatial bidomain convergence on the unit cube | lightweight or full | `standard` |
| [`manufacturedSolutions/bidomain/sphericalShell`](manufacturedSolutions/bidomain/sphericalShell/README.md) | Verify the heart-only bidomain on a spherical shell with sealed curved walls | lightweight or full | `standard` (expected skip on v2312) |

### Bidomain with bath

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`manufacturedSolutions/bathBidomain/slab`](manufacturedSolutions/bathBidomain/slab/README.md) | Verify the heart and bath fields and the heart-bath interface | lightweight or full | `standard` |

### Lead-field ECG

The torso conductor and its lead fields, solved with PETSc; the myocardium is
the source.

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`manufacturedSolutions/leadFieldECG`](manufacturedSolutions/leadFieldECG/README.md) | Verify the lead-field ECG against a direct forward solve of the same heart-torso conductor and against the exact potential of an insulated slab and of a layered torso with a current across it | full only, with PETSc | `standard` (expected skip in lightweight mode) |
| [`manufacturedSolutions/leadFieldECG/sphericalShell`](manufacturedSolutions/leadFieldECG/sphericalShell/README.md) | Verify the lead-field ECG on a curved heart shell between a conducting ball and an outer torso bounded by a sphere, with body-surface electrodes | full only, with PETSc | `standard` (expected skip in lightweight mode) |

### Whole heart, solver chosen at run time

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`idealizedHeart/electroHeart`](idealizedHeart/electroHeart/README.md) | Run whole-ventricle electrophysiology with a Purkinje network and pseudo-ECG, with `monodomainSolver`, `eikonalSolver` or the hybrid combination | lightweight or full | `full` (the three solver variants on the human tree) |
| [`idealizedHeart/idealizedTorso`](idealizedHeart/idealizedTorso/README.md) | Run one beat of the bidomain ventricles inside an idealized torso, with the body-surface ECG from lead fields (`leadField`, one-way) or from a heart-and-torso `phiE` (`bath`, two-way) | full; `leadField` also needs PETSc | `none` (couplings: `bathBidomain/slab`, `leadFieldECG`) |

### Electromechanics

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`electromechanicsProtocols/springSupportedSlab`](electromechanicsProtocols/springSupportedSlab/README.md) | Show spring (Robin) supports: slab shortening from isometric to free as the end-spring stiffness drops | full only | `full` (expected skip in lightweight mode) |
| [`manufacturedSolutions/electromechanics/totalLagrangianEM`](manufacturedSolutions/electromechanics/totalLagrangianEM/README.md) | Verify coupled monodomain and nonlinear solid mechanics | full only | `standard` (expected skip in lightweight mode) |
| [`idealizedHeart/electroMechHeart`](idealizedHeart/electroMechHeart/README.md) | Couple whole-ventricle electrophysiology to the solid on the same anatomy | full only | `full` (expected skip in lightweight mode) |

### Configuration

| Case | Purpose | Build mode | Regression |
| --- | --- | --- | --- |
| [`configurationChecks`](configurationChecks/README.md) | Check which `electroProperties` configurations are accepted, what they write, and which are rejected | lightweight or full | `standard` |

## Regression tiers

A case has a regression only if it exercises something that no other
regression exercises: a code path, a coupling, a boundary condition, a
protocol. The regression is one small check of it. A case runs more than one
configuration only when each one exercises something different, or when the
runs are cheap. A case whose code paths are all covered elsewhere (a
pathology that only changes inputs, a measurement protocol on a solver
verified elsewhere) has none, and its README names the cases that cover it.

| Tier | Runs | Cases |
| --- | --- | --- |
| `standard` | by default, and in CI | the manufactured solutions (`manufacturedSolutions/`), `configurationChecks`, and the two meshless ionic cases `singleCell` and `ionicHeterogeneity` |
| `full` | locally, with `CARDIAC_REGRESSION_SCOPE=full` (the standard cases as well) | the published benchmarks, the tissue protocols and the whole-heart cases: physical models on real fronts and anatomy, which the manufactured solutions do not have |
| `none` | never | the cases the index marks `none` |

The standard tier is the gate for every change: the manufactured solutions
check the operators against exact solutions, the meshless cases the ionic
models and their overrides, and `configurationChecks` the configuration
contract.

A case is in the full tier when its regression folder holds the empty marker
file `regression/full-only`.

## Cross-case regression entrypoint

The aggregate runner discovers every case with a `regression/regressionTest.sh`,
copies it outside the source tutorial tree without its generated output, and
runs the script there. `--list` prints the discovered cases with their
markers and `--only PATTERN` restricts a run. Select the build mode
explicitly:

```bash
./tutorials/Alltest-regression --list
CARDIAC_REGRESSION_BUILD_MODE=lightweight ./tutorials/Alltest-regression
CARDIAC_REGRESSION_BUILD_MODE=with-solids4foam ./tutorials/Alltest-regression
CARDIAC_REGRESSION_SCOPE=full CARDIAC_REGRESSION_BUILD_MODE=with-solids4foam ./tutorials/Alltest-regression
CARDIAC_REGRESSION_BUILD_MODE=lightweight ./tutorials/Alltest-regression --only 'bidomain/cube'
./tutorials/regression/checkIndex.sh      # this index against the discovered cases
```

In the standard tier the runner lists the full-tier cases it leaves out; a
full-tier case selected with `--only` runs only with
`CARDIAC_REGRESSION_SCOPE=full`.

A case that needs cardiacFoam built with solids4foam carries an empty marker
file `regression/needs-solids4foam`; today those are
`electromechanicsProtocols/springSupportedSlab`,
`manufacturedSolutions/electromechanics/totalLagrangianEM`,
`idealizedHeart/electroMechHeart` and the two `manufacturedSolutions/leadFieldECG`
cases, whose conductor solve goes through solids4foam's PETSc layer. In
lightweight mode their regression scripts exit 77 and the runner reports an
expected skip. A case that cannot run on some OpenFOAM versions lists them in
`regression/skips-on-openfoam`; today `manufacturedSolutions/bidomain/sphericalShell`
lists `v2312`, where its Vm and ionic-state error rows differ from the
reference in the CI image. A case that needs tools beyond OpenFOAM has the
marker `regression/needs-tools`; today
`manufacturedSolutions/eikonal/coupled1D3DSphere` and
`manufacturedSolutions/monodomain/coupled1D3DShell`. Any other exit-77 skip
fails the aggregate run. In `with-solids4foam` mode no skip is expected, so a
missing solids4foam build fails these cases instead of skipping them. In CI
only the `with-solids4foam` job that builds solids4foam from the submodule
runs them; every other job runs the regressions in lightweight mode.

### What a regression runs

A regression checks that a case still runs, through every feature it exists
to exercise, and that its results have not changed. It is not the tutorial:
**every standard regression finishes in a few minutes on a 4-vCPU CI runner**,
about five on the slowest, and every full regression in a few minutes on a
workstation. When the tutorial takes longer, the case's
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

One full-tier case does not meet the time bound: `idealizedHeart/electroMechHeart`.
Its ionic model is integrated over the 123k-cell heart and its solid solved
every step, and its stimulus starts at 2 ms, so even the 10 ms the regression
runs (the tutorial runs 20 ms) take 9 to 15 minutes on a 14-core workstation.
