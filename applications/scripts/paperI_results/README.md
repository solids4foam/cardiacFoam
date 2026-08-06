# Verification experiment interface

`applications/scripts/driverFoam/verification_experiments.json` is the single
public map from a verification experiment to the tutorial code that runs it
and the normalized result it produces. The directory name of these aggregation
helpers is historical; experiment identifiers and paths do not depend on a
manuscript layout.

## Naming

Experiment identifiers are lowercase snake case:

- continuum convergence: `<physics>_<mesh_family>`, for example
  `monodomain_cartesian`;
- tetrahedral variants: `<physics>_tet_<family>`, where the family is
  `generic`, `frontal`, or `conformal`;
- coupled operators: `<coupled_physics>_coupled`;
- isolated operator diagnostics: `<physics>_<operator>_<mesh>`, for example
  `eikonal_gradient_tet`.

Every normalized scalar result is
`<tutorial>/setup/results/<experiment_id>.csv`. A numerical reference, when
the experiment has one, is `<tutorial>/reference/<experiment_id>.csv`.
The registered runner is the stable public entry; older runner names remain
only as compatibility implementations.

## Commands

From the repository root:

```bash
applications/scripts/driverFoam/bin/driverFoam experiment-plan
applications/scripts/driverFoam/bin/driverFoam experiment-plan --experiment eikonal_tet_generic
./reproduce_verification.sh --dry-run
./reproduce_verification.sh --skip-run eikonal_tet_generic
./reproduce_verification.sh monodomain_tet_generic
```

The DriverFOAM plan prints the JSON contract and readiness checks without
running OpenFOAM. `experiment-run --experiment <id>` is an equivalent
DriverFOAM entry to the general reproduction script.

`--skip-run` compares existing normalized outputs without requiring archived
raw meshes or solver directories. A full invocation runs the registered test,
normalizes its output, and compares it with the reference. Experiments with
reference `-` are validated for a non-empty result but do not claim a
numerical regression tolerance.

The registry, runner, result, and reference names—not filenames copied into a
paper repository—define the experiment interface.
