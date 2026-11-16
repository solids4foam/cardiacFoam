# Electromechanics manufactured solutions

Verification of the coupled electromechanical model against manufactured
solutions. Every case here needs cardiacFoam built with solids4foam.

| Case | Geometry | What it verifies |
| --- | --- | --- |
| [`totalLagrangianEM`](totalLagrangianEM/README.md) | unit cube | the monodomain potential, an active tension computed from the numerical `Vm` and fibre stretch, and the nonlinear total-Lagrangian solid of solids4foam, against exact `Vm`, `D`, `lambda` and `Ta` |

## Run a case

```bash
cd tutorials/manufacturedSolutions/electromechanics/totalLagrangianEM
./Allrun -h        # options of this case
./Allrun           # build, run, post-process
```

The regression, from the repository root:

```bash
CARDIAC_REGRESSION_BUILD_MODE=with-solids4foam \
    ./tutorials/Alltest-regression --only 'electromechanics/totalLagrangianEM'
```

In `lightweight` mode the regression exits as an expected skip.
