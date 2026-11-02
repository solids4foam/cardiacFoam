# monodomainTotalLagrangianEM

Coupled electromechanics manufactured solution on the unit cube: the
monodomain potential, an active tension computed from the numerical `Vm` and
the numerical fibre stretch, and the nonlinear total-Lagrangian solid of
solids4foam driven by a manufactured displacement on every boundary and a
manufactured body force. It checks `Vm`, the displacement `D`, the fibre
stretch `lambda` and the active tension `Ta` against their exact fields, and
with them the electro-to-tension and tension-to-solid paths and the field
mapping between the two regions. It does not check traction or pressure
boundaries (the displacement is prescribed on every face), non-uniform fibres
(`f0 = (1, 0, 0)`), a physiological active-tension model, or a wall treatment
of `Vm` (zero normal gradient on every face). It needs cardiacFoam built with
solids4foam.

Origin: the `Vm` field is the monodomain problem of the FDA regulatory science
tool (Pathmanathan & Gray 2014); the mechanics problem is this work.

## Problem

```text
chi Cm dVm/dt = div(G grad Vm) - chi Cm Iion + I_s     in [0,1]^3,  n . grad Vm = 0
rho d2D/dt2 = Div P + B                                 in [0,1]^3,  D = D_ex on every face
```

in the reference configuration, with

- `Vm = sqrt(1 + t) cos(pi x) cos(2 pi y) cos(3 pi z)`, `G = diag(0.111453302,
  0.121585420, 0.030396355)`, `chi = 3`, `Cm = 2` and the FDA ionic closure
  (`u1`, `u2`, `u3`)
- `D_ex = (Ax X^2 Y, Ay Y^2 Z, Az Z^2 X) sin(t)`, `Ax = Ay = Az = 0.02`
- `Ta = TaScale Tmax Vm^2/(V0^2 + Vm^2) (1 + gamma (lambda - 1))`,
  `lambda = |F f0|`, `f0 = (1, 0, 0)`, `Tmax = 1000`, `V0 = 1`, `gamma = 1`,
  `TaScale = 1`
- `P = J sigma F^-T`, `sigma` the compressible neo-Hookean stress of solids4foam
  `neoHookeanElastic` (`E = 10000`, `nu = 0.3`) plus the active fibre stress
  `F (Ta f0 f0) F^T / J`; `rho = 1060`
- `B = rho d2D_ex/dt2 - Div P(D_ex)`, in `src/manufacturedSolidForce/B_expr.H`

`verification/generate_B_expr_sympy.py` derives `B` with sympy, prints it at
the test points of `verification/verify_B_expr_harness.cpp` (which compiles
`B_expr.H`) and, given `B_expr.H`, compares the two at random points.

### Stack

- coupling: `electroMechanicalModel sequentialElectroMechanical`, regions
  `electro` and `solid` on the same mesh
- electro: `monodomainSolver`, `solutionAlgorithm implicit`, `ddt(Vm)
  backward`, `Gauss linear` gradient, ionic model `monodomainFDAManufactured`,
  verifier `manufacturedFDAMonodomainVerifier`
- active tension: `ManufacturedElectromechanics`
- solid: `nonLinearGeometryTotalLagrangianTotalDisplacement`,
  `electroMechanicalLaw` with the `neoHookeanElastic` passive law, `d2dt2`
  `Euler`; on every patch the case-local condition
  `manufacturedElectromechanicsDisplacement` (`0/solid/D`); the body force by
  the case-local `fvOption` `manufacturedSolidForce` (`system/solid/fvOptions`),
  which reads `amplitude`, `Tmax`, `V0`, `gamma` and `TaScale` from
  `constant/electroMechanicalProperties` and `E`, `nu` and `rho` from
  `constant/solid/mechanicalProperties`
- verifier: `manufacturedElectromechanicsVerifier`, the `verificationModel`
  entry of `sequentialElectroMechanicalCoeffs` (an entry elsewhere is ignored)

The solid relies on two solids4foam features: `electroMechanicalLaw` reads the
runtime `Ta` and the fibre fields `f0` and `f0f` (`0/solid/f0`, `0/solid/f0f`),
and `nonLinearGeometryTotalLagrangianTotalDisplacement` applies `fvOptions` in
its residual.

## Mesh

`system/blockMeshDict.3D`: `N` hexahedra per unit length in each direction,
copied to both regions.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

`Allrun` compiles the case-local library (`src/`) into
`platforms/$WM_OPTIONS/lib` and stops if `libelectroMechanicalModels` is
missing. Every option writes a file in `local/`, which the dictionaries
include; `Allclean` removes `local/`, the meshes and the case-local library.
`Allrun` stops with `run ./Allclean first` if a mesh log exists.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` per region (6 subdomains), parallel `cardiacFoam` | serial |
| `-N n` | cells per unit length | 20 |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.1 |

`postProcessing/manufacturedElectromechanicsSummary.dat` holds the errors of
`Vm`, `D`, `lambda` and `Ta` at the final time.

Regression (`N = 20`, `-dt 0.00224215`): the final time and the L1, L2 and Linf
errors of `Vm`, `D`, `lambda` and `Ta`. With
`CARDIAC_REGRESSION_BUILD_MODE=lightweight` it exits 77 (skip).

```bash
./regression/regressionTest.sh
```

## Notes on dictionary entries

`0/solid/Ta` is the initial active tension. The electromechanical model
registers and solves `Ta`, but the solid model builds its stiffness while it is
constructed, before that field exists, so `electroMechanicalLaw` reads `Ta`
from this file.

`system/solid/fvOptions` deduces all physics parameters at run time from the
authoritative dictionaries, so they cannot drift from the solver:
`amplitude`, `Tmax`, `V0`, `gamma` and `TaScale` from
`constant/electroMechanicalProperties`, and `E`, `nu` from
`constant/solid/mechanicalProperties` (`passiveMechanicalLaw`).
