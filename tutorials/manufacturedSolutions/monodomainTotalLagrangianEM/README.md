# manufacturedSolutions/monodomainTotalLagrangianEM tutorial

Full coupled electromechanics manufactured-solution verification case built on `cardiacFoam` + `solids4foam` stack.

## Overview

This case verifies the full coupled electromechanical solve, including interior solid momentum balance:

- `Vm` verified with existing monodomain manufactured verifier
- `Ta` computed from numerical `Vm` and numerical fibre stretch
- solid follows manufactured displacement on every outer boundary via tutorial-local `fixedDisplacement` derivative
- tutorial-local `fvOption` (`manufacturedSolidForce`) adds manufactured body force `B = -Div(P)` so interior `D` converges to manufactured field

Stack uses:

- numerical monodomain manufactured electrophysiology
- numerical manufactured active tension (`ManufacturedElectromechanics`)
- numerical nonlinear solid mechanics (`electroMechanicalLaw`)
- exact manufactured displacement on solid boundaries + manufactured mechanical body force

Exercises:

- electro-to-active-tension path
- active-tension-to-solid path
- region coupling and field mapping
- nonlinear solid solve under prescribed manufactured motion and body force

## Technical Implementation

### Manufactured Solution

Manufactured displacement (reference configuration):

```
D = ( Ax x² y,  Ay y² z,  Az z² x ) * sin(t)
```

with fibre `f0 = (1,0,0)`. Active tension:

```
Ta = Tmax * Vm²/(V₀² + Vm²) * (1 + γ(λ - 1))
```

Body force `B` is the symbolic divergence of the first Piola-Kirchhoff stress (compressible neo-Hookean passive part plus active fibre stress), generated in `src/manufacturedSolidForce/B_expr.H`.

All of `Vm`, `D`, `λ`, and `Ta` are rigorous manufactured convergence targets and should converge at formal scheme order.

### Parameter Sources (Single Source of Truth)

Body force `manufacturedSolidForce` does not hard-code physics parameters; deduces them at runtime from solver dictionaries so they never drift:

- `amplitude`, `Tmax`, `V0`, `gamma`, `TaScale` from `constant/electroMechanicalProperties`
- `E`, `nu` from `constant/solid/mechanicalProperties` (`passiveMechanicalLaw`)

Effective active-stress amplitude is `TaScale * Tmax`, matching the `Ta` field the coupler hands to the solid.

### solids4foam Integration

This case requires two `solids4foam` changes:

- `electroMechanicalLaw` looks up runtime `Ta` field and reads fibre fields `f0` and `f0f`. The case ships `0/solid/f0` and `0/solid/f0f` (both uniform `(1 0 0)`) and `0/solid/D`, which carries the manufactured displacement condition on every patch
- `nonLinGeomTotalLagTotalDispSolid` applies `fvOptions` source in its SNES residual (no-rho `fvOptions()(D)` overload), enabling MMS body force

Implementation note: coupled verifier dictionary is the `verificationModel` entry inside `sequentialElectroMechanicalCoeffs`, the same key every other verifier uses. Electromechanical model stores only that `...Coeffs` sub-dictionary, and looks the verifier up by that exact key: a top-level entry, or one under any other name, is silently ignored.

## Usage

### Manual Execution

```bash
./Allrun
./Allrun parallel
```

Local boundary-condition library is compiled from `src/` before case runs. Compiled library is kept case-local under `platforms/$WM_OPTIONS/lib`, so tutorial does not need write access to global `FOAM_USER_LIBBIN`.

Full solids4foam build only: there is no lightweight variant, and `Allrun` stops if `libelectroMechanicalModels` is missing. `regression/regressionTest.sh` covers the case in `Alltest-regression` on a 20³ mesh. With `CARDIAC_REGRESSION_BUILD_MODE=lightweight` it is an expected skip, so it only runs in a `with-solids4foam` regression run; CI currently runs lightweight only.

### Driver-Managed Convergence Sweeps (Suggested)

```bash
driverFoam run --strict --entry manufacturedMonodomainTotalLagrangianEM
```

Refinement sweep follows monodomain manufactured pattern:

- dimensions: `1D`, `2D`, `3D`
- mesh sizes: `10`, `20`, `40`, `80`
- temporal discretization: `dt ~ h²`

Post-processing reports rigorous manufactured convergence for `Vm`, `D`, `λ`, and `Ta`.
