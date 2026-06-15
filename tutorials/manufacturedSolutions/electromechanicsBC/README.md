# manufacturedSolutions/electromechanicsBC tutorial

This tutorial is a full coupled electromechanics manufactured-solution (MMS)
verification case built on the current `cardiacFoam` + `solids4foam` stack.

It uses:

- numerical monodomain manufactured electrophysiology
- numerical manufactured active tension (`ManufacturedElectromechanics`)
- numerical nonlinear solid mechanics (`electroMechanicalLaw`)
- the exact manufactured displacement imposed on the solid boundaries, plus a
  manufactured mechanical body force in the interior so the nonlinear solid
  equations admit the same analytical solution

## What this case is for

This case verifies the full coupled electromechanical solve, including the
interior solid momentum balance:

- `Vm` is verified with the existing monodomain manufactured verifier
- `Ta` is computed from the numerical `Vm` and numerical fibre stretch
- the solid follows the manufactured displacement on every outer boundary
  through a tutorial-local `fixedDisplacement` derivative
- a tutorial-local `fvOption` (`manufacturedSolidForce`) adds the manufactured
  body force `B = -Div(P)` so the interior `D` converges to the manufactured
  field

It exercises:

- the electro-to-active-tension path
- the active-tension-to-solid path
- the region coupling and field mapping
- the nonlinear solid solve under a prescribed manufactured motion and body
  force

## Manufactured solution

The manufactured displacement (reference configuration) is

```
D = ( Ax x^2 y,  Ay y^2 z,  Az z^2 x ) * sin(t)
```

with fibre `f0 = (1,0,0)`. The active tension follows

```
Ta = Tmax * Vm^2/(V0^2 + Vm^2) * (1 + gamma (lambda - 1))
```

and the body force `B` is the symbolic divergence of the first
Piola-Kirchhoff stress (compressible neo-Hookean passive part plus the
active fibre stress), generated in `src/manufacturedSolidForce/B_expr.H`.

So all of `Vm`, `D`, `lambda` and `Ta` are rigorous manufactured convergence
targets and should converge at the formal scheme order.

### Single source of truth for parameters

The body force `manufacturedSolidForce` does not hard-code any physics
parameters. At run time it deduces them from the same dictionaries the solver
uses, so they can never drift:

- `amplitude`, `Tmax`, `V0`, `gamma`, `TaScale`
  from `constant/electroMechanicalProperties`
- `E`, `nu` from `constant/solid/mechanicalProperties` (`passiveMechanicalLaw`)

Note that the effective active-stress amplitude is `TaScale * Tmax`, matching
the `Ta` field the coupler hands to the solid.

## solids4foam changes required

This case requires two small, related `solids4foam` changes:

- `electroMechanicalLaw` looks up the runtime `Ta` field and derives
  `f0f0 = sqr(f0)` from the `f0` field (so the case only ships `f0`)
- `nonLinGeomTotalLagTotalDispSolid` applies the `fvOptions` source in its SNES
  residual (the no-rho `fvOptions()(D)` overload), enabling the MMS body force

## Execution

```bash
./Allrun
./Allrun parallel
```

The local boundary-condition library is compiled from `src/` before the case is
run. The compiled library is kept case-local under `platforms/$WM_OPTIONS/lib`,
so the tutorial does not need write access to the user's global
`FOAM_USER_LIBBIN`.

## Convergence workflow

This case also supports the same driverFoam-style refinement workflow used by
the other manufactured tutorials:

```bash
python3 -m applications.scripts.driverFoam.openfoam_driver all --entry manufacturedElectromechanicsBC
```

The refinement sweep keeps the monodomain manufactured pattern:

- dimensions: `1D`, `2D`, `3D`
- mesh sizes: `10`, `20`, `40`, `80`
- piecewise `dt ~ h^2`

The post-processing reports rigorous manufactured convergence for `Vm`, `D`,
`lambda`, and `Ta`.

Implementation note: the coupled verifier dictionary belongs inside
`sequentialElectroMechanicalCoeffs`. The electromechanical model stores only
that `...Coeffs` sub-dictionary, so a top-level
`electromechanicalVerificationModel` entry would be ignored.
