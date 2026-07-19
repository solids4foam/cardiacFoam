# Future electromechanics MMS mechanical closure

Deferred plan. This is future work, not part of the current active
tension implementation. The electrical MMS should remain closed by the
manufactured ionic model; the missing piece is the solid mechanics MMS
closure.

## Current state

The current implementation has three separate pieces:

- The monodomain electrical problem can be closed by a manufactured
  ionic model.
- `ManufacturedElectromechanics` computes the active tension from the
  electrical coupling signal and fibre stretch:

  ```text
  Ta = f(Vm, lambda)
  ```

- `sequentialElectroMechanical` registers a `Ta` field on the solid
  mesh before `solid().evolve()`.

The solids4foam `electroMechanicalLaw` already consumes that field by
looking up `Ta` in the solid mesh object registry and adding the active
stress contribution:

```text
sigma_active = symm(F . (Ta f0f0) . F^T) / J
```

Therefore, the remaining problem is not passing `Ta` to the mechanics.
The remaining problem is making the chosen manufactured displacement
field satisfy the solid momentum equation.

## Intended MMS split

Keep the closures separated:

```text
electrical equation:
    closed by manufactured ionicModel

active tension law:
    closed by ManufacturedElectromechanics

solid mechanics equation:
    closed by future manufactured mechanical source and/or exact
    mechanical constraints
```

The exact electrical solution `Vm_exact` is used to define the reference
solution and error calculation. During the actual run, the active tension
model receives the numerical electrical coupling signal; for a converged
electrical MMS this signal should approach `Vm_exact`.

## Missing mechanical closure

For a selected manufactured displacement field `D_exact`, future work
should derive a mechanical residual/source such that the solid equation
is satisfied by construction.

For the active-stress mechanical law:

```text
F_exact      = I + grad(D_exact)^T
lambda_exact = |F_exact f0|
Ta_exact     = ManufacturedElectromechanics(Vm_exact, lambda_exact)
sigma_exact  = sigma_passive(D_exact) + sigma_active(Ta_exact)
```

The mechanical MMS closure should enforce:

```text
rho d2D_exact/dt2 = div(sigma_exact) + source_mms
```

or the corresponding sign convention used by the selected solids4foam
solid model.

This `source_mms` is the mechanical equivalent of the manufactured
ionic closure used by the monodomain problem.

## Candidate implementation routes

1. Add a manufactured body-force/source term to the solid displacement
   equation.

   This is the cleanest full MMS route. The source should be computed
   from `D_exact`, `Ta_exact`, the passive material law, and the
   selected solid-model equation form.

2. Use exact displacement boundary conditions.

   This is simpler and useful for a constrained test, but it is weaker
   than a full MMS unless the interior residual is also closed.

3. Use exact traction boundary conditions.

   For traction patches, compute:

   ```text
   traction_exact = sigma_exact n
   ```

   This should be combined with a compatible body force if the goal is a
   full convergence MMS.

4. Combine body force and exact boundary conditions.

   This is the preferred final form for a robust electromechanics MMS
   case.

## Transient history requirement

If `D_exact` is time dependent, initializing only the current `D` field
is not enough. The old-time displacement levels must also be consistent
with `D_exact`, otherwise the first computed `d2dt2(D)` term is not the
manufactured acceleration.

For the first implementation, either:

- use a static/quasi-static mechanical MMS, or
- initialize all required old-time displacement fields consistently.

## Verification model extension

Future `electromechanicalVerificationModel` hooks may include:

- initialize exact solid displacement history,
- apply exact mechanical boundary constraints,
- provide a manufactured mechanical source field,
- post-process `D`, `lambda`, `Ta`, and possibly mechanical residuals.

The current initialize/post-process hooks for `Vm`, `D`, and `Ta` can
remain, but they are not sufficient for a fully closed solid mechanics
MMS.

## Tutorial case requirements

A future tutorial case should combine:

- manufactured ionic model for monodomain closure,
- `ManufacturedElectromechanics` active tension model,
- electromechanics verifier,
- consistent `f0` field,
- exact or manufactured solid boundary conditions,
- future mechanical source/force closure,
- final error report for `Vm`, `D`, and `Ta`.

The expected outcome is a coupled electromechanics MMS where the
electrical equation, active tension law, and solid displacement equation
are all closed by manufactured data.
