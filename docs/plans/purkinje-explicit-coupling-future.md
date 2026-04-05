# Future Approach: Explicit Coupling with Shared Current

## Motivation

The operator-splitting approach (first implementation) advances the 1D Purkinje Network (PN)
and the 3D myocardium sequentially within each timestep. While straightforward, this introduces
a splitting error: the coupling current injected into the 3D domain is computed from the 1D state
at the *start* of the timestep, not the end. The Vergara-Quarteroni 1D-3D formulation is more
naturally expressed as a simultaneous RHS contribution, where the coupling current appears
explicitly in the assembled 3D linear system at the same timestep.

## How It Differs from Operator Splitting

In operator splitting:
1. Advance 1D PN ODE from `t` to `t+dt` using `Vm_3D` at PVJs as input.
2. Store result; at a later sub-step inject it as a source term into the 3D solve.

In explicit coupling with shared current:
1. Advance 1D PN ODE from `t` to `t+dt` using `Vm_3D` at PVJs.
2. Compute `I_coupling` at each PVJ immediately.
3. Write `I_coupling` into `couplingCurrent_` (a `volScalarField` owned by `MonoDomainSolver`).
4. `MonoDomainSolver` assembles and solves the 3D PDE with `couplingCurrent_` on the RHS
   in the **same** `solve()` call — not a subsequent step.

The distinction is subtle but matters: the coupling current is part of the assembled matrix
RHS rather than a deferred injection, which eliminates one layer of temporal lag.

## Data Flow at Each Timestep

```
purkinjeModel::evolve(t, dt, Vm_3D, couplingCurrent_)
  ├── sample Vm_3D at PVJ cell indices → Vm_pvj[]
  ├── advance 1D ODE system on each branch
  ├── compute I_coupling at each PVJ node
  └── scatter I_coupling into couplingCurrent_ at PVJ cell indices

MonoDomainSolver::evolveExplicit()
  └── solve( chi*Cm*ddt(Vm) == laplacian(...) - chi*Cm*Iion + externalStimulus + couplingCurrent_ )
```

## Interface Changes vs Operator Splitting

| Aspect | Operator splitting | Explicit shared current |
|---|---|---|
| `purkinjeModel::evolve` signature | `evolve(t, dt, Vm)` — returns nothing or scalar residual | `evolve(t, dt, Vm, couplingCurrent&)` — writes into field |
| `couplingCurrent_` ownership | Not needed | Owned by `MonoDomainSolver`; non-const ref passed to `purkinjeModel` |
| PDE assembly | Separate injection step | `couplingCurrent_` added directly in `solve()` call |

## Trade-offs vs Operator Splitting

**Advantages:**
- More faithful to the Vergara-Quarteroni coupling: coupling current is part of the assembled RHS.
- Opens the door to subcycling: the 1D PN can be advanced with smaller `dt_1D` within one `dt_3D`
  step, writing the time-averaged `I_coupling` into `couplingCurrent_` before the 3D solve.
- No additional data copies; `couplingCurrent_` is reused across timesteps like `externalStimulusCurrent_`.

**Disadvantages:**
- Slightly more intrusive: `MonoDomainSolver` must own and zero `couplingCurrent_` each step.
- The coupling remains explicit (first-order in time) unless an outer Picard iteration is added.
- Subcycling complicates timestep control.

## Future: VTK Polydata Network Input

For networks derived from imaging pipelines (segmentation software, registration tools),
the OpenFOAM-dict edge list becomes unwieldy. A cleaner approach is to load the Purkinje
network topology from a VTK polydata file (`.vtk`, unstructured grid with `LINES` cells),
consistent with how `greensFunctionECGElectro` loads `torsoSurface` via `triSurface`.

**Dict entry:**
```
purkinjeNetwork
{
    purkinjeModel  purkinjeNetworkModel;
    networkFile    "constant/purkinjeNetwork.vtk";   // LINES polydata
    ...
}
```

**Implementation sketch:**
- Use `vtkUnstructuredGrid` or a lightweight VTK reader (analogous to the `triSurface` path)
  to parse node coordinates and line-cell connectivity.
- Map node coordinates to 3D cell centroids for PVJ identification (nearest-neighbour search
  via `meshSearch` or `pointMesh`), replacing the explicit `pvjCellID` list in the dict.
- Branch parameters (conductance, capacitance, ionic model type) can be stored as VTK
  point/cell data arrays, or fall back to uniform values from the dict.

**Why defer this:**
- The dict edge list (first implementation) is sufficient for idealised test geometries and
  benchmarks (e.g., Niederer et al. slab, single-ventricle phantoms).
- VTK input adds a file-format dependency and a coordinate-mapping step that should only
  be introduced once the core solver is validated.

---

## Call Sequence Sketch

```
evolveExplicit():
    couplingCurrent_ = 0
    if (purkinjeModelPtr_):
        purkinjeModelPtr_->evolve(t0, dt, Vm_, couplingCurrent_)
    solve( chi*Cm*ddt(Vm) == laplacian(conductivity, Vm)
                           - chi*Cm*Iion
                           + externalStimulusCurrent_
                           + couplingCurrent_ )
```
