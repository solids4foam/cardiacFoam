# eikonal/coupled1D3DSphere

## What it does

Eikonal activation from a one-dimensional graph into a sealed tissue sphere, on
gmsh tetrahedra, against an exact radial solution with no manufactured source.
A root and one terminal at the centre set the time `Tk`; through
`eikonalPvjCoupler` the terminal seeds the ball of cell centres within `R` at
`Tk` plus the seed profile, they stay held until the equation reaches them
earlier, and the front spreads to a sealed wall. The case writes the errors
against the exact solution in the bulk, at the wall and beside the ball.

- **Exercises:** the Purkinje graph to tissue junction of the eikonal model
  (`eikonalSolver1D` on the graph, `eikonalPvjCoupler` with the `flat`,
  `front` and `curvedFront` seed profiles, the held ball), and the eikonal
  operators with the curvature of the front and the exit layer at the sealed
  wall. Only this case's regression checks the eikonal junction against an
  exact solution.
- **Use it for:** verifying a change to the eikonal PVJ coupler or its seed
  profiles; measuring how the activation error depends on `Pe_L`, `R/L` and
  `h/L`; comparing gradient, upstream and convection schemes on a curved
  front against an exact solution (`-gradient`, `-upstream`, `-convection`,
  and `-seedProfile exact`, which holds the ball at the exact field with no
  graph).

## Regression

Tier: `standard`.

`regression/regressionTest.sh` runs `./Allrun` three times in serial, once
per seed profile: `front` and `flat` with the case defaults (`N = 10`,
`Pe_L = 46`, `R/L = 0.25`), and `curvedFront` with `-PeL 21 -N 20`, so that
its critical radius is resolved. The three runs are kept because each takes a
different seeding path through the coupler (`flat` deposits the terminal time;
`front` and `curvedFront` add seed offsets from the metric distance and from
the radial solution) and each runs in serial. Besides the `End` line,
`log.cardiacFoam` of each run must contain `PIMPLE: converged`, so a run that
stops at `nOuterCorrectors` fails. Each compares 5 L2 rows of
`postProcessing/radial_errors.dat` (`eps`, `eps_bulk`, `eps_wall`, `eps_ball`,
`grad_bulk`), 15 rows in all, with `regression/coupled1D3DSphere.reference`
(`front`), `coupled1D3DSphere.flat.reference` and
`coupled1D3DSphere.curvedFront.reference`. The expected values come from gmsh
4.15.2 at mesh seed 1, and the tolerances cover the spread over mesh seeds 1
to 4, so they are wider than a frozen-value check. `-seedProfile exact` and
`-axes` are not run. The marker `regression/needs-tools` (gmsh) applies:
without gmsh the script exits 77, an expected skip. With `--check-only` it
compares only the seed profile of the existing run.

```bash
./regression/regressionTest.sh
```

## Problem

Origin: this work.

With isotropic `M = m I` and a field of the radius alone, the eikonal equation
becomes `-m (psi'' + 2 psi'/r) + V psi' = 1`, `V = c0 sqrt(m)`, with `psi(R) = Tk`
and `psi'(Ro) = 0`. The closed form is in
[`radial_eikonal_mms.py`](../../generators/radial_eikonal_mms.py), and
`manufacturedRadialEikonalVerifier` evaluates it inside the solver with no manufactured
source. The surface condition follows the seed profile: `psi(R) = Tk` for `flat`,
`psi(R) = Tk + R/V` for `front`, and for `curvedFront` the same solution measured from
the critical radius `r_c = 2 sqrt(m)/c0`, where `psi(r_c) = Tk`. A time at a fixed point
depends on `R` through `R/L` for `flat`, through a term `(2/Pe_L) ln(R/L)` for `front`,
and only through the far wall for `curvedFront`. The parameters of
the problem are `Pe_L = c0 L/sqrt(m)`, with `L = Ro - R`, and `R/L`; the mesh adds
`h/L`.

## Running

The case needs gmsh and gmshToFoam. Run `./Allclean` before running it again.

    ./Allrun -N 20 -PeL 46

| Option | Meaning |
|---|---|
| `-N` | cells across the wall: gmsh size `L/N` |
| `-PeL` | `Pe_L` |
| `-ballRatio` | `R/L` |
| `-seedProfile` | seed time of the ball cells. Through one terminal of the graph and `eikonalPvjCoupler`: `front` (`Tk + rho/c0`, default), `flat` (all at `Tk`) or `curvedFront` (the exact unbounded solution from `rho = 2/c0`), with `rho` the metric distance to the terminal. `exact`: no graph, the verifier holds the ball at the exact field, which isolates the operators from the ball |
| `-convection` | `linearUpwind` or `upwind` |
| `-axes "ax ay az"` | speed ratios along the principal axes, turned about (1, 1, 1): the domain and the ball become ellipsoids, which the coupler cannot seed: the verifier holds the ball with `-seedProfile flat` (the default here) or `exact`, with no graph |
| `-gradient` | `grad(activationTime)`: `sealed`, `stock`, `gaussLinear` or `pointCells` |
| `-upstream` | gradient of the upstream face value: `stockLimited`, `sealed`, `stock` or `pointCellsLimited` |
| `-conormalStabilisation` | conormal magnitude as the implicit diffusion coefficient |
| `-relaxation x\|none` | equation relaxation of `activationTime` |
| `-seed` | gmsh random seed |

`postProcessing/radial_errors.dat` holds the errors relative to the activation
span and to the RMS exact gradient: `eps` over the free cells and its parts in
the bulk, the wall cells and the cells beside the ball, and the bulk gradient of
the solved field and of the exact field. Every row uses the whole-domain volume,
so the parts add in quadrature.
