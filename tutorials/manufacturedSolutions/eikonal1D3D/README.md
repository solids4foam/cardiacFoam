# eikonal1D3D

Eikonal activation from a one-dimensional graph into a sealed tissue sphere, on
gmsh tetrahedra, against an exact solution with no manufactured source. A root
and one terminal at the centre set the time `Tk`; the terminal seeds the ball of
cell centres within `R` at `Tk` plus the seed profile, and they stay held until the equation reaches them
earlier, and the front spreads to a sealed wall. The case checks
the eikonal operators with the curvature of the front and the exit layer at the
wall, the path from the graph into the tissue, and the held ball.

Origin: this work.

## Problem

With isotropic `M = m I` and a field of the radius alone, the eikonal equation
becomes `-m (psi'' + 2 psi'/r) + V psi' = 1`, `V = c0 sqrt(m)`, with `psi(R) = Tk`
and `psi'(Ro) = 0`. The closed form is in
[`radial_eikonal_mms.py`](../generators/radial_eikonal_mms.py), and
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
| `-RoverL` | `R/L` |
| `-ball constant` | one terminal through the graph and `eikonalPvjCoupler`, the ball seeded by `-seedProfile` |
| `-ball exact` | no graph, the ball at the exact field, which isolates the operators from the ball |
| `-seedProfile` | seed time of the ball cells through the coupler: `front` (`Tk + rho/c0`, default), `flat` (all at `Tk`) or `curvedFront` (the exact unbounded solution from `rho = 2/c0`), with `rho` the metric distance to the terminal; needs `-ball constant` and an isotropic tensor |
| `-convection` | `linearUpwind` or `upwind` |
| `-axes "ax ay az"` | speed ratios along the principal axes, turned about (1, 1, 1): the domain and the ball become ellipsoids, and the verifier holds a constant ball with no graph |
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

`regression/regressionTest.sh` runs the defaults for each seed profile. The marker `regression/needs-tools`
makes it an expected skip where gmsh is missing.
