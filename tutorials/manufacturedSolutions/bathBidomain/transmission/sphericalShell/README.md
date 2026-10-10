# bathBidomain/transmission/sphericalShell

3D heart-in-bath bidomain manufactured solution on a spherical heart shell
around a bath ball, inside an outer bath box, meshed with gmsh tetrahedra:
curved heart/bath interfaces approximated by planar faces, rotated
anisotropic heart conductivities with all six tensor components, spatially
varying bath conductivities and the exact current on the outer boundary.

Origin: this work.

## Problem

In the heart `H` and the bath `B`:

    H:  chi Cm dVm/dt = div(Gi grad phiI) - chi Cm Iion + I_s,   phiI = Vm + phiE
    H:  div(Gi grad phiI) + div(Ge grad phiE) = S_H
    B:  div(sigmaB grad phiE) = S_B
    interface:  [phiE] = 0,  n.Ge grad phiE (heart) = n.sigmaB grad phiE (bath),
                n.Gi grad phiI = 0
    outer:      n.sigmaB grad phiE = exact (Neumann), or phiE = exact (Dirichlet)

The heart is fully immersed and encloses the bath zone `bloodPool`; its whole
boundary is the heart/bath interface. The case shares the verifier, the
interface construction and the error rows of
[obliqueSlab](../obliqueSlab/README.md), and adds what a planar hexahedral slab
cannot show: tetrahedral cells, interfaces whose normal takes every direction,
interface faces that only approximate the analytic spheres, fibres at every
angle to the interface normal, and a pure-Neumann outer boundary as in the
torso of [electroHeartBath](../../../../idealizedHeart/electroHeartBath/README.md).

### Stack

- myocardium solver: `bidomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov`
- ionic model: `bidomainFDAManufactured` (`dimension "3D"`)
- verifier: `manufacturedTransmissionBathBidomainVerifier`, `geometry sphericalShell`
- bath: `bathPotentialDomain` with zones `bloodPool`, `bath`; the outer patch
  `outer` is a `surfaceCurrentPatches` entry whose gradient the verifier sets
  to the exact current over the owner conductivity, with `phiERefPoint` in
  the outer bath; `-outer dirichlet` makes it a `groundPatches` entry with
  the exact potential
- conductivities: `conductivitySource field`; `0/ConductivityIntracellular`,
  `0/ConductivityExtracellular` and `0/bodyAndOrgansConductivity` written by
  `setManufacturedFields`
- interface: `sealedHeartBoundary true`, `sealedWallTrace zeroGradient`,
  `bathHeartPhiETrace global`, `interfaceConductivityInterpolation
  conormalHarmonic`, `bathPredictorCorrector true` (configuration ATI, as in
  electroHeartBath)

## Construction

### Geometry and mesh

`shell.geo`: the box `[-1.5, 1.5]^3` cut by the spheres `|x| = r1 = 0.5` and
`|x| = r2 = 1`; cellZones `bloodPool` (`|x| < r1`), `myocardium` and `bath`
from the gmsh physical volumes, the six box faces in the patch `outer`.
Delaunay tetrahedra of size `1/N`, a fixed random seed and one thread, so a
gmsh version gives the same mesh for the same `N` and seed; tested with gmsh
4.15.2. Interface vertices lie on the spheres and interface faces inside
them, at most about `h^2/(8 r)` from the sphere. `exposed` is a zero-face
patch added to the generated boundary, so the heart sub-mesh puts its
interface faces there.

The verifier checks that every vertex of a cell lies in the closed radial
range of the cell's zone, that every heart/bath face vertex lies on the
sphere of that bath zone, and that the two bath zones share no face. It logs
the largest distance of an interface face centre from its sphere.

| N | element size | cells |
|---|---|---|
| 5 | 0.2 | 1.7e4 |
| 7 | 0.143 | 4.5e4 |
| 10 | 0.1 | 1.3e5 |
| 14 | 0.071 | 3.4e5 |
| 20 | 0.05 | 1.0e6 |

### Conductivities

`Gi = l f f + t s s + n e_c e_c`, `(l, t, n) = (0.17, 0.019, 0.01)`, and
`Ge` likewise with `(0.62, 0.12, 0.08)`, in the fixed frame
`(e_a, e_b, e_c) = Rz(20 deg) Rx(50 deg)`, with `f = cos theta e_a + sin theta e_b`,
`s = -sin theta e_a + cos theta e_b` and
`theta = -20 deg + 90 deg (|x|^2 - r1^2)/(r2^2 - r1^2)`. A fibre field tangent
to a sphere would vanish somewhere on it, so the frame is fixed and the fibre
turns with the radius. The cell tensors use the heart formula everywhere, so
a heart cell whose centre lies just beyond a sphere carries the tensor its
sources use. `sigmaB = sigma0 (1 + 0.2 sin(pi x + 0.2) cos(pi y)(1 + 0.5 sin(pi z)))`,
`sigma0` 0.6 (`bloodPool`) and 0.25 (`bath`).

### Solution

Every field is a polynomial, trigonometric or exponential function of `x`,
`y`, `z`, or a ratio with a denominator that stays positive in its zone and
beyond it, so each zone formula is analytic and the only non-smooth feature
of the exact solution is the normal derivative jump across the spheres.
Without the common mode (`commonMode.amplitude 0`), bath `phiE = f`; heart
`phiE = f + omega g`, `phiI = h + omega m` with
`omega = (|x|^2 - r1^2)(|x|^2 - r2^2)/((r2^2 - r1^2)(r1 + r2))` and `g`, `m`
blended between the spheres, with the sphere normals written as `-x/r1` and
`x/r2`, so the flux continuity and sealing conditions hold exactly. `f` and
`h` are scaled by `S(t) = sqrt(1 + t)`; ionic states follow from `Vm` with
`G = 1 + 0.25 sin(pi x + 0.2) cos(pi y - 1/3) cos(pi z)`. The exact
expressions are in
`src/verificationModels/bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H`,
generated by `generator/shell_mms.py`.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include;
`Allclean` removes `local/` and the mesh. `Allrun` stops with
`run ./Allclean first` if `log.gmsh` exists, and when `gmsh` or `gmshToFoam`
is missing.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar`, parallel `cardiacFoam`, `reconstructPar` | serial |
| `-N n` | cells per unit length: element size `1/n` | 10 |
| `-seed n` | gmsh random seed | 1 |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.02 |
| `-elliptic` | one step of 1e-9 s from the exact fields (overrides `-dt`, `-endTime`) | off |
| `-config 0\|A\|AI\|ATI\|ATBI` | interface configuration, as in [obliqueSlab](../obliqueSlab/README.md) | ATI |
| `-outer neumann\|dirichlet` | exact current or exact potential on `outer` | neumann |
| `-nOuter n` | PIMPLE `nOuterCorrectors` | 1 |
| `-nNonOrth n` | `nNonOrthogonalCorrectors`, a non-negative integer | 1 |
| `-predictor on\|off` | `bathPredictorCorrector` | on |
| `-grad leastSquares\|GaussLinear` | default gradient scheme | leastSquares |
| `-decomposition manual\|scotch` | `manual` reads `constant/cellDecomposition` (4 ranks); `scotch` uses `-np` ranks | manual |
| `-np n` | ranks for `scotch` | 4 |

The manual decomposition gives rank 0 half of the heart plus part of the
adjacent bath, rank 1 the other half of the heart, rank 2 `bloodPool` and
rank 3 the outer bath.

With Neumann data the `phiE` system is singular: the solver projects its
source onto the compatible range and pins the cell at `phiERefPoint`. The
verifier then removes the volume-weighted mean of the `phiE` error over the
whole mesh before taking the `phiE` and `phiI` norms, and writes the removed
constant as `phiEGaugeOffset`.

`postProcessing/3D_<N>_cells.dat` holds the error summary (`field L1 L2
Linf`), with the rows of [obliqueSlab](../obliqueSlab/README.md) and:

- `interfaceFlux_bathSide`: the per-face interface flux error against the
  bath-side exact flux, where `interfaceFlux` uses the heart side.
- `interfaceFluxJumpScale`: the difference of the two exact fluxes through
  the same face. Normal flux is continuous on the sphere, so it measures the
  displacement and tilt of a planar face from the sphere.
- `interfacePartCurrent_1`, `interfacePartCurrent_2`: the current through
  each octant of each sphere; L1 and L2 are means over the octants, Linf the
  largest. The net current through `|x| = r1` is fixed by conservation (the
  `bloodPool` sources, less their share of the compatibility projection), and
  with Neumann data so is the net current through `|x| = r2`; octant
  currents are not.
- `outerCurrent` covers the Dirichlet patches only: on a Neumann patch the
  current is the one the verifier imposes.

The compatibility imbalance logged by the solver is the cell-centre
quadrature of the manufactured sources against the exact outer current. The
sources are about 10 in size while their zone integrals are about 0.2, so it
is about 0.01 to 0.02 up to `N = 14` and varies with the gmsh seed, then
0.006 at `N = 20`; the projection spreads it over every cell.

## Elliptic-only measurements

```bash
./Allrun parallel -elliptic -nOuter 10 -nNonOrth 10 -decomposition scotch -np 8 -N <n>
```

As on the slab, an elliptic-only measurement is one step from the exact
fields and must converge the deferred corrections within it. Iterating the
explicit non-orthogonal corrections to convergence is unstable on some gmsh
meshes (non-orthogonality about 69 degrees: `-N 14 -seed 3`), with the stock
or the jump-aware gradient; the same mesh runs normally with one correction
per step.

## Expected outcomes

Observed L2 order between `N = 14` and `20` and the least-squares order over
`N = 7, 10, 14, 20` (elliptic, `-nOuter 10 -nNonOrth 10`), or between
`N = 7, 10, 14` (full-time, `-endTime 0.02`, default correctors). Errors at
fixed `N` vary by about 3 % with the gmsh seed.

| config / study | field | expected | mechanism |
|---|---|---|---|
| elliptic, ATI, Neumann | `phiE_heart`, `phiE_bath`, `phiE_interfaceBand`, `phiE_boundaryBand` | 1.95 - 2.05; fit 2.3 - 2.4 | second order; faster pre-asymptotic decay at the coarse levels |
| elliptic, ATI, Neumann | `interfacePartCurrent_1`, `interfacePartCurrent_2` | fit 2.2, 2.3 | integrated interface current is second order |
| elliptic, ATI, Neumann | `interfaceFlux`, `interfaceFlux_bathSide` | about 1 (0.7 - 1.5) | per-face flux of a `leastSquares` gradient on irregular tetrahedra; errors cancel in sums |
| elliptic, any | `interfaceFluxJumpScale` | 2.0 | gmsh sphere faces: displacement and tilt O(h^2) |
| elliptic, ATI | `gradJumpAwareTangentialError`; `gradTangentialError` | about 1; constant | `leastSquares` on irregular cells; stock gradient across the kink (S3) |
| elliptic, ATI | `gradHeartTraceError` | grows slowly, about h^-0.25 | zeroGradient `Vm` wall trace with a non-zero normal derivative, and the normal-projected boundary delta of skewed heart cells (S4) |
| elliptic, ATI, Dirichlet | `phiE_heart`, `phiE_bath`; `phiE_boundaryBand` | 1.9; 1.1 | fixedValue face values attributed to the normal foot of skewed boundary cells; Neumann data has no such error |
| ATI, full-time | `phiE`, `phiI` | 2.5 (pre-asymptotic) | as elliptic |
| ATI, full-time | `Vm` | about 0, also at `-dt 0.001`; Linf grows | O(1) `Vm` error in the heart cells on the interfaces: zeroGradient `Vm` wall trace and S4 |
| ATBI, elliptic and full-time | all | diverges | the conormal `Vm` wall trace closes a loop face value - `leastSquares` gradient - wall trace whose gain is 1.2 - 1.4 per correction on these meshes |
