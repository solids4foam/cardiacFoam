# bathBidomain/transmission/obliqueSlab

2D heart-in-bath bidomain manufactured solution on a slab of myocardium between
two bath slabs, with planar heart/bath interfaces oblique to the mesh axes,
rotated anisotropic heart conductivities, spatially varying bath
conductivities and skewed, non-orthogonal cells at the interfaces.

Origin: this work.

## Problem

In the heart `H` and the bath `B`:

    H:  chi Cm dVm/dt = div(Gi grad phiI) - chi Cm Iion + I_s,   phiI = Vm + phiE
    H:  div(Gi grad phiI) + div(Ge grad phiE) = S_H
    B:  div(sigmaB grad phiE) = S_B
    interface:  [phiE] = 0,  n.Ge grad phiE (heart) = n.sigmaB grad phiE (bath),
                n.Gi grad phiI = 0
    outer:      phiE = exact (Dirichlet)

The heart is fully immersed: its whole boundary is the heart/bath interface.
The exact fields satisfy the interface conditions identically (level-set
transmission construction), so the case verifies the coupled heart/bath
discretisation across a sealed, conductivity-jumping interface rather than a
heart-only or bath-only operator.

### Stack

- myocardium solver: `bidomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov`
- ionic model: `bidomainFDAManufactured` (`dimension "2D"`)
- verifier: `manufacturedTransmissionBathBidomainVerifier`
- bath: `bathPotentialDomain` with zones `bathLow`, `bathHigh`, Dirichlet
  `groundPatches` `outerLow`, `outerHigh` (values set by the verifier)
- conductivities: `conductivitySource field`; `0/ConductivityIntracellular`,
  `0/ConductivityExtracellular` and `0/bodyAndOrgansConductivity` written by
  `setManufacturedFields`
- interface: `sealedHeartBoundary true`, `sealedWallTrace conormal`,
  `bathHeartPhiETrace global`, `interfaceConductivityInterpolation
  conormalHarmonic`, `bathPredictorCorrector true` (configuration ATBI)

## Construction

### Geometry and mesh

With `alpha = 25 deg`, `nHat = (cos alpha, sin alpha, 0)` and
`tHat = (-sin alpha, cos alpha, 0)`, a point is `x = s nHat + tau tHat`. The
domain is `0 < s < 2`, one period `0 < tau < 1` (cyclic `sideA`/`sideB`, a
translation by `tHat`), and one `empty` cell of thickness 0.1 in `z`.

| block | s | cellZone | shear |
|---|---|---|---|
| 0 | 0 - 0.2 | `bathLow` | 0 |
| 1 | 0.2 - 0.5 | `bathLow` | `cB` (default -0.4) |
| 2 | 0.5 - 1.5 | `myocardium` | `cH` (default 0.6) |
| 3 | 1.5 - 1.8 | `bathHigh` | `cB` |
| 4 | 1.8 - 2.0 | `bathHigh` | 0 |

A block with shear `c` has its s-direction cell edges along `nHat + c tHat`,
so the heart and inner-bath cells on each side of an interface are skewed in
opposite senses: the interfaces carry non-orthogonality and skewness, while
the outer bath layers are orthogonal to the Dirichlet patches. The interfaces
`s = 0.5` (outward heart normal `-nHat`) and `s = 1.5` (`+nHat`) stay planar.
`N` cells per unit length in `s` and per period in `tau`; every block width
times `N` is an integer for `N` in {10, 20, 40, 80, 160}. `exposed` is a
zero-face patch named by `exposedFacesPatch`, so the heart sub-mesh puts its
interface faces there; the verifier stops if it is missing or has faces.
`alpha` in `system/blockMeshDict` and `geometryCoeffs.alpha` in
`constant/electroProperties` must agree.

### Conductivities

`Gi = R(theta) diag(0.17, 0.019) R(theta)^T (+) 0.019`,
`Ge = R(theta) diag(0.62, 0.12) R(theta)^T (+) 0.12`, with
`theta(s) = -20 deg + 90 deg clip((s - 0.5)/1, 0, 1)`. In each bath zone
`sigmaB = sigma0 (1 + 0.2 sin(2 pi tau)(1 + 0.5 sin(2 pi s)))`, `sigma0` 0.30
(`bathLow`) and 0.45 (`bathHigh`).

### Solutions

`-solution levelSet` (G1, default): smooth `f`, `h` in `s` and `tau` scaled by
`S(t) = sqrt(1 + t)`; bath `phiE = f`; heart `phiE = f + omega g`,
`phiI = h + omega m` with `omega = (s - 0.5)(s - 1.5)` and `g`, `m` blended
between the interfaces so the flux continuity and sealing conditions hold
exactly. Ionic states follow from `Vm` with
`G = 1 + 0.25 sin(2 pi tau) cos(pi s)`. The exact expressions are in
`src/verificationModels/bathBidomainVerification/manufacturedTransmissionBathBidomainReference.H`.

`-solution affine` (G0, patch test): constant `Gi`, `Ge` (`theta = theta0`)
and `sigmaB` (`sigma0`); `phiE` affine in `s` in each zone with continuous
normal flux, `phiI = S(t) C` constant. All sources vanish except the `Vm`
equation's.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include;
`Allclean` removes `local/`. `Allrun` stops with `run ./Allclean first` if
`log.blockMesh` exists.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar`, parallel `cardiacFoam`, `reconstructPar` | serial |
| `-N n` | cells per unit length in `s` and per period in `tau`; a positive multiple of 10 | 10 |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.2 |
| `-elliptic` | one step of 1e-9 s from the exact fields (overrides `-dt`, `-endTime`) | off |
| `-config 0\|A\|AI\|ATI\|ATBI` | interface configuration (below) | ATBI |
| `-solution affine\|levelSet` | G0 or G1 | levelSet |
| `-shear c` | heart and inner-bath shear both `c`; `-shear 0` is the orthogonal slab; applied before `-shearHeart` and `-shearBath`, whatever the order | 0.6 / -0.4 |
| `-shearHeart c` | heart-block shear | 0.6 |
| `-shearBath c` | inner bath-block shear | -0.4 |
| `-nOuter n` | PIMPLE `nOuterCorrectors` | 1 |
| `-nNonOrth n` | `nNonOrthogonalCorrectors`, a non-negative integer; written to `local/coupling` with `-nOuter` | 1 |
| `-predictor on\|off` | `bathPredictorCorrector` | on |
| `-grad leastSquares\|GaussLinear` | default gradient scheme | leastSquares |
| `-decomposition manual\|simple` | `manual` reads `constant/cellDecomposition`; `simple` is (2 2 1) | manual |

Configurations:

| config | `sealedHeartBoundary` | `sealedWallTrace` | `bathHeartPhiETrace` | `interfaceConductivityInterpolation` |
|---|---|---|---|---|
| 0 | false | zeroGradient | zeroGradient | distanceWeightedHarmonic |
| A | true | zeroGradient | zeroGradient | distanceWeightedHarmonic |
| AI | true | zeroGradient | zeroGradient | conormalHarmonic |
| ATI | true | zeroGradient | global | conormalHarmonic |
| ATBI | true | conormal | global | conormalHarmonic |

The manual decomposition gives rank 0 half of the heart plus part of the
adjacent bath (interface faces internal to a rank), rank 1 the other half of
the heart (interface faces on processor patches) and ranks 2 and 3 one bath
slab each (ranks without heart cells).

`postProcessing/2D_<N>_cells.dat` holds the error summary (`field L1 L2
Linf`): `Vm`, `phiI`, `phiE_heart`, `phiE_bath`, `phiE_interfaceBand` (cells
on either side of an interface face), `phiE_boundaryBand` (cells on
`outerLow`/`outerHigh`), `u1`, `u2`, `interfaceFlux` (per-face
extracellular flux through the interface, area-weighted), `interfaceCurrent_1`,
`interfaceCurrent_2` (net current through each interface),
`gradTangentialError` (tangential part of the base-mesh `grad(phiE)` error on
interface cells), `gradJumpAwareTangentialError` (the same for the gradient
used by the face corrections of the `phiE` operator: jump-aware under
`conormalHarmonic`, `grad(phiE)` otherwise), `gradHeartTraceError` (heart
sub-mesh `grad(Vm)` error on the cells of `exposed`) and `outerCurrent`.

Regression (G1, `-N 20`, short end time, parallel with the manual
decomposition):

```bash
./regression/regressionTest.sh
```

## G0: three roles

```bash
./Allrun -solution affine -elliptic -config ATBI -shear 0 -N 10   # exact
./Allrun -solution affine -elliptic -config ATBI -N 10            # exact
./Allrun -solution affine -elliptic -config ATBI -N 10 -grad GaussLinear   # exact
```

- **Orthogonal slab (`-shear 0`): patch test.** The affine field is
  reproduced exactly by the production path, so `phiE_heart`, `phiE_bath`,
  `phiE_interfaceBand`, `phiI` and `interfaceFlux` are at roundoff
  (Linf at most 1e-8 of the largest exact value) at every `N`, serial and
  parallel. A larger error is a defect in the interface coupling, the flux
  diagnostic or the verifier, not a discretisation error.
- **Sheared slab (default shears): patch test of the non-orthogonal
  corrections.** The stock `leastSquares` gradient of a cell next to an
  interface uses neighbours on both sides of the kink in `phiE`, so it is
  wrong by an amount that does not vanish with `h` (S3): `gradTangentialError`
  stays constant with `N`. The `conormalHarmonic` `phiE` operator builds its
  conormal and non-orthogonal corrections from the jump-aware gradient
  instead, so `phiE_*`, `phiI`, `interfaceFlux` and
  `gradJumpAwareTangentialError` are at roundoff (`phiE` about 5e-13,
  `interfaceFlux` about 2e-12 to 5e-12) at every `N`. The stock gradient
  still seeds the jump-aware gradient and sets the heart trace's interface
  value; neither reaches `phiE`. An O(h) error here means S3 has re-entered
  the `phiE` operator.
- **Sheared slab with `-grad GaussLinear`: S3 off.** The Gauss gradient of
  a field that depends on `s` alone is purely normal on the planar interface
  cells, so the whole sheared-slab error vanishes: `phiE_*`, `interfaceFlux`
  and both gradient rows are at roundoff (`phiE` about 4e-13, `interfaceFlux`
  about 1.5e-12, gradient rows about 1e-11). This is the clean control that
  the sheared-slab error is the `leastSquares` gradient of the interface
  cells.

The three G0 runs are a manual patch test: the regression runs G1 only.

## G1: elliptic-only measurements

```bash
./Allrun -elliptic -nOuter 10 -nNonOrth 10 -N <n>
```

An elliptic-only measurement (`-elliptic`) is one step from the exact fields,
and it must converge the deferred corrections within that step, hence
`-nOuter 10 -nNonOrth 10`. With the defaults (`-nOuter 1`, `-nNonOrth 1`) the
explicit non-orthogonal corrections are evaluated at the exact initial field,
so on skewed meshes the result is not the converged discrete solution: the
errors come out up to about 2 times too small and the observed orders are
biased. On the orthogonal slab (`-shear 0`) the same settings give the
orders in the expected-outcomes table below.

## S4, `gradHeartTraceError` and the controls

The heart sub-mesh sees the interface as its non-coupled `exposed` patch.
OpenFOAM's boundary delta there is normal-projected, `n (n.(Cf - Cn))`, so on
a skewed interface cell the heart `leastSquares` gradient of a field with a
tangential gradient is wrong by an amount independent of `N` (S4). It feeds
the conormal `Vm` wall trace, the heart Laplacian corrections and the `Vm`
equation.

`gradHeartTraceError` measures the heart sub-mesh gradient error at the cells
of `exposed`. It is a mixture of two mechanisms: S3, which reaches the heart
through the global `phiE` trace (the trace's interface value uses the stock
`leastSquares` gradient of the global `phiE`), and, for fields with a
tangential gradient, S4 (the normal-projected boundary delta on skewed heart
interface cells). On sheared G0, which has no tangential gradient, it is pure
S3: 3.0e-3 Linf, independent of `N`. With `-shearHeart 0` it is still
7.8e-4 and flat, so setting the heart shear to zero does not isolate S4.

The shear controls are not a clean S3/S4 separation:

```bash
./Allrun -config ATBI -shearHeart 0 -N <n>        # heart cells orthogonal to the interface
./Allrun -config ATBI -shearBath 0 -N <n>         # inner bath cells orthogonal
./Allrun -config ATBI -grad GaussLinear -N <n>    # leastSquares interface-cell gradient off
```

The clean control is `-grad GaussLinear`: on the sheared G0 it is exact to
roundoff in every row.

## Expected outcomes

Observed order, L2 unless stated, G1 hex. Rows marked elliptic are measured
with `-elliptic -nOuter 10 -nNonOrth 10` and give the order between N = 40
and 80, then between N = 80 and 160; the other rows are full-time runs
(`-endTime 0.2`) and give the order between N = 20 and 40, then between
N = 40 and 80. "Orthogonal" is `-shear 0`, "sheared" the default shears. The
temporal and coupling rows are expectations of the scheme, not orders listed
here.

| config / study | field / norm | expected | mechanism |
|---|---|---|---|
| G0, orthogonal, elliptic, ATBI | `phiE_*`, `phiI`, `interfaceFlux` Linf | at most 1e-8 relative, every `N` | exact patch test |
| G0, sheared, elliptic, ATBI | `phiE_*`, `phiI`, `interfaceFlux`, `gradJumpAwareTangentialError` Linf; `gradTangentialError` | roundoff, every `N`; `gradTangentialError` constant | jump-aware gradient in the operator's corrections; S3 confined to the stock gradient |
| G0, sheared, elliptic, ATBI, `-grad GaussLinear` | `phiE_*`, `phiI`, `interfaceFlux`, both gradient rows Linf | roundoff, every `N` | S3 off: planar-interface Gauss gradient is purely normal |
| 0 | `Vm`, `phiE` | at most 0.3 | unsealed interface: O(1) current leak |
| elliptic, orthogonal, 0 | `phiE` | no convergence | unsealed interface: O(1) current leak |
| A, sheared | `phiE_heart`; `phiE_bath`; `phiE_interfaceBand` Linf | 1.4, then 1.3; 1.1; 1.3, then 1.2 | componentwise interface weights: O(1) flux error in an O(h) layer, on top of S3 |
| AI | `Vm` | 1.0 +- 0.3 | zeroGradient heart `phiE` trace |
| ATI | `Vm` | 0.8 +- 0.3 | zeroGradient `Vm` wall trace |
| ATBI, sheared | `Vm`; `phiE_heart`; `phiE_bath` | 1.9 rising to 2.0; 1.75 rising to 1.9; 1.7 rising to 1.85 | rising towards 2: pre-asymptotic heart-shear term |
| ATBI, orthogonal | `Vm`, `phiE` | 1.9 - 2.0 | no S3 or S4 |
| ATBI, sheared | `interfaceFlux` L2, Linf | 1.5, then 1.8 | rising towards 2 |
| ATBI, sheared | `interfaceCurrent_1`, `interfaceCurrent_2` | 1.7, then 1.9; 1.9, then 2.0 | rising towards 2 |
| ATBI, sheared, `-grad GaussLinear` | `gradTangentialError`; `phiE_heart`, `phiE_bath`; `Vm` (full-time) | `gradTangentialError` O(h); `phiE_heart` 1.7 then 1.9, `phiE_bath` 1.6 (elliptic, 40→80 then 80→160); `Vm` 1.6 then 1.4 (full-time) | S3 off; the Gauss skewness error on the skewed cells keeps `phiE_bath` and `Vm` below the `leastSquares` orders |
| any | `phiE_boundaryBand` | 2.6 - 3.0, about one order above the volume rows | the error vanishes on the Dirichlet patch; orthogonal boundary cells |
| elliptic, orthogonal, A | `phiE` | 1.0 | componentwise interface weights |
| elliptic, orthogonal, ATI | `phiE` | 1.2 | zeroGradient `Vm` wall trace |
| elliptic, orthogonal, AI | `phiE` | 1.3 - 1.6, falling with `N` | zeroGradient heart `phiE` trace |
| elliptic, orthogonal, ATBI | `phiE_*`, `interfaceFlux`, `interfaceCurrent_1`, `interfaceCurrent_2` | 2.0 | no S3 or S4 on an orthogonal slab |
| elliptic, sheared, ATBI | `phiE_heart`; `phiE_bath`; `interfaceFlux` | 1.9, then 2.0; 1.85, then 1.95; 1.8, then 1.9 | second order; the heart-shear pre-asymptotic term decays |
| elliptic, sheared, ATBI | `interfaceCurrent_1`, `interfaceCurrent_2` | 1.85, then 1.95; 1.95, then 2.0 | second order |
| elliptic, sheared, ATBI | `gradJumpAwareTangentialError`; `gradTangentialError`; `gradHeartTraceError` | O(h) or better (1.8, then 1.9); constant; constant | jump-aware operator gradient; stock gradient (S3); heart trace gradient (S3 through the trace, S4) |
| temporal (ATBI, fixed N, dt halving) | `Vm`, `phiE` | 1.0 +- 0.15 | Lie splitting with explicit `Iion` |
| coupling (ATBI, `-nOuter`, `-predictor`) | `phiE(nOuter) - phiE(20)` | decreases geometrically with `nOuter` | outer heart/bath coupling iterations |

Decision rule (elliptic and spatial, ATBI): an observed order of G1 between
N = 40 and 80 below 1.5 with `gradJumpAwareTangentialError` not decreasing
means the gradient in the corrections of the `phiE` operator is
order-limiting (S3 in the operator). `gradTangentialError` is the stock
gradient, which the `conormalHarmonic` operator does not use in its
corrections; it stays constant on the sheared slab. `-grad GaussLinear`
removes S3 from the stock gradient (`gradTangentialError` becomes O(h)), but
Gauss linear carries its own skewness error on the skewed cells, so on G1
`phiE_bath` and `Vm` stay below the `leastSquares` orders and `-grad
GaussLinear` alone cannot confirm S3. An order of at least 1.85 with
`gradJumpAwareTangentialError` O(h) means S3 is not active in the operator;
`gradHeartTraceError` stays constant on the sheared slab without limiting the
`phiE` order. Anything else is inconclusive: add N = 160.

On the sheared slab the configurations separate by their elliptic `phiE`
orders between N = 40 and 80: A follows its own row (its operator has no
jump-aware corrections), AI gives about 1.8 (`phiE_heart`) and 1.6
(`phiE_bath`), ATI about 1.2 and 0.8, ATBI the values above.

## Limitation

For configurations 0 and A the solver blends the heart and bath conductivities
componentwise on interface faces. The positivity of `n.S_f.n` for that face
tensor is not checked by the verifier, because the face tensor is internal to
the solver.
