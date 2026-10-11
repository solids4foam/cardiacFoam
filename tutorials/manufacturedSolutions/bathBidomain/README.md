# bathBidomain

Bidomain manufactured solution with a conductive bath: a myocardium slab
`0 <= x <= 1` between two bath slabs, on Cartesian hexahedra (`-mesh hex`,
default) or conformal gmsh tetrahedra (`-mesh tet`). The case checks the
global heart+bath `phiE` equation, the intracellular no-flux and the current
continuity across the two heart-bath interfaces, the bath boundary conditions
(ground and surface currents), the ionic coupling and the predictor-corrector
bath coupling, together with the interface diagnostics of
`bathBidomainInterfaceMetrics`. The exact potentials depend on `x` only and
the interfaces are planar and normal to `x`; the case does not check curved
or oblique interfaces, rotated conductivities or a spatially varying bath
conductivity.

Origin: the bidomain-with-bath problem of the FDA regulatory science tool
(Pathmanathan & Gray 2014).

## Problem

    chi (Cm dVm/dt + Iion) = div(Gi grad(Vm + phiE))      in the myocardium
    div(Gi grad Vm) + div((Gi + Ge) grad phiE) = 0         in the myocardium
    div(sigmaB grad phiE) = 0                              in the bath
    n . Gi grad(Vm + phiE) = 0                             on x = 0, x = 1
    phiE and the extracellular current continuous          on x = 0, x = 1

on `[-1, 2] x [0, 1]^(d-1)` (1D: `y, z` in `[0, 0.1]`; 2D: `z` in
`[0, 0.05]`), with `S = sqrt(1 + t)`, `F = cos(pi x)` and

    Vm   = S F - alpha x/se                                       0 <= x <= 1
    phiE = -k S F + alpha x/se + C(t)                             0 <= x <= 1
    phiE = -k S + alpha x/sigmaB + C(t)                          -1 <= x <= 0
    phiE = k S + alpha/se + alpha (x - 1)/sigmaB + C(t)           1 <= x <= 2

The ionic states are `u1 = (1 + t) G + S F`, `u2 = G^(-1/2)/(1 + t)` and
`u3 = -alpha x/se`, with `G = 1 + x`, `1 + x y^2` or `1 + x y^2 z^3` in 1D,
2D and 3D; the potentials are the same in every dimension. Coefficients:
`Gi = diag(0.111453302, 0.121585420, 0.030396355)`,
`Ge = diag(0.046165469, 0.050362330, 0.012590582)`, `se = Ge_xx`,
`sigmaB = se/2 = 0.0230827345`, `k = 1/sqrt(2)`, `alpha = 0.01`, `chi = 3`,
`Cm = 2`.

Bath boundary (`-variant`):

- `groundElectrode` (default): `phiE = 0` on `xMin` (`C(t) = k S +
  alpha/sigmaB`) and the current `alpha` into `xMax`.
- `electrodePair`: the currents `-alpha` on `xMin` and `alpha` on `xMax`;
  `phiE` is defined up to a constant, pinned at `phiERefPoint`, and the
  verifier compares `phiE` and `phiI` with their volume means removed.

The `sides` patches carry no current.

### Stack

- myocardium solver: `bidomainSolver`, `solutionAlgorithm implicit`,
  `timeCouplingScheme godunov`, `ddt(Vm)` `backward`
- bath: `bathPotentialDomain` (cellZones `myocardium` and `bath`,
  `bathConductivityField bodyAndOrgansConductivity`), `heartAssembly
  matchedSubmesh`, `interfaceConductivityInterpolation`
  `distanceWeightedHarmonic` or `unweightedHarmonic` (`-interface`)
- bath coupling: `bathPredictorCorrector yes` or `no` (`-predictorCorrector`),
  `nOuterCorrectors 1`, `nNonOrthogonalCorrectors 0`
- ionic model: `bathBidomainFDAManufactured`
- verifier: `manufacturedFDABathBidomainVerifier`
- gradient: `leastSquares` or `Gauss linear` (`-gradient`); laplacian and
  snGrad `corrected`, `limited 0.5` or `orthogonal` (`-snGrad`)
- linear solvers: `PCG`/`DIC`, absolute tolerance `1e-15`, `relTol 0`

## Mesh

- `-mesh hex` (`system/blockMeshDict.<d>D`): three blocks of `N` cells per
  unit length along `x` (and along `y`, `z` in 2D and 3D).
- `-mesh tet` (`three_domain_box.geo.template`, 3D): one gmsh model of the
  three boxes, fragmented so that the interfaces `x = 0` and `x = 1` are
  shared internal faces; Delaunay tetrahedra of size `1/N` from
  `gmsh -setnumber lc 1/N -setnumber seed <seed>`, imported with
  `gmshToFoam`. With gmsh 4.15.2 it gives 14792 tetrahedra at `N = 10` and
  111482 at `N = 20`.

`topoSet` (`system/topoSetDict`) builds the cellZones `myocardium` and `bath`
from the cell centres, and `setTorsoOrganConductivityField` writes the bath
conductivity `0/bodyAndOrgansConductivity`.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

Every option writes a file in `local/`, which the dictionaries include;
`Allclean` removes `local/`, `0/` and the mesh. `Allrun` stops with
`run ./Allclean first` if a mesh log exists.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (6 subdomains, scotch), parallel `cardiacFoam`, `reconstructPar -latestTime` | serial |
| `-mesh hex\|tet` | blockMesh boxes or gmsh tetrahedra (3D) | hex |
| `-dim 1\|2\|3` | dimension of the hex box and of the manufactured states | 3 |
| `-N n` | cells per unit length; tet: element size `1/n` | 20 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-variant groundElectrode\|electrodePair` | bath boundary variant | groundElectrode |
| `-predictorCorrector yes\|no` | `bathPredictorCorrector` | yes |
| `-interface distanceWeightedHarmonic\|unweightedHarmonic` | `interfaceConductivityInterpolation` | distanceWeightedHarmonic |
| `-gradient leastSquares\|GaussLinear` | cell-gradient reconstruction | leastSquares |
| `-snGrad corrected\|limited\|orthogonal` | laplacian and snGrad correction (`limited` is `limited 0.5`) | corrected |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.2 |

Outputs, in `postProcessing/`:

- `<d>D_<n>_cells.dat`: `Vm`, `phiE`, `phiI`, `u1`, `u2`, `u3` errors and the
  bulk, boundary and interface rows of
  [Bulk, boundary and gradient rows](../README.md#bulk-boundary-and-gradient-rows);
  on tetrahedra `n = round(cbrt(myocardium cells))`.
- `bathBidomainInterfaceMetrics.csv`: heart and bath `phiE` errors, the trace
  and current-density jumps, the intracellular leakage and the assembled
  current density at `x = 0` and `x = 1` (`bathBidomainInterfaceMetrics
  -latestTime`); `bathBidomainFaceErrors.csv` holds the signed per-face
  assembled-current error.
- `bathBidomainInterfaceMetricsExactField.csv` and
  `bathBidomainFaceErrorsExactField.csv`: the same with `phiE` replaced by
  the exact field (`-exactFields`), which isolates the interface flux
  construction.

Regression (the defaults: 3D hex, `N = 20`, `groundElectrode`, parallel on
six ranks): the summary and every error row of `3D_20_cells.dat`.

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per question, each run in copies of the case
with the commands in its README; `setup/studies/summarize.py` prints the
tables.

| study | question |
|---|---|
| [cartesianConvergence](setup/studies/cartesianConvergence/README.md) | spatial order on hexahedra, 1D/2D/3D, with and without the predictor-corrector, both bath variants, Cartesian exact-field interface current |
| [tetConvergence](setup/studies/tetConvergence/README.md) | spatial order and interface diagnostics on conformal tetrahedra, exact-field interface current, electrode pair |
| [coupling](setup/studies/coupling/README.md) | effect of the predictor-corrector on the same tetrahedral meshes |
| [interfaceCurrentConvergence](setup/studies/interfaceCurrentConvergence/README.md) | interface diagnostics with `distanceWeightedHarmonic` against `unweightedHarmonic` interface conductivity |
| [gradientScheme](setup/studies/gradientScheme/README.md) | gradient and non-orthogonal correction schemes on tetrahedra |
| [temporalConvergence](setup/studies/temporalConvergence/README.md) | temporal order of the complete bath advance on fixed hexahedral meshes |
| [tetTemporalControl](setup/studies/tetTemporalControl/README.md) | share of the time step in the tetrahedral interface diagnostics |
