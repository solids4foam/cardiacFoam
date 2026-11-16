# bathBidomain/sphericalShell

The bath shell: a bidomain heart shell `0.5 < |x| < 1` coupled every step to
the bath ball `bloodPool` (`|x| < 0.5`) and the outer bath `1 < |x| < 1.5`
(`torsoCoupling bath`), verified against the manufactured solution of
[`bath_shell_bidomain_mms.py`](../../generators/bath_shell_bidomain_mms.py)
by `manufacturedShellBathBidomainVerifier`.

The solution carries what a bath adds to the heart-only bidomain on curved
anisotropic interfaces, convex (`|x| = 1`) and concave (`|x| = 0.5`):

- the intracellular space is sealed on both spheres,
  `n.Gi grad(phi_i) = 0`, so `n.Gi grad(Vm) = -n.Gi grad(phi_e)` is non-zero
  (the virtual-electrode current);
- the extracellular current is continuous through `Ge` on the heart side,
  `n.Ge grad(phi_e) = sigma_b dphi_e/dn`, and `phi_e` is continuous;
- `Gi` and `Ge` have unequal anisotropy ratios and are oblique to both
  spheres; the bath conductivity differs between the two bath zones and
  varies in each;
- the outer sphere is insulated (pure Neumann, `phiE` pinned at
  `phiERefPoint`) or grounded at `c0`, with constant data only.

## Problem

`S = sqrt(1 + t)` and the common mode `T = A_T sin(nu t)` (zero in the
default parameter set, on in `alternate`):

```text
heart   phi_e = f + omega g + T p,   phi_i = h + omega m + T p,
        p = f2 + omega q,            Vm = phi_i - phi_e
bath    phi_e = f + T (f2 + beta_j r_j) + (outer closure, outer bath only)
```

with `f = S fhat` and `h = S hhat` the shells' base fields, `omega` zero on
both spheres, `g` the correction that makes the extracellular current
continuous through `Ge`, `m` and `q` the corrections that seal `phi_i`, and
`beta_j r_j` the bath part of the common mode. The construction and its
self-checks are in the generator's docstring and in the P2 bath record.

The ionic model is `monodomainFDAManufactured` (3D); the verifier adds the
applied current `beta Vm - div(Gi grad phi_i)` to the `Vm` equation and the
source of each zone to the `phiE` equation, with its volume mean removed
when the outer sphere is insulated.

## Stack

- myocardium: `bidomainSolver` on the cellZone `myocardium`,
  `sealedHeartBoundary true`, `conductivitySource field` (Gi and Ge written
  by `setManufacturedFields`), `timeCouplingScheme godunov`, `ddt(Vm)`
  `backward`, `bathPredictorCorrector yes`, `bathHeartPhiETrace zeroGradient`;
- bath: the top-level `torso` block with cellZones `bloodPool` and `bath`,
  `bodyAndOrgansConductivity` (written by `setManufacturedFields`),
  `conormalHarmonic` (with `distanceWeightedHarmonic` the cross term of
  the interface faces uses a cell gradient fitted across the kink of
  `phi_e`, and the pool converges at first order in the maximum norm next
  to the concave interface; see the P2 bath record);
- schemes: `sealedWallLeastSquares` for `grad(Vm)` and `grad(phiE)`,
  `Gauss linear corrected`;
- `phiE`: `PCG`/`DIC` to `1e-15`; PIMPLE with one outer corrector.

## Usage

```bash
./Allrun                                  # hex, N 12, insulated
./Allrun -N 24 -outer grounded
./Allrun parallel -mesh tet -N 16 -seed 2
./Allclean
```

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (6 subdomains, scotch), `cardiacFoam` in parallel, `reconstructPar` | serial |
| `-mesh hex\|tet` | blockMesh cubed sphere or gmsh tetrahedra (`shell.geo.template`) | hex |
| `-N n` | cells per unit length; hex: a multiple of 4; tet: element size `1/n` | 12 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-outer insulated\|grounded` | outer sphere insulated, or grounded at the `c0` of the parameter set (0.1 default, -0.05 alternate) | insulated |
| `-parameterSet default\|alternate` | manufactured parameters; `alternate` turns the common mode on | default |
| `-interface distanceWeightedHarmonic\|conormalHarmonic` | `interfaceConductivityInterpolation` | conormalHarmonic |
| `-trace zeroGradient\|global` | `bathHeartPhiETrace` | zeroGradient |
| `-wallGradient sealed\|stock` | `sealedWallLeastSquares` or `leastSquares` for `grad(Vm)` and `grad(phiE)` | sealed |
| `-predictorCorrector yes\|no` | `bathPredictorCorrector` | yes |
| `-outerCorrectors n` | PIMPLE `nOuterCorrectors` | 1 |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.2 |

## Output

`postProcessing/3D_<mesh>_<N>_cells.dat`: the rows `Vm`, `u1`, `u2`, `phiI`
and `phiE` (`phiE_gauge` when insulated: after removing the volume mean of
the error), `phiE*_pool`, `_heart` and `_bath`, the bulk, boundary,
interface-band and gradient rows, and at each electrode the solved `phiE` of
its cell against the exact value at the cell centre
(`phiE*_electrode_<name>`, and `phiE*_bodySurface` over the eight `surface*`
electrodes at `|x| = 1.47`); the comment lines `# electrode` give the
solved value, the exact centroid value and the exact point value.
