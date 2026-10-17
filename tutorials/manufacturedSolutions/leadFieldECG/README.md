# leadFieldECG

Verification of the lead-field ECG (`ecgSolver leadFieldECG` with the lead
fields of `buildLeadFields`). The case checks two things:

- that the lead-field channel equals the converged forward
  passive-conductor solve of the same conductor; the lead field is the
  discrete adjoint `A^-T e` of the finite-volume operator `A`, so the two
  agree to the fixed-point tolerance on hexahedra and tetrahedra; and
- that both converge to an exact torso potential.

The myocardium is a monodomain or a bidomain solver on the cellZone
`myocardium`. The torso exists only as the conductor of the ECG domain; no
bath forward solve is run.

The verifier `manufacturedLeadFieldECGVerifier` evaluates each sample in four
ways:

- **leadField**: the channel `z_e^T b(Vm_h) + c_e`.
- **direct**: a direct forward solve of the same conductor on the same
  `Vm_h`. It uses the same reference and boundary values, and repeats the
  assembly and solve until the relative change between passes is below
  `1e-10`.
- **leadFieldStatic** and **directStatic**: the same two, applied to the
  exact `Vm` at the heart cell centres.

It compares these with the exact potential at the electrode point and at the
centroid of the electrode's cell.

## Problem

### `-problem slab`: insulated heart slab

The heart is a slab `0 <= x <= 1` between torso slabs `-1 <= x <= 0` and
`1 <= x <= 2`, on `[-1, 2] x [0, 1]^2`. Every outer face is insulated. With
`S = sqrt(1 + t)` the fields are:

```text
Vm  = S cos(pi x)                                0 <= x <= 1
phi = -k S cos(pi x) + C                         0 <= x <= 1
phi = -k S + C       (x <= 0),     k S + C       (x >= 1)
k   = Gi_xx/(Gi_xx + Ge_xx)
```

where `phi` solves the forward problem:

- `div((Gi + Ge) grad phi) = -div(Gi grad Vm)` in the heart;
- `div(sigma_b grad phi) = 0` in the torso;
- `phi` and the current are continuous at `x = 0` and `x = 1`;
- the outer boundary carries no current.

These are the fields of the FDA bath slab
([bathBidomain](../bathBidomain/README.md)) with `alpha = 0`. The torso
potential is constant in each slab, so the torso conductivity
`sigma_b = 0.0230827345` does not enter.

`C` is fixed by the pin `phiERefPoint (-0.9013 0.0437 0.0437)` with value 0.
The verifier gauges the exact potential to 0 at the centroid of the pinned
cell. The exact ECG between the two slabs is `2 k S`.

Electrodes:

| name | position | exact potential |
|---|---|---|
| `torsoLeft`, `torsoLeftFar` | `(-0.4567 0.4137 0.5871)`, `(-0.8713 0.2311 0.7129)` | 0 (a null lead) |
| `torsoRight` | `(1.5433 0.6173 0.3311)` | `2 k S` |
| `nearLeft`, `nearRight` | `(-0.0213 0.4521 0.5317)`, `(1.0187 0.3719 0.4613)` | 0, `2 k S` |

`nearLeft` and `nearRight` lie in torso cells next to the heart. No
electrode lies on a hexahedral face at `N = 10, 20, 40, 80`. The exact
potential is constant in each torso slab, so the point value and the
centroid value are equal.

Intracardiac electrodes are not supported: `buildLeadFields` and
`leadFieldECG` stop with a fatal error when an electrode's cell is a
myocardium cell.

Conductivities, as `diag(xx, yy, zz)` in S/m:

| myocardium | conductivity | ECG domain Gi | ECG domain Ge |
|---|---|---|---|
| `monodomain` | `G_m = diag(0.111453302, 0.121585420, 0.030396355)` | `(2 + sqrt 2) G_m` | `sqrt 2 G_m` |
| `bidomain` | `Gi = diag(0.111453302, 0.121585420, 0.030396355)`, `Ge = diag(0.046165469, 0.050362330, 0.012590582)` | the solver's | the solver's |

In both rows `k = 1/sqrt 2`, and the monodomain row has
`G_m = Gi (Gi+Ge)^-1 Ge`. `-conductivity isotropic` uses the `xx` value in
every direction. `-conductivity unequal` keeps Gi and every `xx` value and
sets `Ge_yy = Gi_yy`, `Ge_zz = 0.1 Gi_zz` (`G_m` follows), so the
anisotropy ratios of Gi and Ge differ.

The fields depend on `x` only and the tensors have no `xy` or `xz`
components, so only the `xx` values enter the exact solution, and it is the
same for all three choices. What the static rows measure depends on the
choice:

- `anisotropic` and `isotropic`: `Ge` is proportional to `Gi`. The discrete
  operator is linear in the conductivity, so `phi = -k Vm + c` at the heart
  cell centres satisfies the heart rows exactly, on any mesh, apart from
  their coupling across the interface. The static error is that of the
  interface and the torso; on the hexahedra it is one constant offset per
  torso slab.
- `unequal`: on tetrahedra the `yy` and `zz` parts of `Gi - k (Gi+Ge)` act
  on the `x`-only `Vm` through the cell stencils, so the heart interior
  carries a truncation error and the passive heart operator is exercised.
  On the Cartesian hexahedra those parts give zero and the result equals
  that of `anisotropic`.

The myocardium runs the 1D FDA manufactured solution on the 3D heart mesh:

- ionic models `monodomainFDAManufactured` and `bidomainFDAManufactured`, with
  `dimension "1D"`;
- verifiers `manufacturedFDAMonodomainVerifier` and
  `manufacturedFDABidomainVerifier`;
- `chi = 3`, `Cm = 2`.

So `Vm_h` approximates the `Vm` above.

### `-problem layered`: torso layers and a current across the slab

The slab, with each torso slab split into two layers at `x = -0.5` and
`x = 1.5`: the cellZone `torsoInner` (`-0.5 < x < 0` and `1 < x < 1.5`,
conductivity `sigma_1 = 0.1`) next to the heart and `torsoOuter`
(`x < -0.5` and `x > 1.5`, `sigma_2 = 0.02`) beyond it. The torso faces
between the two zones take the distance-weighted harmonic mean of the two
conductivities. A current enters through the outer faces as the surface
currents of the `torso` sub-dictionary:

```text
surfaceCurrentPatches { xMin -alpha; xMax alpha; }     (-alpha, default 0.01)
```

that is `n . sigma grad phi = -alpha` on `xMin` and `alpha` on `xMax`, the
forward's fixed gradient `I/(n . sigma . n)`. The current density is
`-alpha` along `x` everywhere, and the potential is the slab's plus a
piecewise linear part:

```text
phi = phi_slab + alpha R(x),     R(x) = integral from 0 to x of dx'/sigma(x')
sigma = (Gi + Ge)_xx in the heart, sigma_1 or sigma_2 in the torso layers
```

`Vm` is the slab's. `R` is continuous and `sigma R' = 1` in every layer, so
the potential and the current are continuous at `x = -0.5, 0, 1, 1.5`; the
sides stay insulated because the potential depends on `x` only and the
conductivities have no `xy` or `xz` components. `-alpha 0` writes no
surface currents and leaves the outer boundary insulated: the torso
potential is then constant in each torso slab, and the exact ECG is the
slab's `2 k S` for any `sigma_1`, `sigma_2`.

The verifier (`reference layered`, `layerInterfaces (-0.5 1.5)`) takes
`alpha` from the surface currents and the layer conductivities from the
torso cells. It stops with a fatal error when a torso cell crosses a layer
interface, when a layer has more than one conductivity, when the surface
currents are not along `x` or carry a net current, or when a current is
combined with the bidomain myocardium. With a current the heart slope
`alpha/(Gi + Ge)_xx` at the walls carries the intracellular flux
`alpha Gi_xx/(Gi + Ge)_xx`, which the sealed monodomain `Vm` admits and
the insulated bidomain intracellular potential does not. The pin and the
electrodes are those of the slab: `torsoLeftFar` lies in the pinned outer
layer, `torsoLeft` and `nearLeft` in the inner left layer, `nearRight` in
the inner right layer and `torsoRight` in the outer right layer. With
`alpha != 0` the exact potential varies in the torso, so the point and the
centroid values differ.

### Stack

- myocardium:
  - `monodomainSolver` or `bidomainSolver` on the cellZone `myocardium`;
  - `solutionAlgorithm implicit`, `timeCouplingScheme godunov`;
  - `ddt(Vm)` `backward`;
  - the bidomain heart `phiE` is pinned at `(0.1137 0.0437 0.0437)`.
- ECG domain:
  - `leadFieldECG` with a `torso` sub-dictionary: cellZone `bath`
    (`layered`: `torsoInner` and `torsoOuter`, and the surface currents),
    `bodyAndOrgansConductivity`, `distanceWeightedHarmonic`, the pin above;
  - the verifier `manufacturedLeadFieldECGVerifier` with `reference slab`
    or `layered`.
- schemes:
  - gradient `leastSquares`;
  - laplacian and snGrad `corrected` or `orthogonal` (`-snGrad`); the
    lead-field build refuses nonlinear schemes such as `limited`.
- linear solvers, absolute tolerance `1e-15`:
  - `phiE|phiEFinal`: `PCG`/`DIC`, `relTol 0`, for the direct solves and the
    bidomain heart potential;
  - `leadField`: `PCG` with a `GAMG` (`DICGaussSeidel`) preconditioner,
    `relTol 0.01`, and `leadFieldFinal`: the same with `relTol 0`, for the
    lead fields.
- `nNonOrthogonalCorrectors 0`. The lead fields and the direct solves iterate
  their explicit corrections to convergence themselves.

## Mesh

The meshes are those of [bathBidomain](../bathBidomain/README.md).

- `-mesh hex` (`system/blockMeshDict`): three blocks of `N^3` cells.
- `-mesh tet` (`three_domain_box.geo.template`): gmsh Delaunay tetrahedra of
  size `1/N`, with the interfaces `x = 0` and `x = 1` as shared internal
  faces.

`topoSet` builds the cellZones `myocardium` and `bath`.
`setTorsoOrganConductivityField` writes `0/bodyAndOrgansConductivity`.

`-problem layered` uses the `.layered` dictionaries and
`layered_box.geo.template`:

- `-mesh hex` (`system/blockMeshDict.layered`): five blocks, `N/2` cells
  along `x` in each torso layer and `N` in the heart, so every level has
  faces at `x = -0.5` and `x = 1.5` (`N` even);
- `-mesh tet`: the same gmsh tetrahedra with the planes `x = xLeft` and
  `x = xRight` (`DefineConstant`, `-0.5` and `1.5`) as further shared
  internal faces;
- `topoSet` (`system/topoSetDict.layered`) builds `myocardium`, `torsoInner`
  and `torsoOuter`, and `setTorsoOrganConductivityField`
  (`system/setTorsoOrganConductivityFieldDict.layered`) sets their
  conductivities.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

`Allrun` meshes the case, then runs `buildLeadFields` and `cardiacFoam`. Every
option writes a file in `local/`, which the dictionaries include. `Allclean`
removes `local/`, `0/`, the mesh and the lead fields.

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (2 subdomains, scotch); `buildLeadFields` and `cardiacFoam` in parallel | serial |
| `-problem slab\|layered` | exact passive solution | slab |
| `-alpha x` | surface current of `layered` (`-alpha` on `xMin`, `alpha` on `xMax`; 0: insulated); a current needs the monodomain myocardium | 0.01 |
| `-mesh hex\|tet` | blockMesh boxes or gmsh tetrahedra | hex |
| `-N n` | cells per unit length; tet: element size `1/n` | 10 |
| `-seed n` | gmsh random seed (`-mesh tet`) | 1 |
| `-myocardium monodomain\|bidomain` | myocardium solver | monodomain |
| `-conductivity anisotropic\|isotropic\|unequal` | heart conductivities | anisotropic |
| `-snGrad corrected\|orthogonal` | laplacian and snGrad correction | corrected |
| `-dt x` | time step | 0.01 |
| `-endTime x` | end time | 0.05 |

Outputs in `postProcessing/`:

- `leadFieldECG.dat`: the channels.
- `manufacturedLeadFieldECG_ECG.dat`: for every sample and electrode,
  `leadField`, `direct`, `exactPoint`, `exactCentroid`, `leadFieldStatic` and
  `directStatic`, at full precision.
- `manufacturedLeadFieldECGSummary_ECG.dat`: one row `name L1 L2 Linf` (over
  the samples) per error and electrode, and the maximum over the electrodes
  (no suffix). The rows are:
  - `leadField-direct`, and `leadField-direct_relative`, which is that
    maximum divided by the largest `|direct|`;
  - `leadFieldStatic-directStatic`;
  - `leadFieldStatic-exactPoint`, `leadFieldStatic-exactCentroid`;
  - `directStatic-exactPoint`, `directStatic-exactCentroid`;
  - `localisation` (`exactPoint - exactCentroid`).
- `1D_<cells>_cells.dat`: the myocardium verifier's `Vm` errors.

The regression runs the defaults (hex, `N = 10`, monodomain) and
`-problem layered` on the same mesh, each in parallel on two ranks. For each
it checks:

- the reciprocity gap `leadField-direct_relative`, which must be within
  `1e-9` of 0;
- the static rows against the exact potential
  (`regression/leadFieldECG.reference`,
  `regression/leadFieldECGLayered.reference`). In the layered run the
  current crosses the torso layer interface (left electrodes, exact to
  round-off) and both heart walls (right electrodes, the slab's error).

```bash
./regression/regressionTest.sh
```

## Studies

`setup/studies/` holds one study per question. Each study is run in copies of
the case, with the commands in its README. `setup/studies/summarize.py` prints
the tables.

| study | question |
|---|---|
| [leadFieldConvergence](setup/studies/leadFieldConvergence/README.md) | order of the lead field and of the direct solve against the exact slab potential at the torso electrodes, on hexahedra and tetrahedra, with Ge proportional to Gi and with unequal anisotropy ratios |
| [leadFieldReciprocity](setup/studies/leadFieldReciprocity/README.md) | gap between the lead field and the direct forward solve, on hexahedra and tetrahedra, monodomain and bidomain, serial and parallel |
| [leadFieldOrganJump](setup/studies/leadFieldOrganJump/README.md) | order of the lead field and of the direct solve with torso layers of different conductivity and a current across them, and the independence of the insulated ECG from the layer conductivities |
