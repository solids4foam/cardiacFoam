# Electrophysiology on the idealized biventricular heart

Three solver combinations on the same idealized biventricular mesh, same
anatomy, selected at run time.

## Stack

| variant | myocardium (3D) | conduction system (1D) | coupling |
|---|---|---|---|
| `monodomain` (default) | `monodomainSolver` | `monodomain1DSolver` (full cable) | `reactionDiffusionPvjCoupler` |
| `eikonal` | `eikonalSolver` | `eikonalSolver1D` | `eikonalPvjCoupler` |
| `hybrid` | `monodomainSolver` | `restitutionEikonalSolver1D` | `eikonalMonodomainPvjCoupler` |

`monodomain` is the healthy baseline other tutorials reference and the
regression-covered variant (see `regression/`); `eikonal` is the cheapest
combination; `hybrid` is the middle ground (full 3D myocardium, cheaper
1D Purkinje).

## Mesh and shared anatomy fields

The mesh (`constant/polyMesh`) and the anatomy fields `fiber`, `sheet`,
`sheetNormal`, `tm`, `tv`, `apicobasal`, and `Conductivity` live once in
the sibling `../mesh/` directory, shared with `electroMechHeart` and
`pathos/{conductionBlock,ionicPathology}`. `Allrun` copies them in each run.

The mesh in `../mesh/` is already scaled to adult-heart-scale SI
metres (from a raw idealized geometry spanning x:[0,4.0] y:[-2.0,3.0]
z:[-2.0,2.0], uniformly scaled by `0.02`), spanning approximately
x:[0,0.08] y:[-0.04,0.06] z:[-0.04,0.04] metres. `Allrun` copies it in
as-is — no scaling happens here.

## Setup & Fields

The transmural convention in this mesh — epicardium = 0,
endocardium = 1 — is specific to this mesh, distinct from standard UVC.

The shipped fields are `fiber` and `sheet` (the fibre/sheet direction
basis), `sheetNormal` (the third orthonormal direction), `tm`
(transmural position), `tv` (chamber/intraventricular position,
LV/RV), `apicobasal` (longitudinal apex-to-base position), and
`Conductivity` (the per-cell anisotropic conductivity tensor,
`volSymmTensorField`, the same field `electroMechHeart` uses).
`*SolverCoeffs` sets `conductivitySource field;` to read it.

`0/activationTime` (`uniform -1`, `zeroGradient`) is tracked source data,
not derived — only the `eikonal` variant reads it (`MUST_READ`); the
`monodomain` and `hybrid` variants compute `activationTime` as an output
and never read it.

## Tissue heterogeneity

`*SolverCoeffs.ionicHeterogeneity` classifies cells into
`endocardialCells`/`mCells`/`epicardialCells` regions (`mode namedRegions;`)
from a scalar field that must be 0 at the endocardium and 1 at the
epicardium — the opposite of this mesh's `tm` convention.
`system/setExprFieldsDict` derives `t = 1 - tm` (run by `Allrun` via
`setExprFields`, after preloading `tm` with `readFields`), and
`ionicHeterogeneity` reads `field t;` with `regions` tiling `[0,1]` at the
standard thresholds `endocardialCells { range (0 0.3); }` / `mCells
{ range (0.3 0.7); }` / `epicardialCells { range (0.7 1); }` and a smooth
`transitionWidth 0.1;`/`transitionMode blend;` (all three variants — a hard
cutoff destabilizes `eikonal`'s gradient-dependent advection term, and the
smooth transition is the more physiologically realistic choice for
`monodomain`/`hybrid` anyway).

`monodomain`'s and `hybrid`'s `ionicHeterogeneity` additionally compose a `gradientAxes.apicobasal`
axis for an apex-to-base APD gradient: `apicobasal` is used directly
as the distance field (0=apex, 1=base, already this mesh's convention),
and `variables (tauSi)` scales the Bueno-Orovio-Cherry-Fenton model's
primary APD-determining time constant. `scalingMin 0.97;`/`scalingMax
1.04;` were calibrated against single-cell APD90 measurements to target
the ~20ms apex-to-base gradient commonly reported in the literature —
shorter APD at the apex, longer at the base. `eikonal` has no ionic
model to scale a constant in (activation-time-only solve), so this
doesn't apply there. In `hybrid` it shifts the regression probe's
activation by <0.001ms and only softens the T wave (still inverted in
V2-V6 over a 0.5s beat): the ~20ms gradient is small against the ~70ms
activation spread.

## Purkinje conduction network

`constant/purkinjeGraph` (copied in by `Allrun` from `../mesh/constant/`)
is a biventricular Purkinje tree grown on this mesh by
`generatePurkinjeTree`, from UVC fields derived from `tm`/`tv`/`apicobasal`
(`uvc_transmural = 1-tm`, `uvc_intraventricular = 2*tv-1`,
`uvc_longitudinal = apicobasal`). Seeds were deduced from the AHA
segmentation `setCardiacAnatomy` computes, using the nearest LV/RV
endocardial point to the basal-septal AHA segments (`{2,3}` for LV, `21`
for RV). `electroProperties.*`'s `conductionNetworkDomains.purkinjeNetwork`
reads it via `graphFile purkinjeGraph;`, root node 0.

All three variants conduct along the tree at ~3.3 m/s: `hybrid`'s
`restitutionEikonalSolver1D` at its calibrated restitution maximum (3.33
m/s), `eikonal`'s `eikonalSolver1D` via `purkinjeCV 3.33`, and
`monodomain`'s `monodomain1DSolver` via `purkinjeConductivity 0.4` (measured
3.16 m/s at 0.35 and 5.8 m/s at 1.5; the former 10.0 gave ~8.7 m/s).

Two trees are available, chosen by `Allrun`'s `human`/`pig` argument and
both copied in as `constant/purkinjeGraph`, so no dictionary changes:

| tree | source file | terminals | PVJs (LV/RV) | PVJs inside the wall |
|---|---|---|---|---|
| `human` (default) | `../mesh/constant/purkinjeGraph` | all on the endocardium | 1142 (387/755) | 1% |
| `pig` | `../mesh/constant/purkinjeGraphPig` | transmural insertion | 920 (360/560) | 83% |

The pig tree uses the same His, LV and RV seeds and growth parameters; only
the terminals differ. LV terminals are chosen by `terminalSelectionModel
weightedField` from the Garcia-Bustos subendocardial/intramural weights
`setPurkinjeMorphometry` writes (122 subendocardial, 238 intramural), RV
keeps `allLeaves`, and both march into the wall (`terminalModel transmural`,
`gradientFollow` to a random depth in 0.25-0.4 of wall thickness); 3 of 560
RV marches leave the mesh. 46 PVJs on the basal rim sit up to 0.5 mm past
the open base plane.

Each junction's current is spread over the tissue within `pvjRadius` of it,
so the radius is the smallest the mesh allows at every junction of both
trees (`1.65e-3`; the coarsest junction cell, mid-wall on the pig tree, is
1.60 mm across). A junction is physiologically a point contact, and a larger
sphere dilutes the current: at the former `2.5e-3` only 26% of the pig
tree's intramural junctions captured the tissue within 40 ms under
`hybrid`. At `1.65e-3` every junction of both trees captures; at `rPvj 1000`
intramural junctions capture after ~11 ms against ~5 ms at the surface, so
the pig tree activates slightly more slowly than the human tree (99% of the
myocardium by 76 vs 74 ms under `hybrid`), and faster below `rPvj` 300
(`untracked_development/idealizedHeart/pvjCouplingSweep`).
"Inside the wall" counts PVJs deeper than 0.15 of the wall thickness. The
generation record is `cases/idealizedBivEllipsoidPig` in
cardiacCoreStandalone.

## ECG electrodes

`ecgDomains.ECG.electrodePositions` are placed by angle around the LV long
axis in omnidriver's LV frame (`compute_lv_frame`: `L` apex-to-base, `S` LV
centre to RV centre, anterior `A = L x S`, here `-z`): V1..V6 at 35, 65, 100,
135, 170, 205 deg from `S` toward `A`, each at its original apex-base height
and 25mm from the nearest tissue. V1 thus faces the anterior RV free wall
(which spans about +-54 deg) and V6 the LV lateral wall. The earlier
reference-frame transfer from a real heart put V1 at 75 deg, past the RV:
this z-symmetric anatomy has no posterior QRS component, so that V1 read a
positive QRS with no R-wave progression. With this placement the hybrid
variant gives V1 rS, V2 RS, V3-V6 R. The mesh is mirror-symmetric in `z`, so
the choice of `-z` as anterior is the frame's convention, not anatomy.
All three variants share the same positions.

## The `eikonal` variant's numerical stability

`eikonal`'s 3D solve (`eikonalMyocardiumDomain`,
`eikonalAdvectionDiffusionApproach true`) is a deferred-correction
nonlinear scheme, only marginally stable — it needs both of:

- `system/fvSolution`'s `PIMPLE.nOuterCorrectors` at **50**, not higher.
  Counter-intuitively, more outer iterations make it *worse*: the scheme's
  fixed point isn't strictly contractive, so the error compounds rather
  than converging (measured directly: 1000 iterations diverged to
  `~1e12`, 5000 iterations to `~1e96`; 50 iterations is clean).
- The smooth `transitionWidth 0.1`/`transitionMode blend` heterogeneity
  above, not a hard cutoff — a discontinuity there destabilizes the
  gradient-dependent advection velocity the scheme builds from
  `grad(activationTime)`.
- `system/fvSchemes`'s `grad(activationTime)` overridden to `Gauss
  linear` (not the shared `leastSquares` default `monodomain` needs).

## Execution

```bash
./Allrun              # monodomain (default)
./Allrun eikonal
./Allrun hybrid
./Allrun monodomain parallel
./Allrun hybrid pig   # any variant on the pig Purkinje tree
```

The regression runs every variant on both trees (`regression/injection.<variant>.reference`
for human, `injection.<variant>.pig.reference` for pig).

`constant/electroProperties` and `system/controlDict` are symlinks to
`.monodomain`/`.eikonal`/`.hybrid`, swapped by `Allrun` per variant (same
pattern `pathos/conductionBlock` uses for its own variants).
`controlDict` differs per variant but all share `endTime 0.04`: the
regression probe activates at 28.8 ms (`monodomain`) and 30.3 ms (`hybrid`)
on the human tree, 31.9 and 34.0 ms on the pig tree; `eikonal` overrides its
own time control internally regardless of `controlDict`, per
`eikonalMyocardiumDomain::applyModelTimeControls`.
