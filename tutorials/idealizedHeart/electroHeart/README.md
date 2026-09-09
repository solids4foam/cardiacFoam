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

These are not strictly UVC coordinates. The transmural convention in
this mesh is epicardium = 0, endocardium = 1.

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
`endocardialCells`/`mCells`/`epicardialCells` bands from a scalar field
that must be 0 at the endocardium and 1 at the epicardium — the
opposite of this mesh's `tm` convention. `system/setExprFieldsDict`
derives `t = 1 - tm` (run by `Allrun` via `setExprFields`, after
preloading `tm` with `readFields`), and `ionicHeterogeneity` reads
`field t;` with the standard thresholds `endoMInterface 0.3;` /
`mEpiInterface 0.7;` and a smooth `transitionWidth 0.1;`/`transitionMode
blend;` (all three variants — a hard cutoff destabilizes `eikonal`'s
gradient-dependent advection term, and the smooth transition is the more
physiologically realistic choice for `monodomain`/`hybrid` anyway).

`monodomain`'s `ionicHeterogeneity` additionally composes an `apexBaseBands`
sub-dict for an apex-to-base APD gradient: `apicobasal` is used directly
as the distance field (0=apex, 1=base, already this mesh's convention),
and `variables (tauSi)` scales the Bueno-Orovio-Cherry-Fenton model's
primary APD-determining time constant. `scalingMin 0.97;`/`scalingMax
1.04;` were calibrated against single-cell APD90 measurements to target
the ~20ms apex-to-base gradient commonly reported in the literature —
shorter APD at the apex, longer at the base. `eikonal` has no ionic
model to scale a constant in (activation-time-only solve), so this
doesn't apply there; not yet ported to `hybrid`.

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

## ECG electrodes

`ecgDomains.ECG.electrodePositions` are transferred to this anatomy from a
reference heart's validated V1-V6 placement, not measured on this mesh — a
normalized approximation, not patient-specific placement. 24-46mm from the
epicardium, comparable to the reference case's own spread. All three
variants share the same positions.

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
```

`constant/electroProperties` and `system/controlDict` are symlinks to
`.monodomain`/`.eikonal`/`.hybrid`, swapped by `Allrun` per variant (same
pattern `pathos/conductionBlock` uses for its own variants).
`controlDict` differs per variant because `hybrid`'s coupling needs ~100ms
of simulated time before activation shows (`monodomain` completes within
20ms; `eikonal` overrides its own time control internally regardless of
`controlDict`, per `eikonalMyocardiumDomain::applyModelTimeControls`).
