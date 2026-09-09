# Electrophysiology on the idealized biventricular heart

This is a tutorial case for the idealized model (biventricular
ellipsoid), using an existing mesh and supplied anatomy fields.

## Stack

- electro model: `monodomainSolver`
- ionic model: `BuenoOrovio` (or typical PATHOS defaults)
- conduction system: configured via `constant/electroProperties`'s
  `conductionNetworkDomains` block, read directly by `cardiacFoam`

## Mesh and shared anatomy fields

The mesh (`constant/polyMesh`) and the anatomy fields `fiber`, `sheet`,
`sheetNormal`, `tm`, `tv`, `apicobasal`, and `Conductivity` live once in
the sibling `../mesh/` directory, shared with the `electromechanicalHeart`
case. `Allrun` copies them in each run.

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
`volSymmTensorField`, the same field `electromechanicalHeart` uses).
`monodomainSolverCoeffs` sets `conductivitySource field;` to read it.

## Tissue heterogeneity

`monodomainSolverCoeffs.ionicHeterogeneity` classifies cells into
`endocardialCells`/`mCells`/`epicardialCells` bands from a scalar field
that must be 0 at the endocardium and 1 at the epicardium — the
opposite of this mesh's `tm` convention. `system/setExprFieldsDict`
derives `t = 1 - tm` (run by `Allrun` via `setExprFields`, after
preloading `tm` with `readFields`), and `ionicHeterogeneity` reads
`field t;` with the standard thresholds `endoMInterface 0.3;` /
`mEpiInterface 0.7;`.

An `apexBaseBands` sub-dict composes an apex-to-base APD gradient on top
of the transmural bands above (previously unexercised by any tutorial):
`apicobasal` is used directly as the distance field (0=apex, 1=base,
already this mesh's convention — no inversion needed, unlike `tm`), and
`variables (tauSi)` scales the Bueno-Orovio-Cherry-Fenton model's primary
APD-determining time constant. `scalingMin 0.97;`/`scalingMax 1.04;` were
calibrated against single-cell APD90 measurements (not the library's
generic wide-range defaults) to target the ~20ms apex-to-base gradient
commonly reported in the literature and used in simulation studies —
shorter APD at the apex, longer at the base.

## Purkinje conduction network

`constant/purkinjeGraph` (copied in by `Allrun` from `../mesh/constant/`)
is a biventricular Purkinje tree grown on this mesh by
`generatePurkinjeTree`, from UVC fields derived from `tm`/`tv`/`apicobasal`
(`uvc_transmural = 1-tm`, `uvc_intraventricular = 2*tv-1`,
`uvc_longitudinal = apicobasal`). Seeds were deduced from the AHA
segmentation `setCardiacAnatomy` computes, using the nearest LV/RV
endocardial point to the basal-septal AHA segments (`{2,3}` for LV, `21`
for RV). `electroProperties`'s `conductionNetworkDomains.purkinjeNetwork`
reads it via `graphFile purkinjeGraph;`.

## ECG electrodes

`ecgDomains.ECG.electrodePositions` are transferred to this anatomy from a
reference heart's validated V1-V6 placement, not measured on this mesh — a
normalized approximation, not patient-specific placement. 24-46mm from the
epicardium, comparable to the reference case's own spread.

## Execution

```bash
./Allrun
```
