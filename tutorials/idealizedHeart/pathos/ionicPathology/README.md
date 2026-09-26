# Ionic pathologies on the idealized biventricular heart

Monodomain tissue simulation modelling acute ischemia or Brugada syndrome
type-1 via `ionicConstantOverrides` on the `TNNP` ionic model, with
pseudo-ECG output and a coupled (unmodified) Purkinje network. Two variants
of the same case, selected at run time; no default — see Execution below.

## Stack

- electro model: `monodomainSolver`
- ionic model: `TNNP` (tissue); `BuenoOrovio` (Purkinje network — unchanged
  between variants)
- conduction system: `purkinjeGraphModel` with `monodomain1DSolver`,
  unmodified `graphFile purkinjeGraph` — unlike `../conductionBlock`,
  neither variant here touches the graph
- ECG: `ecgDomains.ECG.ecgSolver pseudoECG`

`constant/electroProperties` is a symlink, set by `Allrun`, to either
`electroProperties.ischemia` or `electroProperties.brugada` — the two
files are identical except for the `ionicConstantOverrides` block:

- **ischemia**: early ischemic zone (10-15 min) — `K_o` 5.4→8.0 mM
  (hyperkalemia), `g_Na`/`g_CaL` ×0.6, `g_K1` ×2.0, `P_NaK` ×0.5, applied
  globally. References: Ferrero et al. 1996, Trenor et al. 2010,
  Rodriguez et al. 2006.
- **brugada**: type-1 pattern — `g_Na` ×0.2 globally (SCN5A
  loss-of-function), `g_to` ×6.0 and `g_CaL` ×0.5 restricted to
  `epicardialCells` (I_to gain-of-function tips the plateau past the
  bifurcation point, causing epicardial dome loss and transmural
  dispersion of repolarisation). References: Coronel et al. 2009,
  Shimizu & Antzelevitch 1999.

Purkinje conduction and junction coupling: `purkinjeConductivity 0.4`
(~3.3 m/s along the tree) and `pvjRadius 1.65e-3`, the smallest junction
radius valid on this mesh.

## Tissue heterogeneity

`monodomainSolverCoeffs.ionicHeterogeneity` (`mode namedRegions;`) classifies
cells into `endocardialCells`/`mCells`/`epicardialCells` regions
(`range (0 0.3)`/`range (0.3 0.7)`/`range (0.7 1)`) from `field t;`,
`system/setExprFieldsDict`'s
`t = 1 - tm` (`tm`: 0 at epicardium, 1 at endocardium — the opposite
orientation `ionicHeterogeneity` requires). `epicardialCells` is also
where the `brugada` variant's `g_to`/`g_CaL` override applies.

## Mesh and anatomy fields

`Allrun` copies from `../../mesh/`:
`fiber`/`sheet`/`sheetNormal`/`tm`/`tv`/`apicobasal`/`Conductivity`/
`polyMesh`/`purkinjeGraph` as-is (no graph modification for either
variant here).

## ECG electrodes

`ecgDomains.ECG.electrodePositions` are placed by angle around the LV long
axis in the LV frame (`L` apex-to-base, `S` LV centre to RV centre, anterior
`A = L x S`, here `-z`): V1..V6 at 35, 65, 100, 135, 170, 205 deg from `S`
toward `A`, each at its original apex-base height and 25 mm from the nearest
tissue.

## Execution

```bash
./Allrun ischemia
./Allrun brugada
./Allrun ischemia parallel
```
