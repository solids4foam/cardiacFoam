# springSupportedSlab tutorial

## What it does

The Niederer et al. (2011) slab in electromechanics, with both fibre-wise
ends resting on spring supports (`solidSpringDashpot`). A plane wave of
monodomain TNNP activation along the fibres drives `LandNiedererTWorld`
active tension in a neo-Hookean solid, and the spring stiffness `kEnds`
decides how much the slab can shorten when it contracts. It writes the
displacement of and force on each end, and `Ta`, `Vm` and `Ca_i` probes.

- **Exercises:** the spring (Robin) support in the coupled electromechanics
  solver, exercised by no other case, including the check that the end force
  from the stress field equals `-kEnds · A0 · <Dx>`.
- **Use it for:** choosing spring stiffnesses for valve-plane or pericardial
  supports, since heart models are cut off at the valves and held by the
  pericardium and surrounding tissue, both usually modelled as springs; one
  number (the slab shortening) against one parameter (`kEnds`) makes the
  effect easy to read. The condition itself (formulation, implementation,
  literature values) is documented in
  `modules/solids4foam/src/solids4FoamModels/solidModels/fvPatchFields/solidSpringDashpot/README.md`.

## Regression

Tier: `full`. Markers: `needs-solids4foam` (exits 77 in `lightweight` build
mode).

`regression/regressionTest.sh` runs `./Allrun parallel` once (6 subdomains
per region) with the default `system/caseParameters` (`kEnds = 1e7`,
`deltaT = 1e-5`) but `endTime = 0.05` instead of the tutorial's 0.25: by then
activation has crossed the slab, `Ta` has risen in its middle and both ends
have moved about 30% beyond their resting preload. The stiffness sweep and
the `deltaT = 1e-6` reference are not part of it. It compares
`regression/springSupportedSlab.reference` (9 rows, values frozen from an
earlier run of the code):

- `Vm` at the far end, resting at 0.02 s and activated at 0.05 s (tolerance
  2e-3 and 3e-3 V)
- `Ta` near the stimulus (M band) and at the slab centre at 0.05 s (tolerance
  200 Pa)
- `Dx` on `xMin` at 0.02 and 0.05 s and on `xMax` at 0.05 s, the shortening
  (tolerance 7e-6 m)
- the spring law, `F_xMin = -kEnds · A0 · <Dx>` with `A0 = 2.1e-5 m²` and
  `kEnds` read from `system/caseParameters`, at 0.02 and 0.05 s (tolerance
  1e-4 N)

```bash
./regression/regressionTest.sh
```

## Folder structure

```text
tutorials/electromechanicsProtocols/springSupportedSlab/
├── 0/solid/
│   ├── D                       solidSpringDashpot ends (xMin, xMax), free lateral faces
│   ├── f0, f0f                 fibres along x
├── constant/
│   ├── physicsProperties       electroMechanicalModel
│   ├── electroMechanicalProperties
│   ├── electro/electroProperties
│   └── solid/                  solidProperties, mechanicalProperties, ...
├── system/
│   ├── caseParameters          kEnds, deltaT, endTime: the only swept entries
│   ├── controlDict             probes, end displacement and end force outputs
│   ├── blockMeshDict           patches xMin, xMax, lateral
│   ├── Taprobes, CaiProbes, Vmprobes
│   ├── electro/, solid/        region fvSchemes, fvSolution, decomposeParDict
├── setup/
│   ├── sweep_springStiffness.json      stiffness sweep and time-step reference
│   └── postprocess_springStiffness.py  summary table and figure
├── regression/
│   ├── regressionTest.sh
│   ├── springSupportedSlab.reference
│   └── full-only, needs-solids4foam
├── Allrun
└── Allclean
```

## Case setup

- **Geometry:** 20 × 3 × 7 mm slab, 40 × 6 × 14 hexahedra (0.5 mm), fibres
  along x.
- **Electrophysiology:** monodomain, TNNP, with endo/M/epi bands along z
  (`ionicHeterogeneity namedRegions` on `t = z/7 mm`). The stimulus covers
  the whole x = 0 layer (box `(0 0 0)`–`(1.5e-3 3e-3 7e-3)`), so activation
  travels along the fibres as a plane wave.
- **Coupling:** `sequentialElectroMechanical` with `LandNiedererTWorld`
  active tension.
- **Mechanics:** total Lagrangian, dynamic (`d2dt2 Euler`),
  `electroMechanicalLaw` with a neo-Hookean passive law (E = 100 kPa,
  ν = 0.3) and active stress along the fibres. The solid solve is converged
  tightly each step (`rTol 0.02`, `sTol 1e-7`). `aTol 1e-9` is kept below the
  per-step residual, which shrinks with `deltaT`, so that the relative
  criteria decide convergence at any time step.
- **Boundaries** (`0/solid/D`):
  - `xMin`, `xMax`: `solidSpringDashpot` with `kNormal = kTangential = $kEnds`,
    no damping
  - `lateral`: traction-free `solidTraction`

## Key dictionary scope

`system/caseParameters` holds every entry a study changes. It is read by
`system/controlDict` and `0/solid/D`:

```cpp
kEnds           1e7;      // end-spring stiffness, normal and tangential [Pa/m]
deltaT          1e-05;    // time step, electro and solid regions [s]
endTime         0.25;     // covers the active-tension peak [s]
```

- `kEnds` sets the support. The transition from isometric to free shortening
  is around 2E/L = 1e7 Pa/m.

| End springs | What the slab does |
|---|---|
| `kEnds → ∞` | ends held in place: isometric twitch, no shortening |
| `kEnds ≈ 2E/L` | springs as stiff as the tissue: partial shortening |
| `kEnds → 0` | ends free: unloaded shortening |

- `deltaT` sets the accuracy of the coupled solution. Lower it to check
  time-step convergence (the sweep includes `deltaT = 1e-6` as reference).

## Running

```bash
./Allrun              # serial, default kEnds = 1e7 Pa/m
./Allrun parallel     # 6 cores, system/decomposeParDict
```

Cost (OpenFOAM v2412, laptop): about 15–25 min serial and about 12 min on
6 cores for the 25,000 steps of the default case. A `deltaT = 1e-6` run costs
about 10× more.

### Stiffness sweep

`setup/sweep_springStiffness.json` lists the study cases as values of the
`system/caseParameters` entries:

| caseId | kEnds [Pa/m] | deltaT [s] |
|---|---|---|
| `k1e10` | 1e10 | 1e-5 |
| `k1e8` | 1e8 | 1e-5 |
| `k1e7` | 1e7 | 1e-5 |
| `k1e6` | 1e6 | 1e-5 |
| `k1e5` | 1e5 | 1e-5 |
| `k1e7_dt1e-6` | 1e7 | 1e-6 (time-step reference) |

Each case is a copy of this tutorial with those entries set, run with
`./Allrun parallel`. When the cases have finished, summarise them with:

```bash
python3 setup/postprocess_springStiffness.py \
    --input-dir <directory holding the case directories> \
    --output-dir <output directory>
```

This writes `springStiffness_summary.csv` and `springStiffness.png`.

## Outputs

| Output | Content |
|---|---|
| `postProcessing/solid/D_xMin/0/surfaceFieldValue.dat`, `D_xMax/...` | area-averaged displacement of each end |
| `postProcessing/0/solidForcesxMin.dat`, `solidForcesxMax.dat` | force on each end, integrated from the stress field |
| `postProcessing/Taprobes/solid/0/Ta` | active tension: transmural column at x = 0.75 mm (cols 2–4), slab centre and far end (cols 5–6) |
| `postProcessing/Vmprobes/electro/0/Vm` | transmembrane potential at x = 0.75, 10, 19.25 mm |
| `postProcessing/CaiProbes/electro/0/Ca_i` | intracellular calcium on the transmural column |

All probes and end outputs are written every 20 time steps. Fields are
written every 10 ms.

## How to check that the spring works

1. **Spring law.** The lateral faces are traction-free, so the only force
   on each end comes from its spring. The force on `xMin` integrated from
   the stress (`solidForces`) must equal `F = -kEnds · A0 · <Dx>` at every
   time. Here A0 = 3 × 7 mm² and `<Dx>` is the area-averaged end
   displacement. The two are computed independently: one from σ, the other
   from D.
2. **Limits.** A very stiff spring (1e10 Pa/m) must give an isometric twitch
   with an end force close to Ta_max · A0. A soft spring (1e5 Pa/m) must give
   the free shortening.
3. **Transition.** Shortening drops from free to zero around
   kEnds ≈ 2E/L = 1e7 Pa/m, where the two end springs are as stiff as the
   tissue.
4. **Time step.** The `deltaT = 1e-6` reference shows how far the default
   `1e-5` is from a time-converged answer.

## Results

The sweep of `setup/sweep_springStiffness.json`, OpenFOAM v2412, 250 ms, with
a solids4foam that includes `solidSpringDashpot` and the `electroMechanicalLaw`
deformation-gradient fix (#393):

| caseId | Peak shortening | Peak end force | Spring-law error* |
|---|---|---|---|
| `k1e10` | 0.03 % (isometric) | 575 mN | 9.9e-7 |
| `k1e8` | 2.19 % | 459 mN | 5.5e-6 |
| `k1e7` | 7.86 % | 165 mN | 1.2e-5 |
| `k1e7_dt1e-6` | 7.85 % | 160 mN | 3.1e-6 |
| `k1e6` | 10.76 % | 22 mN | 2.2e-5 |
| `k1e5` | 11.18 % (≈ free) | 2.3 mN | 2.3e-5 |

\* max \|F_stress − (−kEnds · A0 · ⟨Dx⟩)\| over the run, divided by
Ta_max · A0.

- **Spring law:** the end force from the stress field and the spring law
  agree to within 3e-5 of the active force scale for every stiffness.
- **Limits:** `k1e10` is isometric; its peak end force, 575 mN, is 98 % of
  Ta_max · A0 = 28.0 kPa × 21 mm² = 587 mN (Ta_max is the largest probe
  value; the end force integrates the whole cross-section). `k1e5` gives the
  free shortening.
- **Transition:** peak shortening drops from 11.2 % to 0 between 1e6 and
  1e10 Pa/m, with `k1e7` (kEnds = 2E/L) at 7.9 %.
- **Time step:** `deltaT = 1e-6` changes the peak shortening of `k1e7` by
  0.04 % and its peak end force (at 148 ms) by 3 %.
- **Activation:** the plane wave reaches the far end (x = 19.25 mm) at about
  45 ms, a conduction velocity of about 0.42 m/s along the fibres.
- **solids4foam version:** without #393 the active stress is added in the
  undeformed fibre direction (`Ta f0f0`) instead of pushed forward with the
  deformation gradient (`F (Ta f0f0) F^T / J`). The isometric case is then
  unchanged, but shortening is up to 10 % larger (12.4 % for `k1e5`).

`setup/postprocess_springStiffness.py` writes this table
(`springStiffness_summary.csv`) and a four-panel figure
(`springStiffness.png`): activation, shortening, spring-law check and peak
shortening against `kEnds`.

## Notes

- Use a very stiff `solidSpringDashpot` (1e10 Pa/m) for the isometric limit, not
  `fixedDisplacement` on both ends. With both ends fixed and a plane-wave
  stimulus, the displacement is exactly zero at the start, and the solid
  solver's relative-residual check cannot be met.
- Ta is about 5.7 kPa at t = 0, before any stimulation (the resting output of
  `LandNiedererTWorld`), so the springs are loaded from the first step.
- The twitch (about 200 ms) is slow compared with the spring–mass periods,
  so the response is nearly quasi-static. A dashpot (`cNormal`,
  `cTangential`) has little effect here.

## Notes on dictionary entries

**Mesh.** `system/blockMeshDict` uses 40 x 6 x 14 hexahedra on the
20 x 3 x 7 mm slab, a uniform spacing of 0.5 mm.

**Active tension.** `0/solid/Ta` is the initial active tension. The
electromechanical model registers and solves `Ta`, but the solid model builds
its stiffness while it is constructed, before that field exists, so
`electroMechanicalLaw` reads `Ta` from this file.
`constant/electroMechanicalProperties` selects the scalar `LandNiedererTWorld`
model, integrated with the `ODESolver` settings (`RKF45`, `absTol 1e-8`,
`relTol 1e-6`). The batched versions of the active-tension models are listed
in [`src/activeTensionModels/README.md`](../../../src/activeTensionModels/README.md).

**Solid solver.** In `constant/solid/solidProperties`, `rTol` is the relative
residual (converged when residual < rTol x initial residual), `aTol` the
absolute residual (kept far below the per-step residual so that `rTol` and
`sTol` decide convergence at any `deltaT`), and `sTol` the step norm
(converged when deltaX < sTol x xNorm).

**Probes.** `Taprobes` and `CaiProbes` share the transmural column under the
stimulus footprint at x, y = 0.75 mm: probe 0 (col 2) endo band centre,
t = 0.15, z = 1.05 mm; probe 1 (col 3) M band centre, t = 0.50, z = 3.50 mm;
probe 2 (col 4) epi band centre, t = 0.85, z = 5.95 mm. `Taprobes` adds probes
3-4 (cols 5-6) at the slab centre and far end along the fibres (x = 10.0 and
19.25 mm; y = 1.5 mm, z = 3.5 mm). `Vmprobes` samples the transmembrane
potential along the fibres at y = 1.5 mm, z = 3.5 mm: stimulated layer, slab
centre and far end (x = 0.75, 10.0 and 19.25 mm).

**Controls.** `system/controlDict` records the area-averaged displacement of
each spring-supported end, and the force on each end from the stress field,
independent of the spring law. The top-level `system/fvSchemes` and
`system/fvSolution` are dummies for the parent region and are not used
directly.
