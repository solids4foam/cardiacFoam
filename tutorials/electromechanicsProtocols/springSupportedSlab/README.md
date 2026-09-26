# springSupportedSlab tutorial

The Niederer et al. (2011) slab in electromechanics, with both fibre-wise
ends resting on spring supports (`solidRobin`). The spring stiffness decides
how much the slab can shorten when it contracts:

| End springs | What the slab does |
|---|---|
| `kEnds → ∞` | ends held in place: isometric twitch, no shortening |
| `kEnds ≈ 2E/L` | springs as stiff as the tissue: partial shortening |
| `kEnds → 0` | ends free: unloaded shortening |

## Why this case exists

Heart models are cut off at the valves and are held in the body by the
pericardium and surrounding tissue. Both are usually modelled as spring
(Robin) supports instead of fixed or free boundaries. This case checks the
spring support in the coupled electromechanics solver on a geometry where
the answer is easy to read: one number (the slab shortening) against one
parameter (`kEnds`). The condition itself (formulation,
implementation, literature values) is documented in
`modules/solids4foam/src/solids4FoamModels/solidModels/fvPatchFields/solidRobin/README.md`.

## Folder structure

```text
tutorials/electromechanicsProtocols/springSupportedSlab/
├── 0/solid/
│   ├── D                       solidRobin ends (xMin, xMax), free lateral faces
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
│   └── springSupportedSlab.reference
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
  - `xMin`, `xMax`: `solidRobin` with `kNormal = kTangential = $kEnds`,
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

The sweep of `setup/sweep_springStiffness.json`, OpenFOAM v2412, 250 ms:

| caseId | Peak shortening | Peak end force | Spring-law error* |
|---|---|---|---|
| `k1e10` | 0.03 % (isometric) | 574 mN | 1.0e-6 |
| `k1e8` | 2.25 % | 472 mN | 7.4e-6 |
| `k1e7` | 8.59 % | 180 mN | 1.5e-5 |
| `k1e7_dt1e-6` | 8.58 % | 175 mN | 3.9e-6 |
| `k1e6` | 11.93 % | 25 mN | 2.9e-5 |
| `k1e5` | 12.41 % (≈ free) | 2.5 mN | 3.1e-5 |

\* max \|F_stress − (−kEnds · A0 · ⟨Dx⟩)\| over the run, divided by
Ta_max · A0.

- **Spring law:** the end force from the stress field and the spring law
  agree to within 3e-5 of the active force scale for every stiffness.
- **Limits:** `k1e10` is isometric; its peak end force, 574 mN, is 98 % of
  Ta_max · A0 = 27.8 kPa × 21 mm² = 584 mN (Ta_max is the largest probe
  value; the end force integrates the whole cross-section). `k1e5` gives the
  free shortening.
- **Transition:** peak shortening drops from 12.4 % to 0 between 1e6 and
  1e10 Pa/m, with `k1e7` (kEnds = 2E/L) at 8.6 %.
- **Time step:** `deltaT = 1e-6` changes the peak shortening of `k1e7` by
  0.07 % and its peak end force (at 148 ms) by 3 %. The shortening curves
  differ mainly in the first few milliseconds, where the step load from the
  resting Ta makes the `deltaT = 1e-5` run ring briefly.
- **Activation:** the plane wave reaches the far end (x = 19.25 mm) at about
  45 ms, a conduction velocity of about 0.42 m/s along the fibres.

`setup/postprocess_springStiffness.py` writes this table
(`springStiffness_summary.csv`) and a four-panel figure
(`springStiffness.png`): activation, shortening, spring-law check and peak
shortening against `kEnds`.

## Notes

- Use a very stiff `solidRobin` (1e10 Pa/m) for the isometric limit, not
  `fixedDisplacement` on both ends. With both ends fixed and a plane-wave
  stimulus, the displacement is exactly zero at the start, and the solid
  solver's relative-residual check cannot be met.
- Ta is about 5.7 kPa at t = 0, before any stimulation (the resting output of
  `LandNiedererTWorld`), so the springs are loaded from the first step.
- The twitch (about 200 ms) is slow compared with the spring–mass periods,
  so the response is nearly quasi-static. A dashpot (`cNormal`,
  `cTangential`) has little effect here.

## Regression behavior

The regression is a single run: `regression/regressionTest.sh` runs
`./Allrun parallel` once with the default `system/caseParameters`
(`kEnds = 1e7`, `deltaT = 1e-5`, `endTime = 0.25`; about 12 min on 6
cores). The stiffness sweep and the `deltaT = 1e-6` reference are not part
of it. The run is compared against `regression/springSupportedSlab.reference`:

- `Vm` at the far end (activation has crossed the slab)
- `Ta` at the slab centre
- `Dx` on `xMin` and `xMax` (the slab shortening)
- the spring law, `F_xMin = -kEnds · A0 · <Dx>`, at the same times

It exits 77 (expected skip) in `lightweight` build mode, and is wired into
`tutorials/Alltest-regression`.
