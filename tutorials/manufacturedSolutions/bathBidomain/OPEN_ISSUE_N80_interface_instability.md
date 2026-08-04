# Open issue: localised interface blow-up at N=80, x=0 (bath-bidomain, tetrahedral)

**Status:** unresolved. Two hypotheses tested and refuted. A third is being measured.
**Impact on Paper I:** none if left alone. The manuscript reports the bath ladder
over N = 10, 20, 40 and explicitly excludes N = 80, so no published number depends
on this. Do not "fix" the paper — fix or explain the solver behaviour.

---

## One-paragraph summary

On the conformal tetrahedral bath-bidomain manufactured case at N = 80, the
interface diagnostics at the **x = 0** interface blow up, while the **x = 1**
interface on the same mesh, in the same run, converges normally. The blow-up is
confined to **22 faces out of 14,760**, one of which holds 87% of the total error
energy. Everything static checks out — mesh, decomposition, interface operator,
linear solves — and a 2-step run of the same case is clean. The defect therefore
**develops over the time integration**.

---

## Reproducing it

```bash
cd tutorials/manufacturedSolutions/bathBidomain
RESOLUTIONS=80 FACE_ERRORS=1 ./setup/mesh/tet/run_bath_tet_predictor.sh
```

Roughly 3.3 h on 6 ranks (143 steps, ~81 s/step, 6.83 M cells). The mesh is built
by `setup/mesh/tet/run_mesh_gate.sh 80` and banked at
`setup/mesh/tet/interfaceMeshBank/N80/polyMesh.tar.gz` (210 MB), so meshing need
not be repeated.

Outputs land in
`setup/mesh/tet/interfaceStudy/matchedSubmesh/distanceWeightedHarmonic_predictor/N80/`:

| file | contents |
|---|---|
| `bathBidomainInterfaceMetrics.csv` | solved-field interface norms |
| `bathBidomainInterfaceMetricsExactField.csv` | same, with the manufactured phiE substituted |
| `bathBidomainFaceErrors.csv` | **per-face signed errors with coordinates** — the key file |
| `log.checkMesh`, `log.cardiacFoam` | mesh and solver logs |

`FACE_ERRORS=1` is what produces the per-face dump. Without it you only get norms
and the localisation is invisible.

---

## The evidence

### The blow-up is 22 faces, not a field

From `bathBidomainFaceErrors.csv` at N = 80:

| | x = 0 | x = 1 |
|---|---|---|
| faces | 14,760 | 14,760 |
| RMS error | 1.397e-2 | 2.852e-4 |
| **median \|error\|** | **1.634e-4** | **1.639e-4** |
| max \|error\| | 1.583e+00 | 1.096e-03 |
| top 1 face holds | **87.0%** of error energy | 0.1% |
| top 10 faces hold | **99.8%** | 0.9% |

The medians are identical. The x = 0 interface is healthy everywhere except a
patch of 22 faces with |error| above the physical scale (alpha = 0.01 A/m²).

Those faces occupy **y in [0.0375, 0.0938], z in [0.0038, 0.0399]** — about 0.2%
of an interface spanning y, z in [0.003, 0.997]. Their areas are 1.09–1.20x the
median face area, i.e. **not slivers**.

### The convergence sequences

Effective spacing h = n_cells^(-1/3); cell counts 4874 / 37040 / 289388 / 2278889.

| quantity | N=10 | N=20 | N=40 | N=80 |
|---|---|---|---|---|
| x0 assembled current | 4.442e-4 | 2.814e-4 | 2.283e-4 | **1.514e-2** |
| x1 assembled current | 4.677e-4 | 3.012e-4 | 2.402e-4 | 2.851e-4 |
| x0 intracellular leak L2 | 3.382e-3 | 5.364e-4 | 1.352e-4 | **3.556e-3** |
| x0 intracellular leak Linf/L2 | 1.9 | 1.8 | 2.9 | **38.2** |
| x1 intracellular leak L2 | 3.323e-3 | 5.398e-4 | 1.351e-4 | 5.361e-5 |
| bath potential | 4.662e-3 | 1.517e-3 | 4.724e-4 | 1.783e-4 |

x = 1 keeps converging through the fourth level. The bath potential converges at
1.42 across it. Only x = 1's counterpart at x = 0 misbehaves.

**Intracellular insulation is violated**: the leak should be identically zero, and
its Linf reaches 0.136 — 13.6x the physical current scale — where at N = 40 it was
3.88e-4.

### The exact-field control separates operator from solve

`-exactFields` substitutes the manufactured phiE in cells and boundary faces
before the diagnostics run. The assembled-current probe is then clean at all four
levels and at both interfaces:

| | N=10 | N=20 | N=40 | N=80 | orders |
|---|---|---|---|---|---|
| x0 exact-field | 1.607e-3 | 8.850e-4 | 4.422e-4 | 2.522e-4 | 0.88, 1.01, 0.82 |
| x1 exact-field | 1.640e-3 | 8.685e-4 | 4.519e-4 | 2.515e-4 | 0.94, 0.95, 0.85 |

At N = 80 the two interfaces agree to three significant figures. **The interface
operator on that mesh is sound.** Note the exact-field run does *not* overwrite
phiI, so its intracellular-leak column is still contaminated and should not be
used as a control for that quantity.

---

## What has been ruled out

| hypothesis | test | verdict |
|---|---|---|
| Mesh quality / slivers | `checkMesh`: Mesh OK, max non-orth 70.40 (only 3 faces > 70), max skewness 0.918. Bad faces are 1.09–1.20x median area | **refuted** |
| Domain decomposition, processor-patch coefficient injection | Mapped all 22 faces through `faceProcAddressing` for all 6 processors: **0 of 22** lie on a processor patch; only 3 of 14,760 interface faces touch one at all | **refuted** |
| Algebraic / tolerance | phiE reaches the 1e-15 absolute criterion in ~1200 iterations at every step; Vm in 4–5; no solver diagnostic raised anywhere in the log | **refuted** |
| Interface discretisation | exact-field probe clean and symmetric at N = 80 (table above) | **refuted** |
| Static assembly error | 2-step run at N = 80 is clean: x0 assembled 2.464e-4, x1 2.514e-4, leak Linf/L2 = 4.4. Artefacts in `setup/mesh/tet/interfaceStudy/solveDiagnostic/N80_2step/` | **refuted** |
| Ill-posed manufactured solution | sigma_i is configured as 0.111453302, which is 1.1/pi^2 to nine decimal places, so there is no parameter drift. The Cartesian bath at N = 80 gives a uniform 1.24e-6 at x = 0 and 2.09e-5 at x = 1 -- identical on all 6400 faces, four orders below the physical scale alpha = 1e-2. An inconsistent manufactured description would show a large error there, not a negligible one. Artefacts in `setup/mesh/hex/results/interfaceCartesian/N80/` | **refuted** |

One caveat on the Cartesian control: the manufactured phiE depends on x alone,
so every face of a Cartesian interface is geometrically equivalent and that run
*could not* have shown localisation. It bounds the magnitude, not the spatial
structure. It also settles a smaller question -- the Cartesian solved value's
rise at N = 80 is 1.2434e-6 and uniform, i.e. the tolerance floor, not a
miniature of the tetrahedral effect.

The decomposition test is worth keeping; the script is in the git history of this
investigation and re-derivable from `processorN/constant/polyMesh/faceProcAddressing`
plus the `procBoundary*` entries in each `boundary` file.

---

## Leading hypothesis

**Something in the coupled bath advance crosses a threshold between step 40 and
step 143 and then fails locally.** The predictor–corrector is the first suspect
because it is the only outer iteration in the advance, but note that the growth
curve below argues against the simplest version of that story.

What the evidence supports:

- The failing quantity is intracellular insulation, i.e. the coupled side, not
  the extracellular operator.
- The linear solves are converged throughout, so this is the outer coupling, not
  the inner solve.
- The failure is spatially confined, which fits a local rather than a global
  mechanism.

What it argues against:

- A conditionally convergent fixed-point iteration should degrade monotonically
  from early on. It does not: the error *decays* over the first 40 steps and the
  two interfaces stay symmetric there. Any mechanism has to explain 40 steps of
  improvement followed by a blow-up, not just the endpoint.

Against a pure step-size explanation: dt shrinks by exactly 4.0 per refinement
level (`dt_for_n` in `run_bath_tet_predictor.sh`: 8.92857e-3, 2.24215e-3,
5.60538e-4, 1.401345e-4) while the true h² shrinks by 3.96, so N = 80 is
marginally *more* conservative than N = 40, not less. A naive diffusive stability
limit does not explain it, and the diffusion is implicit in any case.

---

## Next tests, in order

1. **Growth curve — DONE. Onset lies between step 40 and step 143.**

   | steps | x0 leak L2 | x0 leak Linf | x1 leak L2 | x0 assembled flux |
   |---|---|---|---|---|
   | 2 | 1.245e-4 | 5.509e-4 | 1.236e-4 | 2.464e-4 |
   | **40** | **8.127e-5** | **3.146e-4** | **8.109e-5** | **2.739e-4** |
   | 143 | 3.556e-3 | 1.359e-1 | 5.361e-5 | 1.514e-2 |

   At step 40 the case is not merely healthy, it is *improving* — the leak falls
   from 1.245e-4 to 8.127e-5 — and the two interfaces are symmetric to three
   significant figures (8.127e-5 against 8.109e-5). There is no slow drift from
   the start.

   This matters for the mechanism. A conditionally convergent fixed-point
   iteration would degrade from early on. Decay for 40 steps followed by a
   blow-up by step 143 looks instead like a threshold being crossed. Whatever is
   proposed has to reproduce that shape, not just the endpoint.

   Artefacts: `setup/mesh/tet/interfaceStudy/solveDiagnostic/N80_40step/`.

2. **Corrector off — RUNNING at time of writing.** `bath_predictor_corrector`
   is a first-class driverFOAM config key, so this is a sweep spec, not a script:
   ```bash
   driverFoam sweep-run \
     --spec setup/studies/tetConvergence/sweep_tet_generic_correctorOff.json \
     --output-dir setup/studies/tetConvergence/results/sweepRunCorrectorOff
   ```
   ~3.3 h. The spec is a copy of `sweep_tet_generic.json` differing only in
   `bath_predictor_corrector: false` and the archive directory.
   If the baseline one-pass coupling is stable at N = 80, the corrector is the
   mechanism. If it blows up too, the corrector is exonerated and the cause lies
   in the shared bath assembly, which would redirect the search entirely.

3. **Halved timestep.** N = 80 at dt = 7.0e-5 to the same endTime (286 steps,
   ~6.6 h). If the instability disappears at the same physical time, it is a
   step-size stability limit and the coupling needs a stability condition stated.

4. **Localise in the field, not just on the interface.** The per-face dump gives
   the interface trace only. Writing the full phiE/phiI/Vm fields at the onset
   step and inspecting a neighbourhood of (x=0, y~0.06, z~0.02) would show whether
   the corruption is a surface phenomenon or a volume one leaking to the surface.

---

## Where things live

- Case: `tutorials/manufacturedSolutions/bathBidomain/`
- Ladder runner: `setup/mesh/tet/run_bath_tet_predictor.sh` (`FACE_ERRORS=1` for the per-face dump)
- Sweep specs: `setup/studies/tetConvergence/sweep_tet_generic.json` (production)
  and `sweep_tet_generic_correctorOff.json` (corrector disabled). Run either with
  `driverFoam sweep-run --spec <json> --output-dir <dir>`. Truncate a run by
  lowering `end_time` in a copy of the spec rather than writing a shell wrapper.
- Mesh gate: `setup/mesh/tet/run_mesh_gate.sh <N>`
- Metrics utility: `applications/utilities/bathBidomainInterfaceMetrics/bathBidomainInterfaceMetrics.C`
  (`-exactFields`, `-writeFaceErrors`)
- Matched-submesh injection: `src/electroModels/electroDomains/extracellularPotentialDomain/extracellularPotentialDomain.C`,
  around line 920 — the loop that adds heart-submesh boundary coefficients into
  the base-mesh equation. Inspected during this investigation and no defect found,
  but it is the code that owns interface insulation and worth re-reading with the
  onset step in hand.
- Archived data for the paper:
  `data/paperI/solver/tutorials/manufacturedSolutions/bathBidomain/` in the
  manuscript repository, including `exactField/` with all four levels.

## Do not

- Do not widen the paper's bath ladder to N = 80 before this is understood.
- Do not read the N = 80 point as a property of the discretisation; the exact-field
  control on the same mesh shows the operator is first order and symmetric there.
- Do not attribute it to weakening error cancellation. That mechanism predicts the
  solved error rising *to* the exact-field error, which is exactly what x = 1 does
  (ratio 0.285, 0.347, 0.531, 1.13). At x = 0 the ratio reaches 60, and lost
  cancellation cannot push the solved error above the operator's own consistency
  error.
