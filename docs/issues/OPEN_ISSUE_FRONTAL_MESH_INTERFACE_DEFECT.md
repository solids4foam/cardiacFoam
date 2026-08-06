# Open issue: the Frontal/Netgen tet family has a localized, boundary-condition-independent defect at the myocardium-bath interface

**Status: open. Root cause not found. Two independent fix attempts made it
worse, not better. Do not use the Frontal/Netgen mesh family for bath-bidomain
production or paper results until this is understood.**

This is a *different* problem from `OPEN_ISSUE_N80_interface_instability.md`
(the generic-Delaunay N=80 groundElectrode interface asymmetry). That issue
is specific to one mesh, one resolution, one late-onset time threshold. This
one is specific to one *mesh generation algorithm*, reproducible at N=40
already, appears immediately (no decay-then-blowup threshold), and is
completely independent of the electrode boundary condition.

## Why this matters

The Frontal/Netgen family (`Mesh.Algorithm3D=4` + `OptimizeNetgen=1` +
`Smoothing=100`) is not an arbitrary alternative mesh. Per
`14_verification_methodology.qmd` in the paper, it is deliberately calibrated
to match the mesh-quality statistics of real anatomical hearts (mean
non-orthogonality ~14.4-15.8° for the Frontal family vs ~14° for ten real
four-chamber hearts, vs 18.6° for the generic Delaunay family used everywhere
else in the bath-bidomain study). `17_discussion.qmd` explicitly calls N=80
"a mesh resolution typical of anatomical studies (0.125mm)". The Frontal
family was never previously run on the bath-bidomain operator at all --
`15_results.qmd`: "Bidomain and pseudo-ECG were not repeated on the
Frontal/Netgen family." This is the first time it has been.

## The result across four independently generated meshes, N=40, electrodePair unless noted

| mesh | algorithm | cells/dir | mean non-ortho | Vm L2 | Vm Linf | phiE L2 | phiE Linf | interface behavior |
|---|---|---|---|---|---|---|---|---|
| Delaunay (paper default) | `Algorithm3D=1` | 66 | 18.6° | 6.98e-5 | 4.88e-4 | 2.74e-4 | 7.49e-4 | healthy, symmetric |
| HXT (third, unrelated algorithm) | `Algorithm3D=10` | 63 | 22.2° | 8.00e-5 | 5.87e-4 | 3.02e-4 | 7.96e-4 | healthy, symmetric |
| **Frontal (original)** | `Algorithm3D=4` | 70 | ~14.4-15.8° | 5.36e-5 | **4.28e-3** | 1.83e-4 | **4.91e-3** | localized, ~39% x0/x1 asymmetry in assembled flux |
| **Frontal (2nd instance, perturbed lc)** | `Algorithm3D=4` | 72 | ~14.4-15.8° | 5.14e-5 | **1.12e-3** | 1.76e-4 | **1.27e-3** | localized, ~11% x0/x1 asymmetry -- milder, not absent |

Note the `L2` (bulk/average) column tracks mean non-orthogonality exactly as
expected on its own -- worse mean quality gives worse bulk error, monotonic
across all four meshes (HXT worst quality, worst L2; Frontal best quality,
best L2). That relationship is uninteresting and is not what this document
is about. The interesting, and separate, effect is `Linf` (worst-case,
localized): it does **not** track mean quality at all -- Frontal has the
*best* mean non-orthogonality of the four and the *worst* Linf by nearly an
order of magnitude. Comparing HXT's healthy Linf against Frontal's
pathological one is not by itself a quality-matched control, since HXT's
mean quality is also the worst of the four -- see "Mesh-quality-matched
comparison attempt" below for why that matters and what was (and was not)
achieved trying to close that gap.

Delaunay and HXT are both clean at this resolution. Two independently
generated Frontal meshes (same algorithm, different characteristic length,
confirmed genuinely different meshes by cell count and by the fact that
`Mesh.RandomSeed` was tried first and does **not** perturb this pipeline --
see "Fix attempts" below) both show elevated, localized, asymmetric error.
Severity varies by instance; the defect does not.

## The defect is completely boundary-condition-independent

On the same original Frontal mesh (70 cells/dir), groundElectrode and
electrodePair were run and compared directly:

| metric | electrodePair | groundElectrode |
|---|---|---|
| Vm Linf | 4.27968e-3 | 4.27968e-3 (identical to full precision) |
| x0 intracellular leak L2 | 5.5952e-05 | 5.60338e-05 |
| x1 intracellular leak L2 | 7.2553e-05 | 7.16483e-05 |
| x0 assembled flux L2 | 1.9316e-04 | 1.93179e-04 |
| x1 assembled flux L2 | 2.6864e-04 | 2.68574e-04 |
| x0 flux jump L2 | 9.2733e-05 | 9.27383e-05 |
| x1 flux jump L2 | 1.2478e-04 | 1.24591e-04 |

Every column matches to 3-4 significant figures. `Vm` is identical to full
solver precision -- and `Vm`'s exact solution carries no gauge/reference-cell
freedom at all (only `phiE`/`phiI` do, and only for `electrodePair`), so this
cannot be an artifact of the electrode boundary condition or of the
mean-removal gauge correction. It is a property of the (mesh x operator)
combination, full stop.

## Localization

`bathBidomainInterfaceMetrics -latestTime -writeFaceErrors` on the original
Frontal mesh (N=40, electrodePair): one face holds 29.3% of total error
energy across 7420 interface faces; the top 10 hold 40.3%. Concentrated at
interface **x=1** (not x=0 -- see below), clustered near y in [0.12, 0.68], z
in [0.007, 0.09], i.e. near the z=0 domain edge.

**Local mesh quality does not distinguish the flagged region.** Sampled
`checkMesh -writeFields '(nonOrthoAngle skewness)'` for the ~4600 cells in
the flagged bounding box against the whole-mesh population (1,048,451
cells): flagged-region mean non-orthogonality 21.3° vs population mean
19.7°; flagged-region mean skewness 0.182 vs population mean 0.175;
flagged-region max values (49.6°, 0.448) do not exceed the population max
(49.6°, 0.505). Statistically indistinguishable -- the same conclusion
`OPEN_ISSUE_N80_interface_instability.md` reached refuting the mesh-quality
hypothesis for its own (different) case. The usual local-quality
diagnostics do not explain this defect on this mesh either.

**Which interface misbehaves is not fixed.** groundElectrode/Delaunay/N=80's
documented instability is at x=0. This defect, on a completely different
mesh generator and variant, is at x=1. That argues against any explanation
tied to a specific boundary condition or a specific interface's geometry --
something about the Frontal-advancing algorithm's handling of an internal,
fragmented boundary (`BooleanFragments` merging three separately-swept boxes)
appears able to produce a bad local transition at *either* shared face,
depending on the specific mesh instance.

## Fix attempts -- both made it worse, not better

**Attempt 0: different `Mesh.RandomSeed`.** Does not perturb this
`Algorithm3D=4` + `OptimizeNetgen` pipeline at all -- confirmed bit-identical
mesh (same cell count, error values matching to full solver precision) with
`RandomSeed=42` vs the unset default. A different, genuinely distinct second
instance required perturbing the characteristic length instead
(`lc * 1.000037`), which is the "2nd instance" row in the table above.

**Attempt 1: aggressive local refinement at both internal interfaces.**
Gmsh `Distance` + `Threshold` mesh-size fields, `SizeMin = lc/2`,
`SizeMax = lc`, transition width `1*lc`, targeting the x=0 and x=1 internal
surfaces (found via `Surface In BoundingBox`, deliberately not made
`Physical` so `gmshToFoam` still imports them as internal faces). Mesh
generation and `checkMesh` both succeeded. The **solve crashed**: `phiE`'s
PCG initial residual grew every timestep (0.146 -> 0.277 -> 0.789 -> 0.896
over four steps, each converging to tight tolerance but diverging further
between steps), ending in a `SIGFPE` inside the ionic model's ODE solver
(`bathBidomainFDAManufactured::derivatives`). The abrupt size gradation
likely produced worse-shaped transition tets than the uniform mesh it was
meant to improve on.

**Attempt 2: gentler local refinement.** Same interfaces, milder
`SizeMin = 0.75*lc`, wider transition `2.5*lc`, `Sigmoid=1` blending. Mesh
generation and `checkMesh` again succeeded, and the solve did not crash --
but `Vm_Linf` reached **0.0515** and `phiE_Linf` reached **0.0676**, roughly
14x worse than the unrefined Frontal mesh's already-bad 4.28e-3/4.91e-3, with
`x0` degrading far more than `x1` (x0 intracellular leak 7.09e-4 vs x1
3.81e-5 -- an 18x asymmetry, vs ~1.3x on the unrefined mesh). The refined
interface targeted both x=0 and x=1 symmetrically in the script; the
generator did not treat them symmetrically in practice.

**Conclusion: do not keep trying ad hoc Gmsh size-field tuning on this
geometry.** Two principled attempts, in opposite directions (aggressive vs.
gentle gradation), both made the defect worse rather than better, one by
crashing the solve outright and one by amplifying the exact asymmetry being
targeted. That is not noise around a fix that is nearly right; it is evidence
the naive local-refinement approach is fighting something structural in how
this generator handles the fragmented internal boundary, not a tuning
problem solvable by adjusting the same two size-field parameters again.

## Mesh-quality-matched comparison attempt: can HXT be tuned into the anatomical range?

HXT's healthy Linf is not, by itself, a fair same-quality control against
Frontal -- per the table above, HXT's own mean non-orthogonality (22.2°) is
the *worst* of the four meshes tested, well outside the anatomical target
(~14°) Frontal was calibrated against. A mesh that is simply coarser
everywhere would be expected to have worse bulk error (it does, see L2
above) without necessarily reproducing Frontal's specific localized defect.
The open question this raises: if HXT is tightened toward that same
anatomical quality range, does it develop the same pathology, or does it stay
clean at matched quality -- which would be much stronger evidence the defect
is Frontal-algorithm-specific rather than a property of any mesh that
anatomically realistic.

**This was attempted and did not succeed -- not because it disproved
anything, but because HXT's quality could not be moved.** Two attempts:

1. Apply the same post-generation optimization Frontal uses
   (`Mesh.OptimizeNetgen=1; Mesh.Smoothing=100;`) on top of HXT-generated
   (`Algorithm3D=10`) topology, rather than changing the generation
   algorithm. Result: mean non-orthogonality 22.2° -> 22.2° (no change), max
   skewness 0.929 -> 0.989 (slightly worse). `Vm_L2`/`Linf` and
   `phiE_L2`/`Linf` all within noise of the untuned HXT mesh.
2. Force explicit, repeated optimization passes after generation
   (`Mesh 3; For i In {1:5} OptimizeMesh "Netgen"; EndFor`) rather than
   relying on the option flags. Result: **identical to six significant
   figures** to attempt 1 (mean non-orthogonality 22.1976° both times, max
   skewness 0.988817 both times, error norms matching to 4+ significant
   figures). Netgen's smoother converges to the same local optimum after one
   pass; five more do not move it.

**Conclusion: `OptimizeNetgen`/`Smoothing` are effective on meshes Netgen's
own Frontal algorithm generates, and are not effective at improving the
quality of a mesh whose initial topology came from a different generator
(HXT).** This is a plateau, not a partial result -- two attempts at
different intensities landed on the same number. Reaching the anatomical
quality range on HXT topology would need a fundamentally different
technique (e.g. mesh simplification/decimation from a much finer initial
mesh, rather than vertex-smoothing the coarse one), which was not attempted.

**Net effect on the question this section opened with:** unresolved, not
negative. The only mesh generator tested that actually reaches the
anatomically-representative quality range is Frontal, and it has the
defect. Whether a genuinely quality-matched *and* defect-free tet mesh is
achievable for this geometry with some other technique remains open -- this
attempt narrows the "how" (post-generation optimization on non-Netgen
topology is not it) without answering the "whether".

## What this does and does not tell us

- Does not affect or explain `OPEN_ISSUE_N80_interface_instability.md`. That
  remains open, on the Delaunay mesh, specific to N=80.
- Does not implicate the electrodePair variant, the mean-removal gauge fix,
  or `bathBidomainInterfaceMetrics.C`. All confirmed unaffected by direct
  evidence (`Vm` invariance, gauge-free-field argument above).
- Does implicate the Frontal/Netgen mesh generation of this specific
  three-box-fragmented geometry, reproducibly across two independent
  instances, at a resolution (N=40) well below the anatomically-relevant
  N=80 the paper otherwise targets for this family.
- N=80 on the Frontal family was deliberately **not** run for this
  investigation -- each attempt costs several hours, and there is no reason
  to expect refinement to help when it already does not at N=40; the defect
  would very plausibly be at least as bad, following the same pattern as the
  Delaunay N=80 issue (finer resolution, worse outcome).

## Recommended next steps (not attempted here)

1. **Look from outside the solver.** This may need inspecting the Gmsh/Netgen
   mesh directly (Gmsh GUI or a dedicated quality-visualization pass) at the
   flagged region to see the actual element shapes, not just aggregate
   non-orthogonality/skewness statistics that failed to distinguish anything.
2. **Consider a different geometry-construction strategy** for the fragmented
   interface, e.g. meshing the three volumes with explicit conformal
   boundary matching *before* invoking the Frontal algorithm, rather than
   `BooleanFragments` + a single global Frontal pass.
3. **Consider a different coupling/assembly strategy at the interface**
   itself -- if the matched-submesh assembly is sensitive to a specific class
   of transition element the Frontal generator produces, a coupling scheme
   less sensitive to local element shape (e.g. one that doesn't rely on the
   one-sided quadratic face-fit `oneSidedFaceFit` this diagnostic and the
   solver both use) might be worth comparing.
4. Report the specific bad-element characteristics to Gmsh/Netgen upstream if
   step 1 identifies a generator-side degenerate case, since this is
   reproducible and fairly minimal (a 3-box `BooleanFragments` geometry).

## Provenance

Solver repo `noFrontendCardiacFoam_minor_errors`, branch
`no-frontend-minor-errors`, 2026-08-05. Meshes and sweep specs used:

- `setup/mesh/tet/three_domain_box.geo.template.optimised` -- the original
  Frontal family (`Algorithm3D=4`).
- `setup/mesh/tet/three_domain_box.geo.template.hxt` -- third-algorithm
  control (`Algorithm3D=10`).
- `setup/mesh/tet/three_domain_box.geo.template.frontal_seed2` -- second
  Frontal instance (perturbed `lc`, not `RandomSeed`).
- `setup/mesh/tet/three_domain_box.geo.template.frontal_refined_interface`
  -- Attempt 2's interface-refined field. Attempt 1's more aggressive
  version was not preserved as a file; its parameters are recorded above.
- `setup/mesh/tet/three_domain_box.geo.template.hxt_optimised` -- HXT
  generation plus the repeated-optimization-pass attempt described above.
  Kept even though the optimization had no effect, as evidence of what was
  tried (option flags plus explicit `For`-loop `OptimizeMesh` calls, both
  landing on the same plateau).
- Diagnostic-only sweep specs and their raw results are consolidated under
  `setup/studies/tetConvergence/frontalMeshInvestigation/` (this
  investigation's growth-curve, seed2, groundElectrode-frontal,
  refined-interface, hxt, and hxt-optimised specs) rather than left in the
  main study directory, since their results are now fully captured in this
  document.
- The primary electrodePair-vs-Delaunay-vs-Frontal N=10/20/40 comparison
  remains in `setup/studies/tetConvergence/sweep_tet_electrodePair.json` and
  `sweep_tet_electrodePair_frontal.json` (main study directory, not moved --
  these are the primary result this file extends, not one-off diagnostics).

To reproduce the healthy/pathological comparison at N=40:

```bash
python3 -m openfoam_driver sweep-run \
  --spec tutorials/manufacturedSolutions/bathBidomain/setup/studies/tetConvergence/sweep_tet_electrodePair.json \
  --output-dir <scratch-dir> --case-timeout-s 1800
# then with tet_geo_template_relpath overridden to three_domain_box.geo.template.optimised
```
