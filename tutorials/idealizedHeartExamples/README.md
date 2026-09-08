# Idealized heart examples

All cases in this directory share the same idealized biventricular
ellipsoid mesh (26,642 points; raw bounding box x:[0.0, 4.0]
y:[-2.0, 3.0] z:[-2.0, 2.0]; boundary patches EPI, BASE, ENDO_LV,
ENDO_RV) and the same supplied fibre/sheet/transmural anatomy fields.

- `mesh/` — the one canonical copy of the shared mesh (already scaled
  to adult-heart-scale SI metres), anatomy fields (`fiber`, `sheet`,
  `sheetNormal`, `tm`, `tv`, `apicobasal`, `Conductivity`), and a
  generated Purkinje conduction network (`constant/purkinjeGraph`,
  grown by `generatePurkinjeTree` from UVC fields derived from
  `tm`/`tv`/`apicobasal` — see `electrophysiologyHeart/README.md`). Not a
  runnable case on its own
  — each case's `Allrun` copies in whatever subset it needs (renaming
  fields to its own convention where required) rather than checking in
  its own copy, to avoid tracking the same ~10MB mesh twice in git.
- `electrophysiologyHeart/` — monodomain electrophysiology on the
  idealized mesh, including conduction-system configuration for a
  Purkinje network.
- `electromechanicalHeart/` — sequential electromechanical coupling
  (TNNP ionic model, Land-Niederer active tension) on the same mesh.
- `conductionBlock/` — disease-variant sibling of `electrophysiologyHeart`:
  left or right bundle branch block via a structurally modified
  `purkinjeGraph` (one bundle severed, no ionic overrides). Two variants
  of one case, `./Allrun lbbb`/`./Allrun rbbb` — no default.
- `ionicPathology/` — another disease-variant sibling: acute ischemia or
  Brugada syndrome type-1 via `TNNP` `ionicConstantOverrides` (graph
  unmodified). Two variants of one case, `./Allrun ischemia`/
  `./Allrun brugada` — no default.

Scar/ablation is deferred — no scar region is defined for this anatomy
yet, unlike the reference cases this was migrated from
(`tutorials/PATHOS/`).
