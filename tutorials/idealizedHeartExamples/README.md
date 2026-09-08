# Idealized heart examples

Both cases in this directory share the same idealized biventricular
ellipsoid mesh (26,642 points; raw bounding box x:[0.0, 4.0]
y:[-2.0, 3.0] z:[-2.0, 2.0]; boundary patches EPI, BASE, ENDO_LV,
ENDO_RV) and the same supplied fibre/sheet/transmural anatomy fields.

- `mesh/` — the one canonical copy of the shared mesh (already scaled
  to adult-heart-scale SI metres) and anatomy fields (`fiber`, `sheet`,
  `sheetNormal`, `tm`, `tv`, `apicobasal`, `Conductivity`). Not a
  runnable case on its own — each case's `Allrun` copies in whatever
  subset it needs (renaming fields to its own convention where
  required) rather than checking in its own copy, to avoid tracking
  the same ~10MB mesh twice in git.
- `electrophysiologyHeart/` — monodomain electrophysiology on the
  idealized mesh, including conduction-system configuration for a
  Purkinje network.
- `electromechanicalHeart/` — sequential electromechanical coupling
  (TNNP ionic model, Land-Niederer active tension) on the same mesh.

Future work may add disease-variant siblings of `electrophysiologyHeart`
built on this same anatomy (for example a shared conduction-block
variant covering both LBBB and RBBB, and a standalone ischemia case) —
neither exists yet.
