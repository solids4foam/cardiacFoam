# Pathology techniques on the idealized biventricular heart

Two tutorials, each demonstrating a different *technique* for injecting a
specific pathology into the model. Both read the mesh, anatomy fields and
healthy Purkinje graph from `../mesh/`, and keep the healthy graph and ionic
model unless the technique itself changes them.

## `conductionBlock` — structural technique

Modifies the **conduction network itself**: `Allrun` derives
`purkinjeGraph.lbbb`/`.rbbb` from the healthy graph by zeroing a single
bridge edge's conductance, fully disconnecting one ventricle's Purkinje
subtree from the root (it's a tree — no alternate path exists). No ionic
model changes at all; the conduction delay is entirely structural. Use
this technique for anything that blocks or reroutes a conduction
pathway — bundle branch blocks here, but the same edge-zeroing approach
generalizes to any structural lesion of the Purkinje tree.

## `ionicPathology` — ionic-override technique

Modifies **cell-level ionic constants** via `ionicConstantOverrides` on
the TNNP model, leaving the conduction network untouched: `g_Na`/`g_CaL`/
`g_K1`/`P_NaK` scaling (global or region-restricted to a tissue band) for
acute ischemia and Brugada syndrome. Use this technique for anything
that's a channelopathy or a metabolic/ionic derangement rather than a
structural conduction problem.

## Regression

- `conductionBlock`: tier `full`. `lbbb` and `rbbb` run to 35 ms; the
  activation-time probe at the severed ventricle's Purkinje-myocardial
  junction site must still read -1 while the intact ventricle's probe has
  activated.
- `ionicPathology`: tier `none`. Its solver path is that of `../electroHeart`
  (`full`), and `ionicConstantOverrides` is checked by
  `../../electrophysiologyProtocols/ionicHeterogeneity` (`standard`).
