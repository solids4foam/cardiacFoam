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

## Regression coverage

`conductionBlock`'s `lbbb` variant is regression-covered (see
`conductionBlock/regression/`) — it checks the structural technique
actually blocks propagation: an activation-time probe at an LV
Purkinje-myocardial junction site must still be un-activated at t = 35 ms. `ionicPathology`
has no automated regression yet — its pathology signatures only show up
over timescales (its own case runs 700ms) a short automated check
wouldn't meaningfully capture.
