# Pathology techniques on the idealized biventricular heart

Two sibling tutorials, each demonstrating a different *technique* for
injecting a specific pathology into the model — not a catalogue of
diseases, a reference for **how** to build one. Both share the same mesh
and anatomy as `../electroHeart`/`../electroMechHeart` (`../mesh/`) and
default to the same healthy Purkinje graph and ionic model unless the
technique itself says otherwise.

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
actually blocks propagation, using the same injection-probe pattern
`../electroHeart` uses to check propagation *does* happen. `ionicPathology`
has no automated regression yet — its pathology signatures only show up
over timescales (its own case runs 700ms) a short automated check
wouldn't meaningfully capture.
