# NiedererEtAl2011

The Niederer et al. (2011) slab benchmark and variants built on it: a shared
tissue-scale verification geometry, used first as a standalone monodomain
benchmark and then extended to electromechanics.

```text
NiedererEtAl2011/
├── NiedererEtAl2011verification/    the benchmark itself: monodomain slab
└── electroMechanicalNiedererEtAl2011/  + sequential electro-mechanical coupling
```

Purkinje network coupling is exercised in
`../electrophysiologyProtocols/purkinjeRestitution2D`.

## `NiedererEtAl2011verification/`

The Niederer slab verification workflow: `monodomainSolver` on the
benchmark's tissue-scale geometry, checked against the published activation
times. The baseline the variant below extends.

## `electroMechanicalNiedererEtAl2011/`

The same slab setup, sequentially coupled to solid mechanics via
`electroMechanicalModel` — electrophysiology drives an active-tension model,
which drives a solid deformation solve. Requires the `with-solids4foam`
build (`solids4foam` is not needed for the verification case).

## Regression

Both cases are wired into `tutorials/Alltest-regression`; the
electromechanical variant is the one expected skip in `lightweight` build
mode. See `../README.md` for the full canonical table.
