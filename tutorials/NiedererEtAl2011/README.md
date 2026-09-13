# NiedererEtAl2011

The Niederer et al. (2011) slab benchmark and variants built on it: a shared
tissue-scale verification geometry, used first as a standalone monodomain
benchmark and then extended one capability at a time (Purkinje coupling,
electromechanics).

```text
NiedererEtAl2011/
├── NiedererEtAl2011verification/    the benchmark itself: monodomain slab
├── purkinjeNiedererEtAl2011/        + a 1D Purkinje network coupled in
└── electroMechanicalNiedererEtAl2011/  + sequential electro-mechanical coupling
```

## `NiedererEtAl2011verification/`

The Niederer slab verification workflow: `monodomainSolver` on the
benchmark's tissue-scale geometry, checked against the published activation
times. The baseline every variant below extends.

## `purkinjeNiedererEtAl2011/`

The same slab setup, with a small 1D Purkinje graph coupled into the 3D
monodomain tissue via `monodomain1DSolver` and PVJ coupling.

## `electroMechanicalNiedererEtAl2011/`

The same slab setup, sequentially coupled to solid mechanics via
`electroMechanicalModel` — electrophysiology drives an active-tension model,
which drives a solid deformation solve. Requires the `with-solids4foam`
build (`solids4foam` is not needed for the other two variants).

## Regression

All three variants are wired into `tutorials/Alltest-regression`; the
electromechanical variant is the one expected skip in `lightweight` build
mode. See `../README.md` for the full canonical table.
