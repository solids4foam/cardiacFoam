# interpolateFibreField

Turns a cell-centred fibre field `f0` into the unit-length face field `f0f` that the mechanics solver needs.

## Every mechanics run needs `f0f`

`electroMechanicalLaw` reads `f0f` with `MUST_READ`, with no fallback. If `f0f` is missing, the run stops with a fatal error when the mechanical law is constructed.

How you get it depends on where your `f0` came from:

| Your fibres came from | Do this |
|---|---|
| `setFibreField`, which computes rule-based fibres on your mesh | Nothing extra: `setFibreField` writes `f0f` alongside `f0`. |
| Anywhere else: an anatomical dataset, another tool, or a different mesh | Run `interpolateFibreField` once `f0` is in place. |

> **The shipped `tutorials/idealizedHeart` case already includes a precomputed `f0f`** in `mesh/0/`. That is why `electroMechHeart/Allrun` copies it instead of calling this tool. That `f0f` belongs to that mesh: on any other mesh, generate your own. Never copy an `f0f` across meshes.

## Usage

```bash
interpolateFibreField                 # single-region case
interpolateFibreField -region solid   # multi-region case: run it on the solid region
```

- **Input:** `f0` (`volVectorField`) in the start-time directory, usually `0/`. It is required.
- **Output:** `f0f` (`surfaceVectorField`, unit length) in the same directory.

Run it again whenever `f0` or the mesh changes.

## How it differs from `setFibreField`

| | `setFibreField` | `interpolateFibreField` |
|---|---|---|
| Computes `f0` | yes, from a Laplace-solved transmural coordinate and a helix-angle rule | no, it reads the `f0` you provide |
| Writes `f0f` | yes | yes |
| Use it when | you want rule-based fibres | your fibres come from data |

## Source

[`interpolateFibreField.C`](interpolateFibreField.C) builds in both full and lightweight mode, since it links only `finiteVolume` and `meshTools`. Its output is only used in full (solids4foam) mode.
