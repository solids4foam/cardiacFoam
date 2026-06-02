# setFibreField

Computes and writes cardiac fibre, sheet-normal, and transmural-direction
fields for ellipsoidal-ventricle benchmarks following the Rossi-Lassila
rule-based approach.

The fibre helix angle varies linearly from `+60°` at the endocardium to `-60°`
at the epicardium (Land et al. 2015 benchmark defaults).

## What it does

1. Reads a transmural-distance field `t` from the initial time directory.
2. Solves a Laplace equation (`∇²t = 0`) to propagate the transmural gradient
   from the endocardium (`t = 0`) to the epicardium (`t = 1`).
3. Computes the orthonormal frame `(et, en, el)`:
   - `et` — transmural direction (`∇t / |∇t|`)
   - `en` — normal direction (projection of `k = (0,0,1)` onto the plane ⊥ to `et`)
   - `el` — longitudinal direction (`en × et`)
4. Rotates `el` about `et` by the linearly varying helix angle `α(t)` to
   produce the fibre direction `f0`.
5. Writes cell-centred (`f0`, `et`, `en`, `el`, `alphaRadians`) and
   face-centred (`f0f`, `etf`, `enf`, `elf`, `alphaRadiansf`) fields.

## Prerequisites

The initial time directory must contain a `volScalarField t` with:

- `fixedValue 0` on the endocardium patch
- `fixedValue 1` on the epicardium patch
- `zeroGradient` on all other patches

## Usage

```bash
setFibreField [-region <regionName>]
```

## Written fields

| Field | Type | Description |
|---|---|---|
| `t` | `volScalarField` | Solved transmural distance |
| `et` | `volVectorField` | Transmural direction (cell-centred) |
| `en` | `volVectorField` | Sheet-normal direction (cell-centred) |
| `el` | `volVectorField` | Longitudinal direction (cell-centred) |
| `alphaRadians` | `volScalarField` | Helix angle (radians, cell-centred) |
| `f0` | `volVectorField` | Fibre direction (cell-centred) |
| `etf` | `surfaceVectorField` | Transmural direction (face-centred) |
| `enf` | `surfaceVectorField` | Sheet-normal direction (face-centred) |
| `elf` | `surfaceVectorField` | Longitudinal direction (face-centred) |
| `alphaRadiansf` | `surfaceScalarField` | Helix angle (radians, face-centred) |
| `f0f` | `surfaceVectorField` | Fibre direction (face-centred) |

## Notes

- The endocardial and epicardial helix angles (`alphaEndo`, `alphaEpi`) are
  currently hard-coded as `+60°` and `−60°`. Parameterisation is planned.
- The base direction `k = (0, 0, 1)` follows the Rossi-Lassila convention
  for idealized ellipsoidal geometries.
