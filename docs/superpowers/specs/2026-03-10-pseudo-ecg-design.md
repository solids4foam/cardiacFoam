# Design: Pseudo-ECG in greensFunctionECGElectro

Date: 2026-03-10

## Goal

Extend `greensFunctionECGElectro` to simultaneously compute the pseudo-ECG (dipole)
potential alongside the existing Green's function (monopole) potential, writing both
to separate outputs in a single cell loop pass.

---

## Physics

### Green's function (existing)

```
phi_greens(P) = 1/(4*pi*sigmaT) * sum_c [ Is_c * V_c / |C_c - P| ]
```

where `Is_c = -div(Gi * grad(Vm))` is the current source density.

### Pseudo-ECG (new)

```
phi_pseudo(P) = sum_c [ (grad(Vm)_c . r_vec_c) * V_c / |C_c - P|^3 ]
```

where `r_vec_c = C_c - P`. No `1/(4*pi)` or `sigmaT` prefactor — purely dimensional.
A scaling constant can be added later if absolute calibration is needed.

---

## Approach

Extend `greensFunctionECGElectro` in-place. No new class. Both integrals share the
same double loop over heart cells × evaluation points, so the cost is one extra dot
product and accumulator per cell per point.

---

## Data changes

One new private member:

```cpp
//- Output stream for pseudoECG.dat (null if electrodes not configured)
autoPtr<OFstream> outputPseudoPtr_;
```

Initialised in the constructor alongside `outputPtr_` when the `electrodes` subdict
is present. Header written immediately with the same electrode names.

No new `volField` members — `grad(Vm)` is computed locally in `evolve()` and not stored.

---

## `evolve()` changes

Before the cell loop, compute the gradient field once:

```cpp
const vectorField& gradVm =
    fvc::grad(Vm()).ref().primitiveField();
```

Run a single combined loop with two accumulators:

```cpp
List<scalar> localSumsGreens(nAll, scalar(0));
List<scalar> localSumsPseudo(nAll, scalar(0));

forAll(Ctrs, cI)
{
    const scalar IsV     = IsI[cI]*Vols[cI];
    const vector gradVmV = gradVm[cI]*Vols[cI];

    for (label pI = 0; pI < nAll; pI++)
    {
        const vector r_vec = Ctrs[cI] - allPoints[pI];
        const scalar r     = mag(r_vec);
        if (r > VSMALL)
        {
            localSumsGreens[pI] += IsV / r;
            localSumsPseudo[pI] += (gradVmV & r_vec) / (r*r*r);
        }
    }
}
```

After `Foam::reduce`, apply `invCoeff` to `localSumsGreens` as before.
`localSumsPseudo` is used as-is.

---

## Output changes

### Electrode files

| File | Contents | When |
|------|----------|------|
| `postProcessing/ECG.dat` | Green's function values (unchanged) | every timestep |
| `postProcessing/pseudoECG.dat` | Pseudo-ECG values | every timestep |

Both files share the same header format: `# time  LA  RA  LL ...`

### Torso surface VTK

`writeTorsoVtk()` signature changes to accept both fields:

```cpp
void writeTorsoVtk
(
    const scalarList& phiGreens,
    const scalarList& phiPseudo
) const;
```

The `CELL_DATA` block in the VTK file gains a second `SCALARS` section:

```
CELL_DATA N
SCALARS phiGreens float 1
LOOKUP_TABLE default
...
SCALARS phiPseudo float 1
LOOKUP_TABLE default
...
```

Written at `outputTime()` only, same gate as existing torso VTK logic.

---

## Files to modify

| File | Change |
|------|--------|
| `src/electroModels/greensFunctionECGElectro/greensFunctionECGElectro.H` | Add `outputPseudoPtr_`; update `writeTorsoVtk()` signature |
| `src/electroModels/greensFunctionECGElectro/greensFunctionECGElectro.C` | Constructor: open `pseudoECG.dat` + write header; `evolve()`: combined loop + pseudo output; `writeTorsoVtk()`: write two scalar fields |

No other files need to change.

---

## Out of scope

- `1/(4*pi*sigmaT)` prefactor for pseudo-ECG (easy to add if needed)
- Separate torso VTK files per method
- New electroModel class or base class refactor
