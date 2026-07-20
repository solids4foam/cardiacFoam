# Pre-implementation outer-corrector diagnosis

Date: 2026-07-16

Status: historical diagnostic recorded before the equation-level
`correctNonOrthogonal()` loop was implemented. Re-run the harness for current
post-fix results; the no-op findings below intentionally document the original
defect.

Configuration:

- current `staggeredElectrophysicsAdvanceScheme`;
- tetrahedral monodomain MMS, least-squares gradients, `Gauss linear corrected` Laplacian;
- `N = 10, 20, 40`;
- common `endTime = 0.02`;
- resolution-dependent `deltaT ~ h^2` from the spatial verification ladder;
- serial execution, one timing repeat;
- outer-sweep reference: `nOuterCorrectors = 8`;
- predeclared acceptance: field difference/reference MMS L2 error `<= 0.01`, with no more than 1% MMS-error change at the next sweep count.

Raw logs, final fields, dictionaries, and `raw_results.csv` are generated under
`results/` by `run_corrector_study.sh`. That directory is ignored because it
contains generated fields.

## Code-path comparison

| Path | Reaction/ODE updates per physical step | Conduction/PVJ updates | Tissue diffusion solves | Valid interpretation |
|---|---:|---:|---:|---|
| Removed `pimpleStaggeredElectrophysicsAdvanceScheme` | repeated inside the outer loop | repeated, with source accumulation | nested use of the same `pimpleControl` | invalid; removed in `02d5ca87` |
| Current staggered, outer = 1 | 1 | 1 | 1 | single-pass staggered time step |
| Current staggered, outer = `k` | 1 | 1 | `k` | fixed-point reassembly of the tissue diffusion equation only |
| Current staggered, non-orthogonal = `k` | 1 | 1 | controlled only by outer count | dictionary entry is inactive without `correctNonOrthogonal()` |

The removed scheme is not an admissible numerical comparator: it advanced
stateful ODE/graph systems multiple times over the same physical `dt`, reused
one PIMPLE counter in nested loops, and accumulated PVJ deposits. Reintroducing
it to produce a comparison curve would compare the current method with a known
incorrect time integrator. The defensible comparison is therefore among outer
reassembly counts within the current single-pass staggered algorithm.

## Outer-corrector convergence

The table reports the L2 difference from the eight-sweep field divided by the
eight-sweep MMS L2 error.

| N | outer = 1 | outer = 2 | outer = 3 | outer = 4 | first passing count |
|---:|---:|---:|---:|---:|---:|
| 10 | 9.634% | 1.640% | 0.369% | 0.094% | 3 |
| 20 | 4.545% | 0.587% | 0.113% | 0.027% | 2 |
| 40 | 1.145% | 0.115% | 0.021% | 0.006% | 2 |

The correction iteration contracts consistently as the outer count rises.
Two sweeps satisfy the predeclared criterion on `N=20` and `N=40`, but miss it
on the coarsest mesh. Three sweeps are the smallest count satisfying the rule
for every tested resolution.

The raw MMS error is not a valid count-selection criterion. On all three
meshes, one sweep happens to give a slightly smaller MMS error than the
eight-sweep field. This is accidental cancellation between incomplete
cross-diffusion iteration and spatial truncation error. Corrector convergence
must instead be assessed from the field change relative to a converged outer
reference.

## Original `nNonOrthogonalCorrectors` no-op

At fixed `nOuterCorrectors = 2`, setting
`nNonOrthogonalCorrectors = 0, 1, 2` produced:

- identical MMS errors;
- zero L2 and Linf field differences;
- identical SHA-256 hashes of the final serial ASCII `Vm` fields;
- the same number of `Vm` linear solves per physical timestep.

This was direct numerical confirmation that `nNonOrthogonalCorrectors` was a
no-op in the pre-fix solver path. The current implementation now calls
`correctNonOrthogonal()`; these hashes must not be interpreted as post-fix
behaviour.

## Preliminary cost

At `N=40`, one unreplicated wall-time sample gave:

| Outer sweeps | Wall time (s) | Relative to two sweeps |
|---:|---:|---:|
| 1 | 43.09 | 0.96 |
| 2 | 44.98 | 1.00 |
| 3 | 47.08 | 1.05 |
| 4 | 48.34 | 1.07 |
| 8 | 54.20 | 1.20 |

These timings are indicative only; the paper run must use at least three
repeats and report a median and spread. The preliminary data suggest that
moving from two to three sweeps costs roughly 5% for this case, while eight
sweeps are unnecessary as a production setting.

## Post-fix decision

The original outer loop was not performing strong interdomain coupling; it was
reassembling the single corrected tissue equation. After implementing the
equation-level non-orthogonal loop, comparisons at `N=10,20` established that
two outer sweeps with zero additional non-orthogonal assemblies and one outer
sweep with one additional non-orthogonal assembly give bitwise-identical final
fields. The production manufactured cases were therefore migrated to the
second configuration. This preserves the reported numerical result while
assigning the repeated solve to the control that describes its purpose.

Additional assemblies converge the deferred-correction fixed point, but they
do not improve the discretisation order: in the short `N=10,20` comparison the
two-level order changed only from 2.184 with one solve to 2.190 with four
solves, while the absolute MMS error increased slightly as accidental error
cancellation disappeared.
