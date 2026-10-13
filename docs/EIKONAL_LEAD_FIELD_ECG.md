# Eikonal lead-field ECG

## Mathematical model

A lead sensitivity `L_l(x)` is a vector at each myocardial location for one
ECG lead. In production it is built from a reciprocal conductor potential
`z_l` as `L_l = -Gi^T grad(z_l)`, with the normalized test current incorporated
into the sensitivity. It has units m⁻². It includes the physical intracellular
conductivity and electrode reference; the online observer must not multiply
by Gi again. Torso geometry and passive conductor properties enter the
preprocessing problem, not the eikonal propagation solve.

For a spatially uniform action-potential template,

```text
Vm(x,t) = U(t - tau(x)),
grad(Vm) = -U'(t - tau) grad(tau),
y_l(t) = int_heart L_l . grad(Vm) dx
       ~= sum_cells Vcell L_l,cell . [-U'(t - tau_cell) grad_h(tau)_cell].
```

`tau` is activation time in seconds. `U` is voltage in volts, and `U'` is in
V/s. The product with `grad(tau)` gives V/m, and the volume integral gives
volts. Eikonal propagation supplies activation time rather than a solved
voltage trajectory. The observer therefore reconstructs a voltage gradient
from the template; it does not use the eikonal domain's placeholder `Vm`.

The new observer is passive and samples the already-computed steady
activation field on its own time grid. Changing its sampling interval is
not a propagation time-step refinement.

## Classes and integration

| Class/file | Responsibility |
|---|---|
| [eikonalTemplateSource](../src/electroModels/ecgModels/eikonalECG/eikonalTemplateSource.H) | Shared template configuration, generation, regional weights and chain-rule reconstruction, extracted from the existing eikonal ECG. |
| [eikonalECG](../src/electroModels/ecgModels/eikonalECG/eikonalECG.C) | Retains the analytical pseudo-ECG recording kernel and existing sampling interface; delegates template reconstruction to the shared helper. Only the master writes its series. |
| [leadFieldData](../src/electroModels/ecgModels/leadFieldECG/leadFieldData.H) | Adds `projectGradient`: validates a fixed matching mesh and V/m dimensions, then applies the existing imported sensitivities and MPI volume sum. Monodomain/bidomain `project` computes its private voltage gradient and delegates to this method. |
| [eikonalLeadFieldECG](../src/electroModels/ecgModels/eikonalLeadFieldECG/eikonalLeadFieldECG.H) | New runtime-selected observer combining uniform template reconstruction, private activation-gradient construction, imported sensitivities and its own sampling/output loop. |
| [ecgDomain](../src/electroModels/electroDomains/ecgDomain/ecgDomain.C) | Routes the new observer only to `eikonalSolver`; channel names come from the sensitivities, without artificial electrode positions. |
| [manufacturedEikonalLeadFieldECGVerifier](../src/verificationModels/ecgVerification/manufacturedEikonalLeadFieldECGVerifier.H) | Independent surface-voltage oracle for the prescribed constant-sensitivity cube. Requires activation time, without requiring Vm. |
| [setManufacturedLeadFields](../applications/utilities/setManufacturedLeadFields/setManufacturedLeadFields.C) | Adds synthetic eikonal sensitivities, four affine template/unit/state controls and deterministic invalid-input checks. |

The shared helper preserves the existing regional/personalized pseudo-ECG
path. The existing manufactured pseudo-ECG configuration still selects the
same sine template implicitly. Explicit `templateType endo`, `mid`, `epi`
or `manufactured` selects one template across the heart; both observers can
use these explicit selections.

The compiled tables store voltage in mV and time in seconds. Their derivative
is the secant slope of the existing piecewise-linear interpolant, multiplied
by 1e-3 for V/s. Outside the table range the derivative is zero, consistent
with clamped endpoint voltages.

The reciprocal observer creates `eikonalLeadFieldObservationTau`, copies the
activation owner's published wall policy, and selects its gradient scheme
through `eikonalLeadFieldGradTau`. Its gradient identity is separate from the
source's `grad(activationTime)` cache. It neither corrects the source field
nor overwrites that cache. The existing sensitivity geometry, connectivity,
conductivity and provenance checks remain in the importer.

## Configuration

After preparing sensitivities with the existing conductor workflow, add:

```foam
ecgDomains
{
    Reciprocal
    {
        ecgSolver eikonalLeadFieldECG;
        templateType endo;
        importMode originalConductor;
        sampling { start 0; end 1; deltaT 0.001; }
    }
}
```

Use `importMode heartOnly` only for an explicitly exported myocardial mesh
and its matching sensitivities, as for the voltage-based observer. Physical
Gi and Ge remain mandatory in `leadFieldProperties`, through the established
uniform or field input contract. The eikonal propagation tensor does not
uniquely specify that pair. No conductor solve or new potential-domain
subclass was added for this extension.

Output is `postProcessing/eikonalLeadFieldECG.dat`, with seconds and volts.
It is generated once per observer instance; repeated calls return its final
values without repeating the sampling loop. This follows the batch eikonal
observation lifecycle, rather than the voltage-based restart append schedule.
Use one observer of each type per run with the current fixed output filenames.
For simultaneous comparison, give both observers the same uniform template
and sampling block.

## Scope of heterogeneous templates

For `Vm = sum_r w_r(x) U_r(t-tau(x))`, the full spatial gradient is

```text
grad(Vm) = -sum_r w_r U'_r grad(tau) + sum_r U_r grad(w_r).
```

The legacy regional path holds the local weights fixed during differentiation
and computes the first term. Its behavior is preserved. The new reciprocal
observer rejects regional/personalized templates instead of silently
identifying that approximation with the full blended-voltage gradient.
Smooth blending and discontinuous regional interfaces need a separately
specified and verified source model before extending this observer.

## Verification and tutorials

The [existing eikonal cube](../tutorials/manufacturedSolutions/eikonalECG/README.md)
adds `-leadFieldECG yes`; its default and original numerical references are
preserved. The [dual variant](../tutorials/manufacturedSolutions/eikonalLeadFieldECG/README.md)
runs that same maintained case in `localCase`, with 40 cells per direction.
The regression driver copies the shared base when needed, and the tutorial
index includes the new variant.

The fixture utility first uses affine activation time and exact boundary
values to isolate the observation code. It exercises manufactured, endo,
mid and epi templates through the production observer. Its independent
oracle scans table intervals directly, so it checks derivative selection,
mV-to-V conversion, sign and volume weighting without calling the production
interpolator. A provider without Vm, a deliberately stale gradient, unchanged
source/boundary/cache arrays and repeated solve calls check the observer
contract. Separate failures exercise regional/unknown template rejection and
invalid sampling.

The solved-field test uses `tau_exact = exp(x+y+z)` and `U = sin(2*pi*s)`
volts. Prescribed sensitivities are `L_positive=(1,2,3)` and its negative,
in m⁻². For constant L, the divergence theorem converts the reference to face
voltage differences. Cube symmetry gives a two-dimensional face integral,
independent of the production gradient reconstruction and volume quadrature.
The C++ reference uses Gauss-Legendre orders 64/96; the Python control uses
order 80 and checks signs, finite monotone traces, common sampling, written
references and independently normalized recording errors. It also bounds
activation-time, solved-gradient and exact-field reconstruction errors at
the default resolution.

These synthetic sensitivities verify the new source-to-observer path. The
previous conductor recovery controls verify construction separately. The new
cube does not establish anatomical torso accuracy or a discrete-adjoint
identity for the complete bath operator.

## Executed evidence

The first native 40³/six-rank dual run records 201 samples over [0,2] s:

| Quantity | Result |
|---|---|
| Reciprocal lead cohort normalized RMS error | 0.01679165366 (1.68%) |
| Reciprocal lead cohort RMS error | 0.001649418797 V |
| Reciprocal maximum absolute error | 0.002330785578 V |
| Raw pseudo-ECG E3–E5 cohort normalized RMS error | 0.02498087834 (2.50%) |
| Surface reference Q64/Q96 maximum difference | 1.94e-15 V |
| Independent Q80 reference maximum mismatch with written Q96 | 4.82e-10 V, within scientific precision-8 output rounding |

The two cohorts measure different functionals and cannot rank the recording
models. A 20³ midpoint calculation of the exact constant-sensitivity
integrand is under-resolved because of rapid phase variation near the upper
corner. The 40³ result is an accuracy check at a selected resolution, not a
spatial convergence-order claim.

Paired native serial 40³ runs with and without the reciprocal observer write
byte-identical activation-time fields and pseudo-ECG traces. The validation
report records the final suite/platform checks; local logs and paired-run
snapshots are retained outside the repository under
`/tmp/cardiacFoam-lead-field-current-numerics/eikonal`.
