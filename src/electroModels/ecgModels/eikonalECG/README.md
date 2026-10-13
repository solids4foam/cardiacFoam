# eikonalECG

A template-based pseudo-ECG observer for eikonal runs. It reconstructs the
voltage gradient from an action-potential template and activation time, then
applies the analytical recording kernel used by `pseudoECG`. The
`eikonalTemplateSource` helper now supplies this reconstruction to both
eikonal ECG observers.

## What's available

- `eikonalECG`: the surrogate ECG model.
- Action-potential templates for endocardial, mid-myocardial and epicardial cells (`tissueTemplates.H`), and templates generated from a chosen ionic model for each tissue region.
- A manufactured-template mode, used to verify the ECG calculation.

It writes `postProcessing/eikonalECG.dat`, with one column per electrode, sampled on its own output time grid rather than at every solver step.

A worked set-up is in [tutorials/manufacturedSolutions/eikonalECG](../../../../tutorials/manufacturedSolutions/eikonalECG/README.md).

`templateType endo`, `mid` or `epi` selects one compiled template across the
heart; `manufactured` selects the sine fixture. The default `regional` path
preserves the existing regional/personalized template behavior. The
[eikonalLeadFieldECG observer](../eikonalLeadFieldECG/README.md) uses the same
uniform-template reconstruction with imported reciprocal sensitivities.
