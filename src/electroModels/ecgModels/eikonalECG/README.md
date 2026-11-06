# eikonalECG

A surrogate ECG for eikonal runs. Instead of solving for the voltage, it shifts a tabulated action potential by the local activation time and sums the result over the cells with the same lead field as `pseudoECG`. It is not a general-purpose ECG model.

## What's available

- `eikonalECG`: the surrogate ECG model.
- Action-potential templates for endocardial, mid-myocardial and epicardial cells (`tissueTemplates.H`), and templates generated from a chosen ionic model for each tissue region.
- A manufactured-template mode, used to verify the ECG calculation: one action-potential-like event by default, or `template sine;` (`manufacturedTemplate.H`).

It writes `postProcessing/eikonalECG.dat`, with one column per electrode, sampled on its own output time grid rather than at every solver step.

## Several templates

With region weights `w_k` (from `ionicHeterogeneity`), the reconstructed voltage is `Vm = sum_k w_k U_k(t - tau)`. Its gradient has two parts, and both are used:

- `-sum_k w_k U_k'(t - tau) grad(tau)`, from the arrival of activation;
- `sum_k U_k(t - tau) grad(w_k)`, from the change of template between neighbouring cells. This carries the voltage difference between templates during repolarisation, for example endocardial against epicardial cells.

The compiled templates are selected by region name (`endocardialCells`, `mCells`, `epicardialCells`), so a subset such as endocardial and epicardial only is valid; other names need `personalizedTemplates`. `grad(w_k)` is computed once with zero-gradient walls. With `transitionMode hard`, it is concentrated on the region boundaries. With one template, or uniform weights, the second part is zero.

`checkGradVm true;` (default `false`) also evaluates the ECG from `fvc::grad` of the reconstructed `Vm`. It writes `postProcessing/eikonalECG_checkGradVm.dat` and prints the largest difference from the main ECG. The two agree as the mesh resolves the template upstroke and the weight transitions; the check is a diagnostic, not a reference.

A worked set-up is in [tutorials/manufacturedSolutions/eikonalECG](../../../../tutorials/manufacturedSolutions/eikonalECG/README.md).
