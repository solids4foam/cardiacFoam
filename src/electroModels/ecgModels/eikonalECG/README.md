# eikonalECG

A surrogate ECG for eikonal runs. Instead of solving for the voltage, it shifts a tabulated action potential by the local activation time and sums the result over the cells with the same lead field as `pseudoECG`. It is not a general-purpose ECG model.

## What's available

- `eikonalECG`: the surrogate ECG model.
- Action-potential templates for endocardial, mid-myocardial and epicardial cells (`tissueTemplates.H`), and templates generated from a chosen ionic model for each tissue region.
- A manufactured-template mode, used to verify the ECG calculation.

It writes `postProcessing/eikonalECG.dat`, with one column per electrode, sampled on its own output time grid rather than at every solver step.

Worked set-ups are in [tutorials/electrophysiologyProtocols/eikonalECGPersonalized](../../../../tutorials/electrophysiologyProtocols/eikonalECGPersonalized/README.md) and [tutorials/manufacturedSolutions/eikonalECG](../../../../tutorials/manufacturedSolutions/eikonalECG/README.md).
