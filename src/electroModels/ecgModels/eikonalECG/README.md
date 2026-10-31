# eikonalECG

A surrogate ECG for eikonal runs. Instead of solving for the voltage, it shifts a tabulated action potential by the local activation time and sums the result over the cells with the same lead field as `pseudoECG`. It is not a general-purpose ECG model.

## What's available

- `eikonalECG`: the surrogate ECG model.
- Action-potential templates for endocardial, mid-myocardial and epicardial cells (`tissueTemplates.H`), and templates generated from a chosen ionic model for each tissue region.
- A manufactured-template mode, used to verify the ECG calculation: one action-potential-like event by default, or `template sine;` (`manufacturedTemplate.H`).

It writes `postProcessing/eikonalECG.dat`, with one column per electrode, sampled on its own output time grid rather than at every solver step.

## Templates from converged states

`personalizedTemplates` makes one template per `ionicHeterogeneity` region by pacing a single cell of the chosen model.
Each template records Vm and Cai, and its time zero is the steepest upstroke of the last beat, so a cell activated at time τ
depolarises at τ. By default each cell starts from the model's initial state and is paced for `nBeats`.

`regionStates` names a file in the `setIonicRestartState` map format (`regionStates { <region> "<Model>State file"; }`).
Each region's cell then starts from its converged single-cell state, and `nBeats` defaults to 1. The same file seeds a
monodomain run with `setIonicRestartState`, so both use the same cells.

```cpp
personalizedTemplates
{
    ionicModelConfig { ionicModel Gaur; tissue myocyte; solver RKF45; singleCellStimulus { ... } }
    regionStates  "system/setIonicRestartStateDict";
    duration      0.5;
    dt            1e-4;
}
```

A worked set-up is in [tutorials/manufacturedSolutions/eikonalECG](../../../../tutorials/manufacturedSolutions/eikonalECG/README.md).
