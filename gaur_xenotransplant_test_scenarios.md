# Gaur Porcine Ionic Model: Xenotransplantation Test Scenarios

## Purpose

This document defines a small set of test scenarios for demonstrating the capabilities of a porcine ventricular ionic model, based on the Gaur model, in the context of cardiac xenotransplantation.

The focus is on **calcium-driven contractility loss**, since this has been identified by the doctors as a major phenotype of interest.

The goal is not to claim that these scenarios are biologically validated. Instead, they are **synthetic, interpretable perturbations** designed to show that the model can represent different possible mechanisms of graft dysfunction across scales.

---

## Scales to Demonstrate

```text
Cell scale
    -> action potential, calcium transient, SR calcium, ionic currents

Tissue scale
    -> activation spread, conduction velocity, regional calcium loss

Organ scale
    -> pseudo-ECG, activation sequence, global calcium/contraction proxy

Xenotransplant biology scale
    -> rejection-like, ischemia-like, edema-like, and contractility-loss phenotypes
```

---

## Core Variables to Expose

The following scaling factors should be added to the model as simple multipliers.

```cpp
struct ScenarioScales
{
    scalar GNa   = 1.0;
    scalar GCaL  = 1.0;
    scalar GKr   = 1.0;
    scalar GKs   = 1.0;
    scalar GK1   = 1.0;

    scalar INaK  = 1.0;
    scalar INaCa = 1.0;

    scalar Jrel  = 1.0;
    scalar Jup   = 1.0;
    scalar Jleak = 1.0;
    scalar Jtr   = 1.0;
    scalar Jdiff = 1.0;

    scalar D     = 1.0;

    scalar contractility = 1.0;
};
```

### Most Important Calcium / Contractility Variables

```text
GCaL_scale          -> L-type calcium entry
Jrel_scale          -> SR calcium release
Jup_scale           -> SERCA calcium uptake
Jleak_scale         -> SR calcium leak
INaCa_scale         -> sodium-calcium exchanger
contractility_scale -> calcium-to-tension or active-tension scaling
```

---

# Scenario 0: Baseline Porcine Graft

## Purpose

Reference case.

## Scaling Factors

```cpp
ScenarioScales baseline;

baseline.GNa   = 1.00;
baseline.GCaL  = 1.00;
baseline.GKr   = 1.00;
baseline.GKs   = 1.00;
baseline.GK1   = 1.00;

baseline.INaK  = 1.00;
baseline.INaCa = 1.00;

baseline.Jrel  = 1.00;
baseline.Jup   = 1.00;
baseline.Jleak = 1.00;
baseline.Jtr   = 1.00;
baseline.Jdiff = 1.00;

baseline.D     = 1.00;
baseline.contractility = 1.00;
```

## Expected Phenotype

```text
Normal action potential
Normal calcium transient
Normal contraction proxy
```

---

# Scenario 1: Primary Calcium-Release Failure

## Purpose

Model loss of contractility due primarily to impaired SR calcium release.

This is one of the most clinically relevant first tests if contractility loss is the main observed problem.

## Scaling Factors

```cpp
ScenarioScales caReleaseFailure;

caReleaseFailure.GCaL  = 0.95;
caReleaseFailure.Jrel  = 0.60;
caReleaseFailure.Jup   = 0.85;
caReleaseFailure.Jleak = 1.10;

caReleaseFailure.contractility = 1.00;
```

## Expected Phenotype

```text
Calcium transient amplitude decreases strongly
Action potential may change only modestly
Contraction proxy decreases because cytosolic calcium peak is smaller
```

## Interpretation

```text
The graft cell is still electrically excitable, but SR calcium release is impaired.
```

---

# Scenario 2: SERCA / Relaxation Dysfunction

## Purpose

Model impaired calcium reuptake and slower relaxation.

This is useful if the clinicians observe poor relaxation, diastolic dysfunction, or progressive contractile weakening.

## Scaling Factors

```cpp
ScenarioScales sercaFailure;

sercaFailure.GCaL  = 1.00;
sercaFailure.Jrel  = 0.90;
sercaFailure.Jup   = 0.50;
sercaFailure.Jleak = 1.20;

sercaFailure.INaCa = 1.10;

sercaFailure.contractility = 1.00;
```

## Expected Phenotype

```text
Slower calcium transient decay
Higher diastolic calcium
Reduced SR refilling
Possible beat-to-beat instability
Weaker subsequent contractions
```

---

# Scenario 3: Calcium-Entry Reduction

## Purpose

Test whether loss of L-type calcium entry alone can explain contractility loss.

This checks whether the contractility problem could begin at the membrane level rather than inside the SR.

## Scaling Factors

```cpp
ScenarioScales caEntryFailure;

caEntryFailure.GCaL = 0.60;
caEntryFailure.Jrel = 0.85;

caEntryFailure.Jup   = 1.00;
caEntryFailure.Jleak = 1.00;

caEntryFailure.contractility = 1.00;
```

## Expected Phenotype

```text
Reduced plateau calcium entry
Reduced calcium-induced calcium release
Smaller calcium transient
Possibly shorter plateau or altered APD
```

---

# Scenario 4: Myofilament Contractility Loss With Preserved Calcium

## Purpose

Separate electrical/calcium failure from mechanical coupling failure.

This is important because contractility can decrease even if the calcium transient is preserved.

## Scaling Factors

```cpp
ScenarioScales myofilamentFailure;

myofilamentFailure.GCaL  = 1.00;
myofilamentFailure.Jrel  = 1.00;
myofilamentFailure.Jup   = 1.00;
myofilamentFailure.Jleak = 1.00;

myofilamentFailure.contractility = 0.50;
```

## Expected Phenotype

```text
Normal action potential
Normal calcium transient
Reduced active tension or contraction proxy
```

## Interpretation

```text
The cell has enough calcium, but the myofilaments respond less strongly.
```

---

# Scenario 5: Rejection-Like Calcium Injury

## Purpose

Model combined immune/metabolic injury centered on calcium handling.

This is likely the best xenotransplant-relevant cell-scale scenario.

## Scaling Factors

```cpp
ScenarioScales rejectionCaInjury;

rejectionCaInjury.GCaL  = 0.80;
rejectionCaInjury.Jrel  = 0.65;
rejectionCaInjury.Jup   = 0.60;
rejectionCaInjury.Jleak = 1.50;

rejectionCaInjury.INaK  = 0.80;
rejectionCaInjury.INaCa = 1.20;

rejectionCaInjury.GKr   = 0.85;
rejectionCaInjury.GK1   = 0.90;

rejectionCaInjury.contractility = 0.75;
```

## Expected Phenotype

```text
Smaller calcium transient
Slower calcium removal
Higher diastolic calcium
Mild APD prolongation may occur
Strong decrease in contraction proxy
Possible instability under pacing
```

## Biological Proxy

```text
Inflammation
ROS
Complement injury
Energetic stress
```

## Important Note

This is a synthetic rejection-like perturbation, not a validated rejection model.

---

# Scenario 6: Ischemia-Reperfusion / Metabolic Stress

## Purpose

Model ATP/pH-related acute dysfunction.

This is useful for early post-transplant or preservation/reperfusion effects.

## Scaling Factors

```cpp
ScenarioScales ischemiaReperfusion;

ischemiaReperfusion.GNa   = 0.85;
ischemiaReperfusion.GCaL  = 0.75;

ischemiaReperfusion.INaK  = 0.60;
ischemiaReperfusion.INaCa = 1.15;

ischemiaReperfusion.Jrel  = 0.80;
ischemiaReperfusion.Jup   = 0.55;
ischemiaReperfusion.Jleak = 1.40;

ischemiaReperfusion.GK1   = 0.85;

ischemiaReperfusion.contractility = 0.80;
```

## Expected Phenotype

```text
Reduced excitability
Smaller calcium transient
Slower calcium recovery
Possible diastolic calcium elevation
Weaker contraction
```

---

# Scenario 7: Repolarization Reserve Loss

## Purpose

Test APD/QT-type behavior without making calcium failure the main driver.

## Scaling Factors

```cpp
ScenarioScales repolarizationReserveLoss;

repolarizationReserveLoss.GKr   = 0.50;
repolarizationReserveLoss.GKs   = 0.75;
repolarizationReserveLoss.GCaL  = 1.10;

repolarizationReserveLoss.Jrel  = 1.00;
repolarizationReserveLoss.Jup   = 1.00;
repolarizationReserveLoss.Jleak = 1.00;

repolarizationReserveLoss.contractility = 1.00;
```

## Expected Phenotype

```text
APD90 increases
Plateau duration increases
Calcium transient may broaden
Possible EAD-like behavior if pushed further
```

---

# Scenario 8: Conduction-Slowing Tissue Phenotype

## Purpose

Show tissue-scale capability.

This is not the main contractility-loss hypothesis, but it demonstrates conduction and activation effects.

## Scaling Factors

```cpp
ScenarioScales conductionSlowing;

conductionSlowing.GNa = 0.75;
conductionSlowing.GK1 = 0.90;

conductionSlowing.D   = 0.60;

conductionSlowing.GCaL  = 1.00;
conductionSlowing.Jrel  = 1.00;
conductionSlowing.Jup   = 1.00;
conductionSlowing.Jleak = 1.00;

conductionSlowing.contractility = 1.00;
```

## Expected Phenotype

```text
Slower activation
Reduced conduction velocity
Wider activation time
Possible QRS-like widening in pseudo-ECG
Spatial delay in calcium transient timing
```

---

# Scenario 9: Regional Graft Injury

## Purpose

Show that the framework can model heterogeneous injury.

Apply these modifiers only in an injured tissue region. The rest of the tissue remains baseline.

## Scaling Factors

```cpp
ScenarioScales regionalInjury;

regionalInjury.GCaL  = 0.75;
regionalInjury.Jrel  = 0.55;
regionalInjury.Jup   = 0.60;
regionalInjury.Jleak = 1.50;

regionalInjury.D     = 0.75;

regionalInjury.contractility = 0.65;
```

## Expected Phenotype

```text
Regional calcium transient loss
Regional contraction loss
Delayed or weakened activation-contraction pattern
Possible mechanical dyssynchrony
```

## Suggested Maps

```text
Activation time
Calcium peak
Active tension peak
APD90
```

---

# Scenario 10: Severe Calcium-Collapse Stress Test

## Purpose

Show the failure boundary of the model.

This is a stress test and should not be presented as the default realistic case.

## Scaling Factors

```cpp
ScenarioScales severeCalciumCollapse;

severeCalciumCollapse.GCaL  = 0.60;
severeCalciumCollapse.Jrel  = 0.40;
severeCalciumCollapse.Jup   = 0.40;
severeCalciumCollapse.Jleak = 2.00;

severeCalciumCollapse.INaK  = 0.60;
severeCalciumCollapse.INaCa = 1.30;

severeCalciumCollapse.contractility = 0.50;
```

## Expected Phenotype

```text
Very small calcium transient
Possible high diastolic calcium
Poor recovery
Large contractility loss
Possible numerical or physiological instability
```

---

# Recommended Minimal Test Set

For the first demonstration, run only these scenarios:

```text
0. Baseline porcine graft
1. Primary calcium-release failure
2. SERCA / relaxation dysfunction
3. Myofilament contractility loss
4. Rejection-like calcium injury
```

If tissue or organ-scale capability is also needed, add:

```text
5. Regional graft injury
6. Conduction-slowing tissue phenotype
```

---

# Suggested Cell-Scale Outputs

Run a single Gaur cell paced at 1 Hz.

## Pacing Protocol

```text
BCL = 1000 ms
nBeats = 20 to 100 for a quick test
nBeats = 500 or more for better steady-state metrics
dt = solver default
stimulus = identical for all cases
```

## Save These Traces

```text
Vm(t)
Cai(t)
CaSR(t)
ICaL(t)
INaCa(t)
Jrel(t)
Jup(t)
Jleak(t)
activeTension(t), if available
```

## Compute These Metrics

```text
APD50
APD90
dV/dtmax
Ca transient peak
Ca transient amplitude
Diastolic calcium
Ca transient decay time
SR calcium load
Peak active tension
Time-to-peak tension
Relaxation time
```

---

# Suggested Tissue-Scale Outputs

Run a 2D sheet or 3D block.

## Best Scenarios

```text
Baseline
Conduction-slowing tissue phenotype
Regional graft injury
```

## Save These Maps

```text
Activation time
APD90
Conduction velocity
Calcium transient peak
Diastolic calcium
Active tension peak
```

## Main Demonstration

```text
The same Gaur cell model can be embedded in tissue to show spatial propagation and regional graft dysfunction.
```

---

# Suggested Organ-Scale Outputs

Run a simple ventricular geometry, if available.

## Best Scenarios

```text
Baseline
Global rejection-like calcium injury
Regional graft injury
Conduction-slowing tissue phenotype
```

## Save These Outputs

```text
Activation sequence
Pseudo-ECG
Global calcium peak
Global active tension or pressure proxy
Regional calcium peak maps
Regional tension maps
```

## Main Demonstration

```text
Cell-scale calcium impairment can propagate to organ-scale contraction loss.
```

---

# Clinically Focused Comparison

Given the doctors' focus on calcium and contractility loss, the most important comparison is:

```text
A. Baseline
B. Primary calcium-release failure
C. Myofilament contractility failure
D. Rejection-like calcium injury
```

## Expected Interpretation Table

| Scenario | Vm | Calcium Transient | Contractility |
|---|---:|---:|---:|
| Baseline | Normal | Normal | Normal |
| Primary calcium-release failure | Near-normal | Low | Low |
| Myofilament failure | Normal | Normal | Low |
| Rejection-like calcium injury | Altered | Low / slow | Very low |

This comparison is useful because it separates different possible causes of contractility loss:

```text
Less calcium available
vs
same calcium but weaker myofilament response
vs
combined electrical, calcium, and metabolic injury
```

---

# Suggested Figure Set

## Figure 1: Cell-Scale Traces

Overlay:

```text
Vm(t)
Cai(t)
activeTension(t)
```

for:

```text
Baseline
Primary calcium-release failure
SERCA / relaxation dysfunction
Myofilament failure
Rejection-like calcium injury
```

## Figure 2: Metrics Table or Bar Plot

Show:

```text
APD90
Ca transient amplitude
Diastolic calcium
Ca decay time
Peak active tension
```

## Figure 3: Tissue / Organ Map

For regional injury, show:

```text
Activation time map
Calcium peak map
Active tension peak map
```

---

# Suggested Caption

```text
A porcine ventricular Gaur-model cell was paced under baseline and synthetic xenotransplant-relevant perturbations. Conductance and calcium-handling scaling factors were used to mimic primary calcium-release failure, impaired SERCA-mediated relaxation, myofilament contractility loss, and combined rejection-like calcium injury. The test demonstrates that the model framework can represent distinct mechanisms of contractility loss and propagate them from cell-scale electrophysiology to tissue- or organ-scale function.
```

---

# Implementation Notes

## Keep the Scenarios Synthetic

Use wording such as:

```text
synthetic perturbation
capability demonstration
rejection-like phenotype
calcium-injury proxy
```

Avoid wording such as:

```text
validated rejection model
confirmed biological mechanism
clinically calibrated parameter set
```

## Recommended Code Pattern

Use the original Gaur model parameters as the baseline and multiply only the selected conductances or fluxes.

Example:

```cpp
scalar GCaL_eff = GCaL_base * scenario.GCaL;
scalar Jrel_eff = Jrel_base * scenario.Jrel;
scalar Jup_eff  = Jup_base  * scenario.Jup;
```

Then use the effective values inside the existing current and flux equations.

## Avoid Using postTxDay Directly in Fast Ionic Equations

For these tests, do not simulate days of transplantation with a millisecond ionic timestep.

Instead, use each scenario as a fixed state:

```text
baseline
early injury
calcium-release failure
rejection-like injury
regional injury
```

The biological time scale can be added later as a slow parameter-selection layer.

---

# Final Recommendation

For the first presentation to doctors, make the headline result:

```text
Calcium transient amplitude and active tension peak
```

rather than APD alone.

The key story should be:

```text
The model can distinguish whether contractility loss comes from:
1. reduced calcium entry,
2. impaired SR calcium release,
3. impaired SERCA relaxation,
4. reduced myofilament response to calcium,
5. combined rejection-like calcium injury.
```
