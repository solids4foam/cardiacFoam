# Purkinje S1-S2 Calibration and Restitution

This document outlines the procedure, findings, and electrophysiological phenomena observed while calibrating the 1D Purkinje cable tutorial and extracting its functional S1-S2 restitution curve for use in `driverFOAM` eikonal sweeps.

## 1. Initial Calibration (Resting CV)

The goal was to calibrate the `monodomain1DCableCV` tutorial (using the Stewart ionic model) to a baseline Conduction Velocity (CV) of **3.0 m/s**.

Through iterative testing, the `conductivity` in `constant/electroProperties` was tuned to **2.3 S/m**. When the tissue is paced from a mathematically perfect resting state (the very first beat at $t=0$), this conductivity yields a CV of exactly **3.03 m/s**, perfectly hitting the calibration target.

## 2. S1-S2 Protocol Automation (Smart Branching)

To generate a steady-state functional restitution curve, we needed to run a standard S1-S2 pacing protocol.

We automated this using a "smart branching" bash script (`run_smart_restitution.sh`) that leverages OpenFOAM's native restart capabilities to save massive amounts of compute time:

1. **Phase 1 (S1 Drive Train)**: We fired 5 S1 beats at a Basic Cycle Length (BCL) of 1000 ms. This simulated from `0.0s` to `4.25s` using `mpirun` across 6 cores. We saved the OpenFOAM field state at `t = 4.25s`.
2. **Phase 2 (S2 Branches)**: For each Diastolic Interval (DI) we wanted to test, we restored the `4.25s` checkpoint and resumed the parallel simulation, injecting the S2 premature beat and only simulating the brief ~25ms window required for the wave to propagate.

## 3. Electrophysiological Findings

During the generation of the restitution curve, two critical, physiologically accurate behaviors emerged from the Stewart model:

### A. Supernormal Conduction (Velocity Peaking)

While the 1st S1 beat traveled at the calibrated `3.03 m/s`, subsequent beats in the drive train sped up:

- **Beat 1**: 3.03 m/s
- **Beat 2**: 3.22 m/s
- **Beat 5**: 3.17 m/s
- **S2 (DI = 0.700)**: 3.37 m/s

**Why this happens:** When pacing the tissue at 1 Hz, the resting membrane potential ($V_m$) does not perfectly return to its absolute minimum (e.g., -85 mV) before the next beat arrives due to ionic memory (e.g., slight accumulation of extracellular $K^+$). Because $V_m$ sits slightly higher (less negative), the membrane is closer to the excitation threshold. It therefore requires less depolarizing current to trigger adjacent cells, resulting in a faster conduction velocity known as **supernormal conduction**.

### B. APD Prolongation and ERP Shift

The nominal resting Action Potential Duration (APD) of the Stewart model is `~290 ms`. However, pacing it 5 times at 1.0 Hz caused the APD to physiologically lengthen to **`~450 ms`**.

Because the APD prolonged, the Effective Refractory Period (ERP) pushed significantly outward. Any premature S2 beats with a DI below `0.330` fell inside the Absolute Refractory Period and naturally failed to propagate.

## 4. Final Restitution Curve

By sweeping the S2 intervals that successfully propagated outside the ERP, we isolated the steep gradient of the restitution curve. These values were hardcoded into `src/electroModels/conductionSystemModels/restitutionEikonalSolver1D/restitutionTemplates.H`:

**Diastolic Intervals (s):**
`{ 0.330, 0.350, 0.400, 0.450, 0.500, 0.700 }`

**Conduction Velocities (m/s):**
`{ 2.03,  2.43,  2.95,  3.18,  3.29,  3.37 }`

These functional values are now actively used by the eikonal solver during 3D Purkinje sweeps.
