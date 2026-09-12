# Removal of the GoktepeKuhl active tension model

**Status:** removed from `src/activeTensionModels/` (scalar + batched).
**Reason:** unusable as shipped, and not calibratable in the form the paper specifies.

Read this before re-adding a Göktepe–Kuhl model. The reasoning below is why it
went, not a claim that the model is worthless.

## What was in the tree

`GoktepeKuhl_2004.H` implemented Göktepe & Kuhl (2010), *Electromechanics of the
heart: a unified approach to the strongly coupled excitation–contraction
problem*, Comput Mech 45:227–243 — Eq. (46)–(47):

```
sigma_dot = eps(Phi) [ k_sigma (Phi - Phi_r) - sigma ]
eps(Phi)  = eps_0 + (eps_inf - eps_0) exp[ -exp( -xi (Phi - Phi_bar) ) ]
```

The paper's Table 4 and Fig. 3 caption fix these in **mV / ms**:

| parameter | paper | shipped code | error |
| --- | --- | --- | --- |
| argument of `eps` | `Phi` in mV | `u = (Vm+80)/100` in [0,1] | argument compressed 100x |
| `xi` | 1 mV^-1 | 1.0, applied to `u` | effectively 0.01 mV^-1 — **switch destroyed** |
| `Phi_bar` | 0 mV (u = 0.80) | 0.0, applied to `u` | threshold sits at −80 mV, i.e. at rest — never crossed |
| `eps_0` | 0.1 ms^-1 | 1.0 s^-1 | 100x too slow |
| `eps_inf` | 1 ms^-1 | 10.0 s^-1 | 100x too slow |
| `k_sigma` | 5 kPa/mV (0.005 MPa/mV) | `kTa/100` = 0.479 kPa/mV | 10.4x too small |

`timeScaleFactor()` returned `1.0` (the OpenFOAM second clock) rather than
`1000.0` (ms), which is where the 100x came from.

Net effect of the destroyed switch: `eps` varied only 4.31 → 7.23 s^-1 (1.68x)
across the whole action potential, instead of the paper's 0.1 → 1 ms^-1 (10x,
effectively Heaviside at 1/xi = 1 mV transition width).

Corroborating evidence that this was a mV→u transplant rather than a deliberate
re-parameterisation: `AV_Vm` was written in three places and never read by the
kernel, and `/100.0` was hardcoded in four places while `AC_Vr` stayed
dict-overridable (so overriding `AC_Vr` moved the numerator but not the
denominator). `NashPanfilov` does the same mapping correctly via `(Vp - Vr)`.

## What it was NOT

An earlier pass through this concluded that a correctly-implemented
Göktepe–Kuhl "collapses onto NashPanfilov". **That was wrong** and should not be
repeated. It came from comparing the two only at the paper's own fast rates
(`eps_inf` = 1 ms^-1, tau ≈ 1 ms), where both models degenerate to
`Ta ≈ kTa·u` — a filter that fast merely copies the action potential, so the
switch cannot express itself. The agreement was a degenerate-limit artifact.

Away from that regime the two are clearly distinct. Holding the excited rate and
the 10x ratio equal, so that only the threshold location and smoothness differ:

- waveforms differ by up to **~21 kPa, about 45% of full scale**, concentrated in
  the repolarisation phase;
- sweeping `Phi_bar` from −75 → +10 mV moves **RT50 from 199 → 386 ms**.

`NashPanfilov` cannot express that: Eq. (23) hard-codes its threshold at
V = 0.05. So Göktepe–Kuhl is a genuinely richer model with an independent
relaxation-onset knob.

## Why it was still removed

That extra freedom is precisely the problem. `Phi_bar` and `xi` are two extra
parameters with no data in this project to constrain them, and the calibration
difficulty is reported elsewhere for this family too. Restoring paper fidelity
gives a ~498 kPa rectangle (k_sigma = 5 kPa/mV against Göktepe–Kuhl's own demo
material, lambda = 0.5 MPa / mu = 0.2 MPa — not calibrated to real myocardium),
which is further from anything usable here than what it replaced.

The shipped numerics were only reachable *because of* the 100x error, so the
model could not be kept honestly labelled and useful at the same time.

## If you re-add it

1. Work in **mV and ms**: `timeScaleFactor()` = 1000.0, `eps` evaluated on `Vm`,
   not on normalised `u`.
2. Take `xi` and `Phi_bar` in mV^-1 / mV, and `eps_0`/`eps_inf` in ms^-1.
3. Do not reuse Nash–Panfilov's `kTa = 47.9 kPa`. It is not a physiological peak
   tension — Nash & Panfilov state it was "adjusted to produce a maximal
   shortening of approximately 25%" for their c1 = 2 / c2 = 6 kPa passive law.
   Göktepe–Kuhl's `k_sigma` and Land's `Tref = 120 kPa` are three mutually
   incompatible amplitude conventions.
4. Derive `u` from `(Vp - Vr)`, never a hardcoded 100.
5. Bring a calibration target with it. Without one, the extra knobs are the
   reason it was deleted.

## Related

- Nash & Panfilov (2004), Prog Biophys Mol Biol 85:501–522 — Eq. (22c)/(23),
  Table 1. `NashPanfilov` in this repo is a faithful implementation of these.
- Chaste's `Nash2004ContractionModel.hpp` documents the same time-scaling
  ambiguity from the other direction ("the paper suggests using t0=25.9 (?),
  which gives Ta growing too quickly at the beginning"; it ships t0 = 100).
- The Physiome CellML for Nash–Panfilov
  (`models.physiomeproject.org/workspace/nash_panfilov_2004`) has the `eps`
  branches **swapped** relative to the paper's Eq. (23) — it applies `10*e0`
  when `u < 0.05`. Do not use it as a cross-check.
