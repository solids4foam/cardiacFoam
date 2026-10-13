# Concurrent pseudoECG and leadFieldECG controls

Validated on native OpenFOAM v2412, lightweight, against PR59 base 8bb2e420
in the managed branch codex/lead-field-ecg-integration. No production kernel
or source equations were changed for this coexistence test, and no new build
was needed. No paper files, full campaigns or PRs were changed.

## Configuration and maintained coverage

Both new monodomainLeadFieldECG/bidomainLeadFieldECG runners accept
`-pseudoECG yes|no` (default no). Allrun copies a bare observer fragment into
local/ when requested; the tracked electroProperties includes it inside
`ecgDomains`. The existing lead observer retains the name ECG; Pseudo is a
second independent observer. Their trace and verifier filenames differ.

```bash
./Allclean
./Allrun -N 20 -pseudoECG yes
./Allclean
./Allrun fine -pseudoECG yes
./regression/regressionTest.sh
```

Both maintained regression paths now enable both observers in parallel,
retain all existing source/lead numerical bounds, require pseudo output and
invoke tutorials/regression/checkDualECG.py. The checker requires matching
actual reference times and sample counts, finite nonzero reference cohorts,
functional smoke bounds and subordinate 48/96 quadrature discrepancy. It
accepts arbitrary matching sample counts, so other dt/endTime runner controls
are still usable. The two regressions passed 6/8 numerical reference checks
plus these coexistence checks; check-only validation also passed.

## What is being compared

The lead fixture is the existing synthetic two-channel signed FDA vector
functional. The pseudo fixture uses E3(1.2,.23,.61), E4(1.35,.74,.28) and
E5(1.55,.41,.83), avoiding analytically zero symmetry-plane references. Its
existing independent FDA cosine verifier uses anisotropic=false and
Gauss-Legendre orders 96 and 48.

Explicit sigmaExtracellular=0 selects the established **unnormalised case
functional** for pseudoECG. The current verifier does not apply the optional
production factor 1/(4*pi*sigmaExtracellular). Raw pseudo output must therefore
not be labelled voltage or directly compared with the lead fixture's volts.
This control does not change that normalization contract or the verifier.

The conductivity contracts also differ: pseudo uses the existing provider's
source conductivity (the monodomain propagation tensor, or the bidomain
intracellular tensor). These particular coefficients nearly coincide:
(.111453302,.121585420,.030396355) versus the bidomain native rounded
(.111453,.121585,.0303964). Lead uses its explicitly configured physical Gi/Ge
and prescribed conductivity-weighted synthetic vectors. Similar mono/bidomain
pseudo values are fixture-specific, not evidence of model equivalence.

For each observer separately, use the dimensionless cohort metric

\[
E_{rel}=\sqrt{
 \frac{\sum_{t,e}(y_{numeric,e}(t)-y_{reference,e}(t))^2}
 {\sum_{t,e}y_{reference,e}(t)^2}}.
\]

The sums include 20 common timesteps, the two signed channels for lead and
three electrodes for pseudo. Different cohorts and kernels preclude a
relative-accuracy ranking: this normalization allows their own verification
trends to be reported side by side.

## Bounded results

All four serial runs use current sealed-wall gradients, dt 1e-4 s and final
T=.002s, with 20^3 and 40^3 meshes.

| Source | Observer | Erel20^3 | Erel40^3 | Two-level observed rate |
|---|---|---:|---:|---:|
| Monodomain | lead | .03594335 | .00818636 | 2.13443 |
| Monodomain | raw pseudo | .33386692 | .06400328 | 2.38306 |
| Bidomain | lead | .03595414 | .00818907 | 2.13439 |
| Bidomain | raw pseudo | .33385223 | .06400044 | 2.38306 |

Both error cohorts decrease with refinement. These two levels do not prove
an asymptotic order, functional equivalence or equal accuracy. The pseudo
cohort's 33.4% coarse and 6.4% fine error is substantial and remains visible;
the 40% checker threshold is a smoke bound, not publication acceptance.

Per-electrode pseudo relative-RMS errors also remain visible:

| Source | Point | Erel20^3 | Erel40^3 |
|---|---|---:|---:|
| Monodomain | E3 | .470990 | .089182 |
| Monodomain | E4 | .237683 | .048241 |
| Monodomain | E5 | .460105 | .076570 |
| Bidomain | E3 | .470973 | .089179 |
| Bidomain | E4 | .237669 | .048238 |
| Bidomain | E5 | .460089 | .076567 |

The 48/96 reference differences are about 7e-17 in the raw functional,
approximately 3e-14 relative to the reference cohort and over 12 orders below
the numerical functional error. The check uses the verifier's explicit
quadrature-delta columns, since subtracting rounded published reference values
would erase those tiny differences.

For each of the four matched serial controls, adding Pseudo left native Vm,
the lead trace and the full lead reference/error series **byte-identical**.
Bidomain native phiE and phiI were also byte-identical. The cardiac dictionary
excluding observers, gradient configuration, mesh/time inputs and physical
lead conductivity inputs matched. SHA256 values and all comparisons are
recorded in the evidence archive. Only the missing bidomain 20^3 serial baseline
needed a short lead-only repeat; other baselines were reused.

At 20^3, both source families' dual MPI runs produced byte-identical recorded
lead/pseudo traces and their full verifier series versus their serial runs.
This statement applies to recorded stream precision; it does not assert
bitwise MPI equivalence of every internal native scalar.

## Evidence and paper use

The local, unversioned `docs/validation/dual-ecg-pr59` archive contains matched series,
independent quadrature summaries, source errors, effective dictionaries,
per-electrode/cohort metrics, invariance checks and implementation/data hashes.
Generated native meshes/fields and full logs remain outside the source tree.

The supported paper claim is that both passive observers coexist, preserve
source solutions, and retain their own independent verification trends with
controlled quadrature. The two-level tables are supporting evidence. A paper
should distinguish their kernels, conductivity inputs, references and units;
a matched physiological observation comparison would require the same lead
combination and conductor assumptions and is separate future work.

For the final format-2/input-guard implementation, see
[the implementation guide](LEAD_FIELD_ECG_IMPLEMENTATION.md) and
[PR validation](LEAD_FIELD_ECG_PR_VALIDATION.md). The earlier two-level
measurements concern the same source/projection kernels; regenerated format-2
fixtures and final regression evidence are reported separately.
