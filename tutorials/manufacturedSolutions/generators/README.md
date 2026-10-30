# Sealed-shell manufactured-solution generators

These scripts belong to the sealed-shell cases (`monodomainPseudoECG`, `bidomain` and
`eikonalECG` under `sphericalShell`). They are not part of the solver build: the
headers they write are committed under `src/verificationModels`, and the cases
compile those headers in. Run a generator only to regenerate or check a header.

Three sympy generators build exact solutions on a spherical heart shell
`r1 < |x| < r2` whose two spheres are sealed walls, with a full anisotropic
conductivity whose fibre and sheet rotate through the wall in a fixed frame,
so the conductivity is oblique to both walls. Each checks its construction and
writes a C++ reference header for its verifier; `mms_codegen.py` holds the C++
emission they share.

| generator | problem | header | verifier |
|---|---|---|---|
| `shell_monodomain_mms.py` | monodomain, `n . G grad V = 0` on both spheres | `../../../src/verificationModels/monodomainVerification/manufacturedShellMonodomainReference.H` | `manufacturedShellMonodomainVerifier` |
| `shell_bidomain_mms.py` | heart-only bidomain, `phiI` and `phiE` each sealed by its own tensor | `../../../src/verificationModels/bidomainVerification/manufacturedShellBidomainReference.H` | `manufacturedShellBidomainVerifier` |
| `shell_eikonal_mms.py` | eikonal activation time, `n . M grad psi = 0` on both spheres | `../../../src/verificationModels/eikonalVerification/manufacturedShellEikonalReference.H` | `manufacturedShellEikonalVerifier` |
| `shell_mms.py` | heart shell between a bath ball and an outer bath: `phiE` and the extracellular current continuous, `n . Gi grad phiI = 0` on both spheres, rotated heart tensors, varying bath conductivity | `../../../src/verificationModels/bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H` | `manufacturedLeadFieldECGVerifier` (`reference sphericalShell`) |

Every field is a polynomial, trigonometric or exponential function of x, y, z,
or a ratio with a positive denominator, so the formulas extend smoothly into
the band beyond each sphere where polyhedral mesh faces put cell and face
centres. A sealed field is `F = hhat + omega m`: `omega` vanishes on both
spheres and `m` blends the per-sphere corrections that make the conormal
derivative zero there exactly. The fibre frame is fixed and the fibre angle
varies with the radius: a fibre field tangent to a sphere would vanish
somewhere on it.

## Generators

`shell_monodomain_mms.py`: `V = sqrt(1 + t) F` with the FDA ionic closure and
the positive function `g`; the conductivity `G` has three distinct eigenvalues.

`shell_bidomain_mms.py`: `phiI = sqrt(1 + t) Fi`, `phiE = sqrt(1 + t) Fe`,
`Vm = phiI - phiE`, with intracellular `Gi` and extracellular `Ge` of unequal
eigenvalue ratios, so `Gi + Ge` is not proportional to `Gi`. `Fi` equals the
field of `shell_monodomain_mms.py`. The `phiE` equation carries the source
`sqrt(1 + t) (div(Gi grad Fi) + div(Ge grad Fe))`, whose integral over the
shell is zero.

`shell_eikonal_mms.py`: `psi` is the smoothed anisotropic distance from a
mid-wall stimulus point plus the sealing correction of the monodomain shell,
with `M = G/(chi cm)` and `G` the tensor of `shell_monodomain_mms.py`; the
source is the exact residual of
`-div(M grad psi) + c0 sqrt(grad psi . M grad psi) = 1 + S`. The sealed
problem fixes `psi` up to a constant, so the header carries the stimulus point
and the radius `rPin` of the ball of pinned cells.

## Checks

Run before a header is written, by all three: sealing on both spheres (30
digits), sources and gradients against finite differences, and positivity of
every denominator. The monodomain and bidomain generators also check the
richness of the tangential conormal flux on the walls and the positivity of
`g`. The bidomain generator also checks that its `Fi` equals the monodomain
field and that the `phiE` source integrates to zero (Gauss-Legendre
quadrature). The eikonal generator also locates every critical point of `psi`
(`--report` prints the table): on each sphere the Morse count (minima -
saddles + maxima) must be 2 and every point must be non-degenerate; the
source stays finite there, where the metric is not differentiable.

Each entry of `selfTestSamples` in a header holds a function index, a
parameter set (0 `defaultParams()`, 1 `alternateParams()`), the solution, the
zone, a component and the reference value; the verifier compares its compiled
functions against them at start-up.

## Usage

`--write` runs every generator check first and writes nothing if one fails.
A relative `--write` or `--check` path is taken against this directory.

Regenerate a header (from this directory):

```bash
python shell_monodomain_mms.py --write ../../../src/verificationModels/monodomainVerification/manufacturedShellMonodomainReference.H
```

```bash
python shell_bidomain_mms.py --write ../../../src/verificationModels/bidomainVerification/manufacturedShellBidomainReference.H
```

```bash
python shell_eikonal_mms.py --write ../../../src/verificationModels/eikonalVerification/manufacturedShellEikonalReference.H
```

```bash
python shell_mms.py --write ../../../src/verificationModels/bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H
```

Check that a tracked header matches its generator (exit status 1 on a
difference): the same commands with `--check`.

Run the generator checks:

```bash
python -m pytest test_shell_monodomain_mms.py test_shell_bidomain_mms.py
```

```bash
python test_shell_eikonal_mms.py
```

```bash
python test_shell_mms.py
```

Requires sympy, mpmath, numpy and scipy; each header records the sympy
version that produced it. `--check` and the header comparisons in the tests
compare text exactly, including that version line and the common
subexpression output, so they hold only for the recorded sympy version.
