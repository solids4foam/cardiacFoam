# Transmission manufactured-solution generators

Two sympy generators build exact solutions of the bath-bidomain transmission
problem: fields that satisfy potential continuity, extracellular flux
continuity and intracellular sealing on every heart/bath interface, with
rotated anisotropic heart conductivities and spatially varying bath
conductivities. Each checks its construction and writes a C++ reference
header for `manufacturedTransmissionBathBidomainVerifier`; `mms_codegen.py`
holds the C++ emission they share.

| generator | geometry | header |
|---|---|---|
| `transmission_mms.py` | oblique slab: a heart slab between two bath slabs, planar interfaces oblique to the mesh axes, periodic along them | `../manufacturedTransmissionBathBidomainReference.H` |
| `shell_mms.py` | spherical shell: a heart shell around a bath ball, inside an outer bath | `../manufacturedShellTransmissionBathBidomainReference.H` |

Checks, run before a header is written: interface conditions at 30 digits,
sources and gradients against fourth-order finite differences, boundedness,
richness of the tangential conormal flux on the interfaces, common-mode
invariance of Vm, and the ionic closure. The slab also checks periodicity and
its affine solution. The shell also checks its sources across the extension
band beyond each sphere, where polyhedral mesh faces put cell and face
centres, and the positivity of every denominator there.

Every shell field is a polynomial, trigonometric or exponential function of
x, y, z, or a ratio with a positive denominator, so each zone formula is
analytic; the only non-smooth feature of the exact solution is the normal
derivative jump across the interfaces. The fibre frame is fixed and the fibre
angle varies with the radius: a fibre field tangent to a sphere would vanish
somewhere on it. Where the sphere normal is an eigenvector of a conductivity
tensor the tangential conormal flux vanishes, so the shell's richness area
fraction is 0.7 rather than the slab's 0.8.

Each entry of `selfTestSamples` holds a function index, a parameter set
(0 `defaultParams()`, 1 `alternateParams()`), the solution (0 affine,
1 levelSet), the zone, a component and the reference value. For
`interfaceFlux` the component `c` is evaluated with the unit normal along
axis `c`. Functions without a time argument (`Gi`, `Ge`, `sigmaB`, `G`) are
sampled at the first time only; the others at two times.

`--write` runs every generator check first and writes nothing if one fails.
A relative `--write` or `--check` path is taken against this directory.

Regenerate a header (from this directory):

```bash
python transmission_mms.py --write ../manufacturedTransmissionBathBidomainReference.H
```

```bash
python shell_mms.py --write ../manufacturedShellTransmissionBathBidomainReference.H
```

Check that a tracked header matches its generator (exit status 1 on a
difference):

```bash
python transmission_mms.py --check ../manufacturedTransmissionBathBidomainReference.H
```

```bash
python shell_mms.py --check ../manufacturedShellTransmissionBathBidomainReference.H
```

Run the generator checks:

```bash
python test_transmission_mms.py
```

```bash
python test_shell_mms.py
```

Requires sympy, mpmath and numpy; each header records the sympy version that
produced it. `--check` and the header comparisons in the tests compare text
exactly, including that version line and the common subexpression output, so
they hold only for the recorded sympy version.

`shell_monodomain_mms.py` writes `monodomainVerification/manufacturedShellMonodomainReference.H`:
the heart-only monodomain solution on the spherical shell, sealed on both spheres with the
intracellular tensor of `shell_mms.py`. Its F equals the `shell_mms.py` intracellular potential at
t = 0. Tests: `test_shell_monodomain_mms.py`.

`shell_bidomain_mms.py` writes `bidomainVerification/manufacturedShellBidomainReference.H`:
the heart-only bidomain solution on the spherical shell, with the intracellular and extracellular
tensors of `shell_mms.py`, each sealed on both spheres by its own potential: Fi equals the field of
`shell_monodomain_mms.py`, and Fe is built on the time-independent part of the `shell_mms.py` heart
extracellular potential. Besides the sealing, source, denominator and richness checks of the
monodomain shell, it checks that the integral of the phiE source div(Gi grad Fi) + div(Ge grad Fe)
over the shell vanishes by Gauss-Legendre quadrature. Tests: `test_shell_bidomain_mms.py`.
