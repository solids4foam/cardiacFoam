# Manufactured-solution generators

One sympy generator per manufactured solution of the monodomain, bidomain,
eikonal-diffusion and lead-field ECG cases under
`tutorials/manufacturedSolutions/`, the four documents that describe how each
one is built, and the NumPy reference modules other codes can import. Each
generator checks its construction before it writes anything and emits:

- `reference/<case>.py`: a NumPy-only module with the exact fields, sources
  and, for the ECG and junction cases, the reference quadratures (all
  generators except `shell_mms.py`);
- `latex/<case>_*.tex`: the formula fragments the documents print;
- for cardiacFoam only, the generated regions (or the whole) of a C++
  reference header under `src/verificationModels/`.

The documents are `docs/monodomain/monodomain_mms.pdf`,
`docs/eikonal/eikonal_mms.pdf`, `docs/bidomain/bidomain_mms.pdf` and
`docs/leadField/leadField_mms.pdf`. They state the problems only; they contain
no solver output.

## Files

| file | what it is |
|---|---|
| `fda_monodomain_mms.py` | FDA monodomain field on the unit cube, d = 1, 2, 3, and the manufactured ionic closure |
| `rotated_monodomain_mms.py` | `sin^2` field on the unit cube with a full constant conductivity and an applied source |
| `shell_monodomain_mms.py` | monodomain field on the sealed spherical shell with a rotating anisotropic conductivity |
| `pseudo_ecg_mms.py` | pseudo-ECG integrand and Gauss-Legendre reference of the cube fields |
| `graph_monodomain_mms.py` | reaction-diffusion Purkinje graph field and the bidirectional terminal term |
| `coupled_1d3d_mms.py` | continuum 1D-3D junction: half-ball kernel sampling, junction current, deposition |
| `cartesian_eikonal_mms.py` | exponential activation time on the unit cube, d = 1, 2, 3 |
| `eikonal_ecg_mms.py` | eikonal-driven ECG integrand, templates and Gauss-Legendre reference |
| `shell_eikonal_mms.py` | activation time on the sealed shell, parametric in `c0`; critical-point report |
| `fda_bidomain_mms.py` | FDA bidomain potentials on the unit cube, d = 1, 2, 3: the FDA monodomain field with phi_e a multiple of V_m |
| `shell_bidomain_mms.py` | heart-only bidomain potentials on the sealed shell with unequal anisotropy ratios |
| `slab_leadfield_mms.py` | lead-field ECG potentials of the insulated slab and the layered torso |
| `shell_mms.py` | lead-field ECG transmission field of a heart shell between a bath ball and an outer bath (LaTeX fragments and C++ header, no NumPy module) |
| `template_mms.py` | passive monodomain field on the unit square: the starting point of a new generator (see "Creating a new manufactured solution") |
| `mms_codegen.py` | C++ emission and the shared command line |
| `mms_common/` | NumPy and LaTeX emitters (`numpy_emit.py`, `latex_emit.py`) and Gauss-Legendre quadratures (`quadrature.py`) |
| `reference/` | emitted NumPy modules, one per case |
| `latex/` | emitted LaTeX fragments |
| `docs/_common/` | document preamble, nomenclature, bibliography and figure style (`figstyle.py`) |
| `docs/monodomain/`, `docs/eikonal/`, `docs/bidomain/`, `docs/leadField/` | Quarto entry file `<equation>_mms.qmd` (front matter), its includable text `_<equation>_body.qmd`, rendered PDF, figure script `figures.py` and `figures/` |
| `docs/bidomain/bidomain_reference.py`, `docs/leadField/leadfield_reference.py` | checks that the tutorial dictionaries match the generators, and comparison of the exact fields with finished runs (`--run-check`, `--check-run`); the lead-field module also transcribes the hand-written C++ slab functions |
| `test_*.py` | generator tests; `test_release.py` holds the release checks |
| `requirements.txt` | Python dependencies |
| `parity/` | cardiacFoam-only: an OpenFOAM application comparing the C++ headers with the NumPy references |

## Reproduce

From a fresh clone, regenerate everything and check that it matches the
committed files (Python 3, Quarto and a TeX distribution; see "Setup"):

```bash
cd tutorials/manufacturedSolutions/generators
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
for g in fda_monodomain rotated_monodomain shell_monodomain pseudo_ecg \
         graph_monodomain coupled_1d3d cartesian_eikonal eikonal_ecg \
         shell_eikonal fda_bidomain shell_bidomain slab_leadfield shell \
         template; do
    .venv/bin/python ${g}_mms.py --self-check && .venv/bin/python ${g}_mms.py --emit || break
done
git status --short reference latex     # empty: the emitted files are unchanged
.venv/bin/python -m pytest -m "not slow"
.venv/bin/python -m pytest -m slow test_release.py
quarto render docs/monodomain/monodomain_mms.qmd --to pdf
quarto render docs/eikonal/eikonal_mms.qmd --to pdf
quarto render docs/bidomain/bidomain_mms.qmd --to pdf
quarto render docs/leadField/leadField_mms.qmd --to pdf
```

The slow test repeats the generator runs and the renders in a copy of this
folder with a new environment. The sections below explain each step.

## Setup

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

Every command below runs from this folder with that environment
(`.venv/bin/python`, or `source .venv/bin/activate`). `requirements.txt` pins
the versions that produced the committed outputs.

Outside Python:

- `quarto render` and the slow release test need [Quarto](https://quarto.org)
  and a TeX distribution with `pdflatex` (TeX Live or MiKTeX);
- the mesh figures need OpenFOAM and gmsh, to mesh a tutorial case; the
  parity harness needs OpenFOAM and a cardiacFoam source tree (see
  "cardiacFoam only").

The generators, their self-checks, the NumPy references and the analytic
figures need neither.

## One command per case

`--emit` runs the self-checks and writes `reference/` and `latex/` only if
they pass; `--self-check` runs the checks alone. The two options are exclusive.

```bash
python fda_monodomain_mms.py --emit
python rotated_monodomain_mms.py --emit
python shell_monodomain_mms.py --emit
python pseudo_ecg_mms.py --emit
python graph_monodomain_mms.py --emit
python coupled_1d3d_mms.py --emit
python cartesian_eikonal_mms.py --emit
python eikonal_ecg_mms.py --emit
python shell_eikonal_mms.py --emit
python fda_bidomain_mms.py --emit
python shell_bidomain_mms.py --emit
python slab_leadfield_mms.py --emit
python shell_mms.py --emit
python template_mms.py --emit
```

and the same with `--self-check`. The shell eikonal checks and the
critical-point table take a `c0`:

```bash
python shell_eikonal_mms.py --self-check --c0 5
python shell_eikonal_mms.py --report --c0 5
```

## Documents and figures

```bash
quarto render docs/monodomain/monodomain_mms.qmd --to pdf
quarto render docs/eikonal/eikonal_mms.qmd --to pdf
quarto render docs/bidomain/bidomain_mms.qmd --to pdf
quarto render docs/leadField/leadField_mms.qmd --to pdf
```

The figures are committed; regenerate them with

```bash
python docs/monodomain/figures.py
python docs/eikonal/figures.py
python docs/bidomain/figures.py
python docs/leadField/figures.py
```

(`--draft` for a fast low-resolution pass, `--outdir DIR` to write elsewhere).
Mesh figures need a meshed case: `--mesh-case NAME=PATH` reads
`PATH/constant/polyMesh`; `--help` lists the names. Mesh a copy of the
tutorial outside the repository with only the meshing steps of its `Allrun`
(write `N n;` to `local/mesh`, then `blockMesh`, or `gmsh` with
`-setnumber lc 1/n -setnumber seed 1` and `gmshToFoam`). The committed mesh
figures use the hex meshes at `N = 8` and the tet meshes at element size
`1/5` of `monodomainPseudoECG`, `monodomainPseudoECG/sphericalShell`,
`eikonalECG`, `eikonalECG/sphericalShell`, `bidomain` and
`bidomain/sphericalShell`; the lead-field mesh figures use the meshes of
`leadFieldECG` (`-problem slab|layered`) at `N = 10` and of
`leadFieldECG/sphericalShell` at `N = 8`, as made by their `Allrun`:

```bash
python docs/monodomain/figures.py --mesh-case fda_cube_hex=CUBE_HEX \
    --mesh-case fda_cube_tet=CUBE_TET --mesh-case shell_hex=SHELL_HEX \
    --mesh-case shell_tet=SHELL_TET
python docs/eikonal/figures.py --mesh-case cartesian_cube_hex=CUBE_HEX \
    --mesh-case cartesian_cube_tet=CUBE_TET --mesh-case shell_eikonal_hex=SHELL_HEX \
    --mesh-case shell_eikonal_tet=SHELL_TET
python docs/bidomain/figures.py --mesh-case fda_cube_hex=CUBE_HEX \
    --mesh-case fda_cube_tet=CUBE_TET --mesh-case shell_hex=SHELL_HEX \
    --mesh-case shell_tet=SHELL_TET
python docs/leadField/figures.py --mesh-case slab_hex=SLAB_HEX \
    --mesh-case slab_tet=SLAB_TET --mesh-case layered_hex=LAYERED_HEX \
    --mesh-case layered_tet=LAYERED_TET --mesh-case shell_hex=SHELL_HEX \
    --mesh-case shell_tet=SHELL_TET
```

A run without `--mesh-case` leaves the committed mesh figures in place.

## Tests

```bash
python -m pytest -m "not slow"
python -m pytest -m slow test_release.py
```

The slow test copies this folder without `parity/` and `.venv/` to a temporary
directory, builds a fresh environment from `requirements.txt`, runs every
generator with `--emit` and `--self-check` and renders the four documents.

## Creating a new manufactured solution

Start from `template_mms.py`, the smallest complete generator (symbolic
field, source derived as the residual, self-checks, `--emit`, no C++
header), and its test `test_template_mms.py`:

```bash
cp template_mms.py mycase_mms.py
cp test_template_mms.py test_mycase_mms.py
python mycase_mms.py --self-check
python mycase_mms.py --emit
python -m pytest test_mycase_mms.py
```

Edit every step marked `REPLACE`. The guide, from the choice of a field
whose boundary condition holds by construction to the cardiacFoam verifier,
parity entry and tutorial case, with the checklist of traps and the check
that catches each, is `docs/_common/creating.qmd`, printed as "Creating a new
manufactured solution" in every document.

## cardiacFoam only

`parity/` and the `--check PATH` / `--write PATH` options of the generators
that own a C++ header need a cardiacFoam source tree and, for `parity/`, an
OpenFOAM environment. Users of the release do not need them. `--write` runs
every self-check first and writes nothing if one fails; a relative path is
taken against this directory; `--check` exits 1 on a difference.

| generator | header (under `../../../src/verificationModels/`) |
|---|---|
| `fda_monodomain_mms.py` | `monodomainVerification/manufacturedFDAReference.H` (generated regions) |
| `rotated_monodomain_mms.py` | `monodomainVerification/manufacturedAnisotropicMonodomainReference.H` (generated regions) |
| `pseudo_ecg_mms.py` | both headers above (generated regions) |
| `cartesian_eikonal_mms.py`, `eikonal_ecg_mms.py` | `eikonalVerification/manufacturedEikonalReference.H` (generated regions) |
| `shell_monodomain_mms.py` | `monodomainVerification/manufacturedShellMonodomainReference.H` |
| `fda_bidomain_mms.py`, `slab_leadfield_mms.py` | none: `bidomainVerification/manufacturedFDABidomainReference.H` and `ecgVerification/manufacturedLeadFieldECGReference.H` are hand-written without generated regions; `test_slab_leadfield_mms.py` compares the slab module with `docs/leadField/leadfield_reference.py`, a branch-by-branch transcription of the latter |
| `shell_bidomain_mms.py` | `bidomainVerification/manufacturedShellBidomainReference.H` |
| `shell_mms.py` | `bathBidomainVerification/manufacturedShellTransmissionBathBidomainReference.H` |
| `shell_eikonal_mms.py` | `eikonalVerification/manufacturedShellEikonalReference.H` |

A generated region lies between `// BEGIN GENERATED <name>` and
`// END GENERATED <name>`; the rest of such a header is hand-written. The
shell headers are written whole. Each records the sympy version that produced
it, and `--check` compares text exactly, so it holds only for that version.

```bash
python fda_monodomain_mms.py --check ../../../src/verificationModels/monodomainVerification/manufacturedFDAReference.H
python -m pytest parity/test_parity.py
```

`parity/test_parity.py` fails, rather than skips, without OpenFOAM; the
per-generator parity tests skip with "parity harness absent (cardiacFoam-only)"
when `parity/` or the OpenFOAM environment is missing.

## Sealed-shell generators

`shell_monodomain_mms.py`, `shell_bidomain_mms.py` and `shell_eikonal_mms.py`
build exact solutions on the heart shell `r1 < |x| < r2` whose two spheres are
sealed walls, with a full anisotropic conductivity whose fibre and sheet
rotate through the wall in a fixed frame, so the conductivity is oblique to
both walls.

| generator | problem | verifier |
|---|---|---|
| `shell_monodomain_mms.py` | monodomain, `n . G grad V = 0` on both spheres | `manufacturedShellMonodomainVerifier` |
| `shell_bidomain_mms.py` | heart-only bidomain, `phiI` and `phiE` each sealed by its own tensor | `manufacturedShellBidomainVerifier` |
| `shell_eikonal_mms.py` | eikonal activation time, `n . M grad psi = 0` on both spheres | `manufacturedShellEikonalVerifier` |

Every field is a polynomial, trigonometric or exponential function of x, y, z,
or a ratio with a positive denominator, so the formulas extend smoothly into
the band beyond each sphere where polyhedral mesh faces put cell and face
centres. A sealed field is `F = hhat + omega m`: `omega` vanishes on both
spheres and `m` blends the per-sphere corrections that make the conormal
derivative zero there exactly. The fibre frame is fixed and the fibre angle
varies with the radius: a fibre field tangent to a sphere would vanish
somewhere on it.

The fields, sources and conductivities of `shell_monodomain_mms.py`,
`shell_bidomain_mms.py` and `shell_eikonal_mms.py` are written out in the
"Sealed spherical shell" section of `docs/monodomain/monodomain_mms.pdf`,
`docs/bidomain/bidomain_mms.pdf` and `docs/eikonal/eikonal_mms.pdf`.

`shell_bidomain_mms.py`: `phiI = sqrt(1 + t) Fi`, `phiE = sqrt(1 + t) Fe`,
`Vm = phiI - phiE`, with intracellular `Gi` and extracellular `Ge` of unequal
eigenvalue ratios, so `Gi + Ge` is not proportional to `Gi`. `Fi` equals the
field of `shell_monodomain_mms.py`. The `phiE` equation carries the source
`sqrt(1 + t) (div(Gi grad Fi) + div(Ge grad Fe))`, whose integral over the
shell is zero. Its command line is that of `shell_monodomain_mms.py`:
`--emit`, `--self-check`, `--write PATH` and `--check PATH`.

`shell_eikonal_mms.py`: the sealed problem fixes `psi` up to a constant, so the
header carries the stimulus point and the radius `rPin` of the ball of pinned
cells.

Checks, run by all three before a header is written: sealing on both spheres
(30 digits), sources and gradients against finite differences, and positivity
of every denominator. The monodomain and bidomain generators also check the
richness of the tangential conormal flux on the walls and the positivity of
`g`. The bidomain generator also checks that its `Fi` equals the monodomain
field and that the `phiE` source integrates to zero (Gauss-Legendre
quadrature), that `Gi` and `Ge` share one frame with unequal eigenvalue
ratios, and that the `Vm` equation and the ionic rates hold symbolically. The
eikonal generator also locates every critical point of `psi`
(`--report` prints the table): on each sphere the Morse count (minima -
saddles + maxima) must be 2 and every point must be non-degenerate; the
source stays finite there, where the metric is not differentiable.

Each entry of `selfTestSamples` in a shell header holds a function index, a
parameter set (0 `defaultParams()`, 1 `alternateParams()`), the solution, the
zone, a component and the reference value; the verifier compares its compiled
functions against them at start-up.

```bash
python shell_bidomain_mms.py --write ../../../src/verificationModels/bidomainVerification/manufacturedShellBidomainReference.H
python -m pytest test_shell_bidomain_mms.py
```

## Transmission-shell generator

`shell_mms.py` builds the exact transmission solution of a heart shell
`0.5 < |x| < 1` between the bath ball `bloodPool` and the outer bath: `phiE`
and the extracellular current are continuous and `n . Gi grad phiI = 0` on
both spheres, with rotated full heart tensors and a spatially varying bath
conductivity. Each zone formula extends smoothly across the spheres. Its
header serves `manufacturedLeadFieldECGVerifier` (`reference sphericalShell`,
case `leadFieldECG/sphericalShell`). Its fields, sources and conductivities
are written out in the "Curved conductor" section of
`docs/leadField/leadField_mms.pdf`; the geometry, the tensors and the smooth
fields are the shared fragments of the sealed shells, and the self-checks
assert that they are the same. Checks: `python -m pytest test_shell_mms.py`;
`--emit`, `--self-check`, `--check` and `--write` as for the other shell
generators; it emits no NumPy module.
