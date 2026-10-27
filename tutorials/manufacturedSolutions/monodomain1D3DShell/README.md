# monodomain1D3DShell

The 1D-3D junction of [monodomain1D3D](../monodomain1D3D/README.md) with each
support cut by a sealed curved wall, on tetrahedra. The myocardium is the
sealed-shell monodomain problem of
[monodomainPseudoECG/sphericalShell](../monodomainPseudoECG/sphericalShell/README.md)
on gmsh tetrahedra; a Purkinje graph crosses the cavity on the `x` axis and
ends at two terminals on the endocardium, `(-0.5, 0, 0)` and `(0.5, 0, 0)`.
The junction forcing is the continuum one (`reference continuum`, see
[monodomain1D3D](../monodomain1D3D/README.md#continuum-reference)), with the
continuum normalisation taken over the part of the support ball outside the
inner sphere, `mu >= -r/(2 r1)` about the outward radius at the terminal.

## Problem

```text
chi Cm dVm/dt = div(G grad Vm) - chi Cm Iion + I_s + S_pvj      r1 < |x| < r2
n . G grad Vm = 0                                                on |x| = r1 and |x| = r2
chi Cm dV1D/dt = d/ds(g dV1D/ds) - chi Cm Iion1D + I_s1D + I_app1D  on the graph
```

`Vm = sqrt(1 + t) F(x)` is the sealed-shell field of
`manufacturedShellMonodomainReference.H` (`r1 = 0.5`, `r2 = 1`, rotating full
conductivity tensor `G`, written cell-wise by `setManufacturedFields`).
`V1D = sqrt(1 + t) cos(pi (x + 0.5))` (`xOrigin -0.5` in both the graph and
the coupled verifier), whose slope vanishes at both terminals, so the free
graph ends need no boundary forcing. `G1D`, `g`, `chi`, `Cm`, the FDA ionic
closure and `purkinjeFibreRadius 1/sqrt(pi)` are those of `monodomain1D3D`.

The graph has its root at the centre and `ceil(N/2)` equal edges per branch
(edge length `1/(2 ceil(N/2))`), written to `constant/purkinjeGraph` by
`Allrun`. The support radius is 0.2 (`-radius`): it must exceed the local cell
length, and `r1 + R <= r2` keeps the outer sphere out of the support.

Coupling is bidirectional by default. The tissue side of bidirectional
transfer is always explicit, so the continuum forcing is the gathered
operator `I_c w/W`.

## Usage

```bash
./Allrun [parallel] [options]
./Allclean
```

| option | effect | default |
|---|---|---|
| `parallel` | `decomposePar` (4 subdomains, scotch), parallel `cardiacFoam`, `reconstructPar -latestTime` | serial |
| `-N n` | gmsh element size `1/n`; graph edge length `1/(2 ceil(n/2))` | 10 |
| `-seed n` | gmsh random seed | 1 |
| `-coupling unidirectional\|bidirectional` | PVJ `couplingMode` | bidirectional |
| `-rPvj x` | junction resistance | 1 |
| `-kernel uniform\|linear\|gaussian\|cosine` | PVJ `pvjKernel` | linear |
| `-radius x` | PVJ `pvjRadius` | 0.2 |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.1 |

Outputs: `postProcessing/3D_tet_<N>_cells.dat` (shell errors, as in
sphericalShell), `postProcessing/graph_1D_<k>_nodes.dat` (graph errors), and
the continuum-reference files `verification/coupled1D3DMonodomain_map.csv`,
`_steps.csv` and `_local.csv` described in
[monodomain1D3D](../monodomain1D3D/README.md#usage).

There is no regression test; the study and its results are in
[continuumCoupling](../monodomain1D3D/setup/studies/continuumCoupling/README.md).
