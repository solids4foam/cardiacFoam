# monodomain/coupled1D3DShell

## What it does

The 1D-3D junction of [monodomain/coupled1D3DCube](../coupled1D3DCube/README.md) on
tetrahedra, with the supports cut by a sealed curved wall, inside the wall, or
overlapping. The myocardium is the sealed-shell monodomain problem of
[monodomain/sphericalShell](../sphericalShell/README.md) on gmsh tetrahedra,
and the junction forcing is the continuum one; the case writes the shell and
graph errors and the junction errors of each terminal against that reference.

- **Exercises:** the Purkinje graph to tissue junction where the support is
  cut by the sealed curved wall (the default, terminals on the endocardium),
  lies inside the wall (`-graphY`) or overlaps another support (`-graphY2`,
  which also adds a graph branch point of degree three), with bidirectional
  transfer by default, together with the sealed-wall operators on
  tetrahedra. Only this case's regression checks the junction against its
  continuum definition, including the charge balance across the junctions.
- **Use it for:** verifying a change to the PVJ coupler, its kernels or its
  normalisation where a support meets a wall or another support; measuring
  the junction error against the continuum definition (the
  [continuumCoupling](../coupled1D3DCube/setup/studies/continuumCoupling/README.md)
  study runs here and on the cube).

## Regression

Tier: `standard`.

`regression/regressionTest.sh` runs `./Allrun` in serial with the case
defaults only (gmsh tetrahedra `N = 10`, seed 1, two terminals on the
endocardium, bidirectional transfer, continuum reference, linear kernel,
`pvjRadius 0.2`, `T = 0.1`); `-graphY`, `-graphY2` and the other kernels are
not run. It compares 41 rows: `regression/coupled1D3DShell.reference` (32
rows) holds the cells per direction, the final time and 13 shell error norms
(`Vm`, `u1`, `u2`, the `Vm` bulk and boundary rows and two wall-gradient
rows); from `verification/coupled1D3DMonodomain_local.csv`, for each terminal
the support error, the graph terminal-voltage error, the junction-current
error and its largest step value, and the largest charge imbalance (expected
0 within `1e-12`); from `verification/coupled1D3DMonodomain_map.csv`, for
each terminal the linear-kernel support cell count, normalisation, gather and
deposit errors. `regression/coupled1D3DShell.graph.reference` (9 rows) holds
the L1, L2 and Linf errors of `Vm1D`, `u1` and `u2` on the graph. The expected
values are frozen from an earlier run on the gmsh mesh named in
`shell.geo.template` (gmsh 4.15.2, seed 1, one thread), with absolute
tolerances of about `1e-5` of each value (the cell counts are exact). The
marker `regression/needs-tools` (gmsh) applies: without gmsh the script
exits 77, an expected skip.

```bash
./regression/regressionTest.sh
```

## Problem

By default a Purkinje graph crosses the cavity on the `x` axis and ends at
two terminals on the endocardium, `(-0.5, 0, 0)` and `(0.5, 0, 0)`; `-graphY`
and `-graphY2` move them (below). The junction forcing is the continuum one
(`reference continuum`, see
[monodomain/coupled1D3DCube](../coupled1D3DCube/README.md#continuum-reference)), with the
continuum normalisation taken over the part of the support ball outside the
inner sphere, `mu >= -r/(2 r1)` about the outward radius at the terminal.

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
closure and `purkinjeFibreRadius 1/sqrt(pi)` are those of `monodomain/coupled1D3DCube`.

The graph has its root at the centre and `ceil(N/2)` equal edges per branch
(edge length `1/(2 ceil(N/2))`), written to `constant/purkinjeGraph` by
`Allrun`. The support radius is 0.2 (`-radius`): it must exceed the local cell
length, and `r1 + R <= r2` keeps the outer sphere out of the support.

Two variants move the terminals off the endocardium:

- `-graphY y` puts the graph on the line `(x, y, 0)`, with terminals at
  `(-0.5, y, 0)` and `(0.5, y, 0)` inside the wall, `r = sqrt(0.25 + y^2)`.
  When `r - R >= r1` and `r + R <= r2` the support is a full ball (the
  verifier's `interior` support): no facet of either sphere enters it, so the
  continuum normalisation and the mesh agree up to the cell clipping.
  `y = 0.5590169943749474` gives `r = 0.75`, the centre of the wall.
- `-graphY2 y2` adds a second chain on `(x, y2, 0)` with its own two
  terminals, joined to the first root by one edge along `y` at `x = 0`, where
  `V1D` vanishes. `V1D` depends on `x` only, so it is zero along the whole
  joining edge, which carries no current; at each joined root the two
  branches of its chain keep the flux balance of the straight chain, so the
  node, a branch point of degree three, satisfies the graph equations with the
  same field. With `|y2 - y| < R` each terminal lies
  in its neighbour's support and the supports overlap.

`V1D` depends on `x` only, so both variants keep the graph field and its
closure unchanged.

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
| `-kernel uniform\|linear\|gaussian\|cosine` | PVJ `pvjKernel` (`uniform`: top-hat; `gaussian`: truncated Gaussian) | linear |
| `-radius x` | PVJ `pvjRadius` | 0.2 |
| `-graphY y` | graph on the line `y`: terminals inside the wall | 0 (on the endocardium) |
| `-graphY2 y` | second chain on the line `y`, joined at `x = 0` | none |
| `-dt x` | time step | `endTime/ceil(endTime/(0.00892857 (10/N)^2))` |
| `-endTime x` | end time | 0.1 |

Outputs: `postProcessing/3D_tet_<N>_cells.dat` (shell errors, as in
sphericalShell), `postProcessing/graph_1D_<k>_nodes.dat` (graph errors), and
the continuum-reference files `verification/coupled1D3DMonodomain_map.csv`,
`_steps.csv` and `_local.csv` described in
[monodomain/coupled1D3DCube](../coupled1D3DCube/README.md#usage).

The study and its results are in
[continuumCoupling](../coupled1D3DCube/setup/studies/continuumCoupling/README.md).
