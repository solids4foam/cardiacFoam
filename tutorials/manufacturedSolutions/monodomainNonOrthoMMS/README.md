# manufacturedSolutions/monodomainNonOrthoMMS tutorial

Non-orthogonal-mesh variant of `monodomainPseudoECG`'s 3D manufactured
monodomain MMS. Reuses the exact same exact solution, diagonal
conductivity tensor, and ionic-model forcing — only the mesh geometry
changes.

## Purpose

Verifies that OpenFOAM's non-orthogonal `Gauss linear corrected`
Laplacian scheme sustains second-order spatial convergence on a
smoothly skewed mesh, at several distortion severities. Closes the
"all MMS runs are on orthogonal Cartesian boxes" gap from the Paper I
methods review — see `[[project_paperI_methods_review]]` memory.

## How the mesh is distorted

`setup/distort_mesh.py` displaces every interior blockMesh vertex with
a smooth, boundary-vanishing sinusoidal perturbation (Salari-Knupp
style):

    dx = A*h*sin(pi*x)*sin(pi*y)*sin(pi*z)
    dy = A*h*sin(2*pi*x)*sin(pi*y)*sin(pi*z)
    dz = A*h*sin(pi*x)*sin(2*pi*y)*sin(pi*z)

`A` is a dimensionless severity knob (0 = untouched mesh), `h=1/N` is
the nominal cell size at that resolution. The perturbation is
identically zero on all six domain faces, so the box shape and the
no-flux boundary compatibility with `V_ex` are preserved exactly.

## Running the sweep

```bash
AMPLITUDES="0.0 0.05 0.10 0.15 0.20" bash setup/run_nonortho_sweep.sh
```

`AMPLITUDES` is a space-separated list of distortion severities (default:
`0.0 0.05 0.10 0.15 0.20` if unset). Requires OpenFOAM sourced first
(`source /Volumes/OpenFOAM-v2412/etc/bashrc`). Results land in
`setup/results/<A>/`, and a combined `setup/results/summary.csv` is written
at the end via `setup/summarize_results.py`.
