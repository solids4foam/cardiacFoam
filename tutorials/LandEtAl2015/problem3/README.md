# Idealised Ventricle

## Overview
This 3-D problem consists of the inflation and active contraction of an idealisde ventricle, as
proposed in problem 3 of the Land et al. (2015) benchmark article.
### Mathematical Model
The ventricles undeformed geometry is defined as the ellipsoid equation seen below. 
The $\alpha$ parameter is the fibre angle which is set as $-90^\circ$ on the 
epicardinal surface to $90^\circ$ on the endocardinal surface.
$$
\vec{x} = \begin{pmatrix} x \\ y \\ z \end{pmatrix} =
 \begin{pmatrix} r_{s}\sin u\cos v \\ r_{s}\sin u\sin v \\ r_{l}\cos u \end{pmatrix}
$$
where $r_{s} = 7 + 3 t$ and $r_{l} = 17 +3t$ are derived from the transmural distance 
which varies linearly from the endocardinal surface
to the epicardinal surface.

The fibre direction is then defined by the following equation:
$$
f(u, v) = \frac{1}{|\frac{d\vec{x}}{du}|} \frac{d \vec{x}}{du} \sin(\alpha) + 
\frac{1}{|\frac{d\vec{x}}{dv}|} \frac{d \vec{x}}{dv} \cos(\alpha)
$$

There is a singularity for this problem at the apex of the ventricle.
This is dealt with by setting $v = \mathrm{atan2}(-y, -x)$ and 
setting $u= \mathrm{acos}\left(\max\left(\min\left(\frac{z}{r_l}, 1.0\right), -1.0\right)\right)$
when $u>0$ we set $u=-u$.

The code for this mathematical set up can be seen in the /applications/utilities/setFibreFields/setFibreFields.c
file. Specifically it is within the xyz function.

Active contraction: for this problem, we have an active stress given by a constant, homogenous, second
order piola kirchoff stress in the fibre direction of 60kPa.:
$$
\mathbf{T} = \mathbf{T}_{p} T_{a}\mathbf{f}\mathbf{f}^{T}
$$
Where $T_{a} = 60kPa$, $\mathbf{f}$ is the vector of fibre direction and $\mathbf{T}_{p} = 
\frac{\partial \mathbf{W} }{\partial \mathbf{E}}$
with $\mathbf{W}(\mathbf{E})$ being the strain energy function. 

The code for this can be found in the solids4foam repository under the Guccione law.
/## Instructions
### Compile `extractIdealisedVentricleResults` Utilty
The `extractIdealisedVentricleResults` utility is used to extract results from
the idealised ventricle case.
The `extractIdealisedVentricleResults` utility can be compiled with
```bash
wmake extractIdealisedVentricleResults
```

### Run the Cases
The `Allrun` runs the a mesh study for various mesh and solution procedure
configurations. The configurations are defined near the top of the `Allrun` script:
```bashopen -a "Visual Studio Code" README.md
configs=(
    "BASE=base/snes NAME=poly.hypre PETSC_FILE=petscOptions.hypre"
    "BASE=base/segregated NAME=poly.seg"
e"
)
```
where various flags are used to specify meshing and solution procedure options.
The `Allrun` script is executed as
```bash
./Allrun
```
which creates a directory for the cases called `run_<CPU_NAME>_<DATE_TIME>`, for
example, `run_Apple_M1_Ultra_20250118_151956`. When the `Allrun` script
completes, pdf plots will be available in the directory, if `gnuplot` is
installed.