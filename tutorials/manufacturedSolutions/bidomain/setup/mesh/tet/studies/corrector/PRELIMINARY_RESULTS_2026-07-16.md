# Standalone bidomain iteration-control results

The lightweight solver was run on the `N=10,20,40` unit-cube Delaunay meshes
used by the tetrahedral monodomain study. The short window uses respectively
2, 9, and 36 steps with `dt ~ h^2`. Errors below are changes relative to one
outer `phiE -> Vm` sweep and no additional non-orthogonal assembly on the same
mesh.

| N | control | equation solves/step | change in $L_2(V_m)$ | change in $L_2(\phi_e)$ | change in $L_2(\phi_i)$ |
|---:|---|---:|---:|---:|---:|
| 10 | second outer sweep | 4 | +19.60% | -46.08% | -62.31% |
| 20 | second outer sweep | 4 | +15.76% | -0.31% | -32.71% |
| 40 | second outer sweep | 4 | +5.06% | -1.02% | -13.39% |
| 10 | one non-orthogonal reassembly | 4 | +5.51% | -4.63% | -4.32% |
| 20 | one non-orthogonal reassembly | 4 | +5.74% | +1.09% | -5.97% |
| 40 | one non-orthogonal reassembly | 4 | +2.63% | +0.97% | -3.30% |
| 10 | both controls | 8 | +24.05% | -44.48% | -67.97% |
| 20 | both controls | 8 | +18.95% | +2.97% | -33.21% |
| 40 | both controls | 8 | +5.95% | -0.41% | -12.88% |

The equal-cost four-solve controls are not equivalent. An outer sweep repeats
the coupled `phiE -> Vm` block and materially reduces the reconstructed
intracellular-potential error; a non-orthogonal reassembly updates each
corrected equation without alternating the two fields and has a smaller effect
on these norms. Combining them doubles cost again but does not improve the
finest-mesh manufactured errors relative to the second outer sweep alone.

This is an iteration-sensitivity result, not a replacement for the four-level
Cartesian convergence ladder. A lower MMS error is also not a universal
criterion for accepting a deferred correction: reassembly can remove
fortuitous cancellation between discretisation terms. The result supports
keeping the controls distinct in code and discussing them together in the
paper.
