# Eikonal tetrahedral convergence

The normalized generic experiment is `eikonal_tet_generic`. It is the
axis-aligned conductivity case with the advection-diffusion option disabled.
Its intended matrix is:

- gradient scheme: `gauss_linear`, `least_squares`;
- nominal resolution: `10`, `20`, `40`, `80`;
- observables: `activationTime` and the worst-electrode `Phi_e` error.

The archived `sweepCases` directory currently contains the complete four-level
least-squares subset. Those eight normalized observable rows agree with the
stored reference. The four Gauss-linear case directories are absent, so their
eight observable rows cannot currently be regenerated from raw solver output;
they exist only in the consolidated reference CSV.

The normalized runner now executes both gradient schemes across all four
levels. A completed default run therefore regenerates the missing Gauss-linear
subset and produces the full 16-row field/ECG convergence result; partial runs
do not overwrite that normalized result.

The Frontal least-squares study is a separate experiment,
`eikonal_tet_frontal`, with output
`setup/results/eikonal_tet_frontal.csv`.

The isolated gradient reconstruction utility is not an eikonal convergence
run. It is registered separately as `eikonal_gradient_tet`; its runner is
`setup/studies/gradient_reconstruction/run.sh` and its normalized output is
`setup/results/eikonal_gradient_tet.csv`.
