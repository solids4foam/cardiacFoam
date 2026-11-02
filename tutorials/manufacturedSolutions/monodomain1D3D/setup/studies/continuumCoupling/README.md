# continuumCoupling - monodomain1D3D

The junction against its continuum definition (`-reference continuum`): the
tissue and graph forcings come from the continuum kernel normalised over the
tissue part of the support and averaged over each cell, not from the solver's
own map, so the discrete map's error reaches the solution. Five questions:

1. **Map convergence.** How fast does the discrete map (support of cell
   centres within `R`, weights normalised by `sum w|K|`) converge to the
   continuum map, for each kernel? (`coupled1D3DMonodomain_map.csv`, every
   kernel from every run.)
2. **Coupled order.** At what order do the errors local to the supports, the
   graph terminal voltage and the junction current converge, for the linear,
   cosine and uniform kernels?
3. **Suppression control.** Does suppressing the junction (`R_pvj = 1e6`)
   change the support-local tissue error by order one? Under the discrete
   reference it changed the whole-domain tissue L2 error by 0.4%
   ([weakCouplingControl](../weakCouplingControl/README.md)).
4. **Charge balance.** Does the charge the graph loses through its junctions
   equal the charge the tissue receives, every staggered step?
5. **Time-step order.** At a fixed mesh, what is the observed order in `dt`
   of the coupled solution, with the junction active and suppressed?

Item 2 to 4 are repeated with each support cut by the sealed curved
endocardium of [monodomain1D3DShell](../../../../monodomain1D3DShell/README.md)
on tetrahedra.

## Ladders

Box (`monodomain1D3D`), bidirectional explicit transfer, `R_pvj = 1`, the
case's `dt ~ h^2`, `endTime 0.1`, graph `nodes<N+1>`:

- regime (a), fixed `R = 0.11`: two interleaved halving ladders
  `N = 16, 32, 64` and `24, 48, 96` (`R/h` 1.8 to 10.6), kernels uniform,
  linear, gaussian and cosine; controls on the same ladder: junction
  suppressed (linear, `R_pvj = 1e6`) and the map-based forcing
  (`-reference discrete`, linear)
- regime (b), `R = c h`: `c = 2` on `N = 16-64`, `c = 4` on `N = 32-96`
  (`R <= 1/6` keeps the support off the `y = 0` face), `c = 1` on `N = 16, 24`;
  uniform, linear, gaussian
- `R/h` scan at `N = 40`: `c = 1, 1.5, 2, 3, 4, 6`; uniform, linear, gaussian
- time step at `N = 32`: `dt = 0.1/(11 2^k)`, `k = 0..4`, linear, junction
  active and suppressed

- `R/h` guideline: the scan again at `N = 64` (`c = 1.5, 2, 3, 4, 6, 8`),
  where the terminal has the same position within its cell as at `N = 40`;
  uniform, linear, cosine (and cosine at `N = 40`)

- `R/h` guideline on tetrahedra: shell, `R = c lc` (`c = 1.5, 2, 3, 4, 6`),
  `N = 14, 20`, uniform, linear, cosine; linear also on seeds 2 and 3 at
  `N = 20`

Shell (`monodomain1D3DShell`), bidirectional, `R = 0.2`, `R_pvj = 1`: gmsh
element size `1/N`, `N = 10, 14, 20, 28`, kernels linear, cosine, uniform,
gaussian, and linear with `R_pvj = 1e6` on seed 1; linear, cosine and the
suppressed control also on seeds 2 and 3; `N = 14` once more on 4 ranks.

- Terminals inside the wall, tetrahedra (`monodomain1D3DShell -graphY
  0.5590169943749474`: terminals at `r = 0.75`, full-ball supports, `R = 0.2`),
  `N = 10, 14, 20, 28`: unidirectional linear, suppressed (`R_pvj = 1e6`) and
  bidirectional linear on seeds 1-3; unidirectional cosine, uniform and
  gaussian on seed 1; time step at `N = 14`, `dt = 0.1/(22 2^k)`, `k = 0..4`,
  junction active and suppressed
- Endocardial terminals, unidirectional linear, seed 1
- Overlapping supports and a branch point of degree three (`-graphY 0.49
  -graphY2 0.61 -radius 0.15`: four terminals, each inside its neighbour's
  support), unidirectional linear on seeds 1-3, bidirectional linear on seed 1
- Box, unidirectional linear, regime (a), `N = 16-96`

## Execution

Each run is a copy of a case in `results/` (ignored by git), one at a time.
From this folder:

```bash
mkdir -p results
box=../../..
shell=../../../../monodomain1D3DShell
run() # run <case dir> <name> [Allrun options]
{
    src=$1; name=$2; shift 2
    rsync -a --exclude setup $src/ results/$name/
    results/$name/Allrun "$@" > results/$name/log.Allrun 2>&1 || echo "FAILED $name"
}
cont="-coupling bidirectional -reference continuum"
for N in 16 24 32 48 64 96
do
    run $box lin_N$N  -N $N $cont -kernel linear
    run $box cos_N$N  -N $N $cont -kernel cosine
    run $box uni_N$N  -N $N $cont -kernel uniform
    run $box gau_N$N  -N $N $cont -kernel gaussian
    run $box off_N$N  -N $N $cont -kernel linear -rPvj 1e6
    run $box disc_N$N -N $N -coupling bidirectional -reference discrete -kernel linear
done
for k in uniform linear gaussian
do
    for N in 16 24 32 48 64; do run $box b2_${k}_N$N -N $N $cont -kernel $k -radius $(awk -v N=$N 'BEGIN{printf "%.17g", 2/N}'); done
    for N in 32 48 64 96;    do run $box b4_${k}_N$N -N $N $cont -kernel $k -radius $(awk -v N=$N 'BEGIN{printf "%.17g", 4/N}'); done
    for N in 16 24;          do run $box b1_${k}_N$N -N $N $cont -kernel $k -radius $(awk -v N=$N 'BEGIN{printf "%.17g", 1/N}'); done
    for c in 1 1.5 2 3 4 6;  do run $box scan_${k}_c$c -N 40 $cont -kernel $k -radius $(awk -v c=$c 'BEGIN{printf "%.17g", c/40}'); done
done
for k in 0 1 2 3 4
do
    dt=$(awk -v k=$k 'BEGIN { printf "%.17g", 0.1/(11*2^k) }')
    run $box dton_k$k  -N 32 -dt $dt $cont -kernel linear
    run $box dtoff_k$k -N 32 -dt $dt $cont -kernel linear -rPvj 1e6
done
for N in 10 14 20 28
do
    run $shell shlin_N$N -N $N -kernel linear
    run $shell shcos_N$N -N $N -kernel cosine
    run $shell shuni_N$N -N $N -kernel uniform
    run $shell shgau_N$N -N $N -kernel gaussian
    run $shell shoff_N$N -N $N -kernel linear -rPvj 1e6
done
for s in 2 3
do
    for N in 10 14 20 28
    do
        run $shell shlin_s${s}_N$N -N $N -seed $s -kernel linear
        run $shell shcos_s${s}_N$N -N $N -seed $s -kernel cosine
        run $shell shoff_s${s}_N$N -N $N -seed $s -kernel linear -rPvj 1e6
    done
done
run $shell shpar_N14 parallel -N 14 -kernel linear
for N in 20 14
do
    for k in linear uniform cosine
    do
        for c in 1.5 2 3 4 6; do run $shell tscan_${k}_N${N}_c$c -N $N -kernel $k -radius $(awk -v c=$c -v N=$N 'BEGIN{printf "%.17g", c/N}'); done
    done
done
for s in 2 3
do
    for c in 1.5 2 3 4 6; do run $shell tscan_linear_N20_c${c}_s$s -N 20 -seed $s -kernel linear -radius $(awk -v c=$c 'BEGIN{printf "%.17g", c/20}'); done
done
for c in 1 1.5 2 3 4 6; do run $box scan_cosine_c$c -N 40 $cont -kernel cosine -radius $(awk -v c=$c 'BEGIN{printf "%.17g", c/40}'); done
for k in uniform linear cosine
do
    for c in 1.5 2 3 4 6 8; do run $box scan64_${k}_c$c -N 64 $cont -kernel $k -radius $(awk -v c=$c 'BEGIN{printf "%.17g", c/64}'); done
done
Y=0.5590169943749474
for N in 10 14 20 28
do
    for s in 1 2 3
    do
        run $shell inlin_s${s}_N$N -N $N -seed $s -graphY $Y -coupling unidirectional -kernel linear
        run $shell inoff_s${s}_N$N -N $N -seed $s -graphY $Y -coupling unidirectional -kernel linear -rPvj 1e6
        run $shell inbi_s${s}_N$N  -N $N -seed $s -graphY $Y -coupling bidirectional  -kernel linear
        run $shell ovuni_s${s}_N$N -N $N -seed $s -graphY 0.49 -graphY2 0.61 -radius 0.15 -coupling unidirectional -kernel linear
    done
    for k in cosine uniform gaussian
    do
        run $shell in${k:0:3}_s1_N$N -N $N -seed 1 -graphY $Y -coupling unidirectional -kernel $k
    done
    run $shell enduni_s1_N$N -N $N -seed 1 -coupling unidirectional -kernel linear
    run $shell ovbi_s1_N$N -N $N -seed 1 -graphY 0.49 -graphY2 0.61 -radius 0.15 -coupling bidirectional -kernel linear
done
for k in 0 1 2 3 4
do
    dt=$(awk -v k=$k 'BEGIN { printf "%.17g", 0.1/(22*2^k) }')
    run $shell indton_k$k  -N 14 -seed 1 -graphY $Y -coupling unidirectional -kernel linear -dt $dt
    run $shell indtoff_k$k -N 14 -seed 1 -graphY $Y -coupling unidirectional -kernel linear -rPvj 1e6 -dt $dt
done
for N in 16 24 32 48 64 96
do
    run $box unilin_N$N -N $N -coupling unidirectional -reference continuum -kernel linear
done
python3 summarize_interior.py results
python3 summarize.py results lin cos uni gau off disc b2_uniform b2_linear b2_gaussian \
    b4_uniform b4_linear b4_gaussian shlin shcos shuni shgau shoff
python3 summarize.py results --scan scan_uniform scan_linear scan_gaussian scan_cosine \
    scan64_uniform scan64_linear scan64_cosine
python3 summarize.py results --dt dton dtoff
```

`summarize.py` prints, per prefix, the support-local table of terminal 0 with
orders between successive levels and a least-squares fit over the four finest
levels (with its standard error), the map table of every kernel, the `R/h`
scan, and the time-step self-convergence orders. A run whose error is not
finite or exceeds 1e3 is reported as diverged.

## Measured behaviour

Results of the runs of this study (27 October 2026; terminals inside the wall, overlapping supports and unidirectional box: 2 November 2026). Orders are
least-squares fits over the four finest levels; values in brackets are the
errors at the finest level; terminal 0.

### Map against its continuum definition (box, regime (a), N = 96)

| kernel | deposit L1 per unit current | `Wh/Wc - 1` | `<F>h - <F>c` |
|---|---|---|---|
| uniform | 0.92 +- 0.07 (7.5e-2) | 7.4e-3, erratic | 3.8e-4, erratic |
| linear | 1.80 +- 0.10 (8.1e-3) | 1.4e-5 | 1.6e-6 |
| gaussian | 1.45 +- 0.06 (1.1e-2) | 6.0e-4, erratic | 5.3e-5 |
| cosine | 1.99 +- 0.01 (7.1e-3) | 3.7e-5 | 3.2e-6 |

### Coupled solution (box, regime (a))

| kernel | support L2 | junction current | graph `V1D` at the terminal | far field L2 |
|---|---|---|---|---|
| uniform | 2.32 +- 0.12 (2.3e-2) | 0.76 +- 0.44 (6.6e-3) | 1.10 +- 0.38 (3.0e-3) | 2.00 +- 0.13 (2.6e-5) |
| linear | 2.09 +- 0.13 (9.8e-3) | 1.93 +- 0.12 (3.3e-3) | 1.91 +- 0.11 (1.5e-3) | 2.00 +- 0.01 (2.1e-5) |
| gaussian | 1.93 +- 0.03 (2.4e-2) | 1.81 +- 0.07 (5.2e-3) | 1.80 +- 0.06 (2.4e-3) | 2.00 +- 0.01 (2.1e-5) |
| cosine | 2.05 +- 0.01 (1.8e-2) | 2.02 +- 0.00 (4.1e-3) | 2.00 +- 0.01 (1.9e-3) | 2.00 +- 0.00 (2.1e-5) |

The junction quantities follow the map: erratic for the uniform kernel, near
second order for linear and gaussian, second order for cosine. The uniform
kernel's final-time current error changes sign during the run at `N = 24`
(-2.0e-2 at mid-run, +3.5e-3 at `T`), which produces its pairwise orders 8.3
and -4.8; the maximum over time (`maxStepCurrentError`) gives 1.59 +- 0.23
(uniform), 1.94 (linear), 1.84 (gaussian), 2.02 (cosine). Fitted separately,
the interleaved ladders differ for the irregular kernels (uniform current:
16-32-64 1.95, 24-48-96 -0.30; linear deposit: 1.80 and 2.03): at `N = 24,
48, 96` the terminal sits on a mesh vertex, at 16, 32, 64 inside a cell. The
gaussian kernel's jump at `R` (0.011 of its peak) gives no visible first-order
term over this range. The far field converges at second order; its constant
is the same for linear, gaussian and cosine and 20% higher for uniform
(2.6e-5 against 2.1e-5). The two terminals are mirror images and give the
same orders to 0.02.

### Suppression control

Tissue error with the junction (`R_pvj = 1`) over the error without it
(`R_pvj = 1e6`), linear kernel:

| case | within R | within 2R | whole domain | far field |
|---|---|---|---|---|
| box, N = 16 / 96 | 437 / 337 | 248 / 223 | 61 / 54 | 1.05 / 1.01 |
| shell, N = 10 / 28 | 22 / 16 | 10 / 7.3 | 2.5 / 1.8 | 1.0 / 1.0 |

With the map-based forcing (`-reference discrete`) the junction lowers the
box support error by 33-35% and the whole-domain error by 0.38%: the map's
errors cancel against the forcing, and only the junction's feedback on the
discretisation error remains.

### Charge balance

The per-step difference between the charge the graph loses and the charge
the tissue receives is 1e-14 to 6e-13 of the charge crossing the junctions in
every bidirectional run with the junction active (box and shell, serial, and
`shpar_N14` on 4 ranks). With the junction suppressed the crossing charge is
six orders smaller and the ratio reaches 5e-7. Under unidirectional transfer
the graph loses no charge by construction, so the ratio only reports the net
current of the terminals and is not a conservation check; there the deposit
integrates to the current through the normalisation `sum w|K|`.

### Time step (box, N = 32, linear)

Self-convergence orders over dt = 9.1e-3 halved four times: junction current
0.93, 1.01, 1.01; gathered tissue voltage 0.94, 1.02, 1.01; graph terminal
voltage 0.58, 0.84, 0.93. Suppressed: 1.00 throughout. The junction triples
the temporal error constant of the gathered tissue voltage (the graph sees the
tissue one step late).

### Regime (b): R = c h

| | support L2 | junction current | graph `V1D` at the terminal | far field |
|---|---|---|---|---|
| c = 2, uniform / linear / gaussian | -1.8 / -2.1 / -1.6 | 0.37 / 0.03 / 0.40 | 0.33 / -0.06 / 0.14 | 2.03 / 2.00 / 2.00 |
| c = 4, uniform / linear / gaussian | -1.1 / -1.5 / -1.5 | erratic / 0.39 / 0.35 | erratic / 0.25 / 0.27 | 2.23 / 2.03 / 2.02 |

The support-local error grows under refinement, the junction quantities stall
and the far field converges: the mesh dependence is confined to the junction.
At `N = 40` (support L2 / junction current error):

| R/h | uniform | linear | gaussian |
|---|---|---|---|
| 1 | diverged | diverged | diverged |
| 1.5 | 3.6 / 0.19 | 2.5 / 0.13 | 8.2 / 0.19 |
| 2 | 1.8 / 0.079 | 1.2 / 0.077 | 2.6 / 0.106 |
| 3 | 0.65 / 0.021 | 0.29 / 0.036 | 0.54 / 0.048 |
| 4 | 0.23 / 0.023 | 0.089 / 0.023 | 0.19 / 0.032 |
| 6 | 0.054 / 0.0105 | 0.019 / 0.0100 | 0.038 / 0.0152 |

`R = h` diverged for every kernel. The continuum estimate of the explicit
tissue-side limit, `dt q/(R^3 R_pvj chi Cm) < 2` with `q` = 0.48, 0.76, 1.29,
1.12 for the uniform, linear, gaussian and cosine kernels on a half-ball, is
not predictive at `R/h = 1.5`: `scan64_linear_c1.5` (2.16) and
`scan_gaussian_c1.5` (2.28) exceed it and stayed bounded, while
`scan64_cosine_c1.5` diverged. With about six cells in the support the
discrete `sum w^2 |K|/W^2` is not the continuum `q`. At `R = h` and `N = 24`, where
the terminal lies on cell faces, the errors are three orders smaller than at
`N = 16`: one cell per support makes the result depend on where the terminal
sits in its cell.

### R/h guideline

Junction current error relative to the current, `|I_h - I_c|/|I_c|`, at
`N = 40` / `N = 64` (same terminal position within its cell):

| R/h | uniform | linear | cosine |
|---|---|---|---|
| 1.5 | 12% / failed (a) | 8.3% / 7.5% | 11% / diverged |
| 2 | 5.1% / 4.6% | 4.9% / 4.4% | 6.1% / 5.5% |
| 3 | 1.4% / 1.2% | 2.3% / 2.0% | 3.0% / 2.6% |
| 4 | 1.6% / 1.3% | 1.5% / 1.3% | 1.8% / 1.6% |
| 6 | 0.74% / 0.67% | 0.68% / 0.63% | 0.90% / 0.80% |
| 8 | - / 0.41% | - / 0.38% | - / 0.49% |

(a) the order-one local error drove the FDA ionic state out of its domain
one step before the end; the junction was not unstable.

The relative current error, and likewise the graph terminal voltage error,
is set by `R/h`: the two meshes agree to about 10%. The support-local error
at fixed `R/h` roughly doubles from `N = 40` to `N = 64` (the shrinking
support sharpens the source), and the far field scales with `h^2` alone
(1.2e-4 to 4.7e-5). For the linear kernel the current error falls roughly as
`(h/R)^1.8`: about 4.5% at `R = 2h`, 2% at `3h`, 1.3-1.5% at `4h`, 0.65% at
`6h`, 0.4% at `8h`.

### R/h guideline on tetrahedra

Sealed shell, gmsh tetrahedra, `R = c lc` with `lc = 1/N` the gmsh element
size (`N = 14, 20`; the verifier's mean cell spacing is `0.6 lc`, so
`R/h = 1.67 c`), bidirectional transfer, `t = 0.1`; linear kernel also on
seeds 2 and 3 at `N = 20`. Indexed by the number of cells in the support,
the hexahedral and tetrahedral errors overlap; tetrahedra add a scatter
between terminals and meshes. Junction current error relative to the
current, linear kernel, range over terminals and meshes:

| cells in the support | hexahedra (cube) | tetrahedra (shell) |
|---|---|---|
| 38-59 | 2.0-2.3 % (59) | 1.4-7.6 % (38-46) |
| 90-134 | 1.3-1.5 % (134) | 1.6-3.6 % (90-99) |
| 311-452 | 0.63-0.68 % (452) | 0.57-1.5 % (311-337) |
| 739-1071 | 0.38 % (1071) | 0.17-0.60 % (739-800) |
| 2640-2824 | - | 0.02-0.12 % |

The uniform and cosine kernels follow the same pattern (uniform 2.0-7.5 %,
cosine 1.9-4.6 % at 38-46 cells). Under refinement at fixed `c` the
tetrahedral current error grows (linear, `c = 4`: 0.17-0.25 % at `N = 14`,
0.42-0.60 % at `N = 20`), as for the point-source limit on the cube; the
far field stays at the shell operator's error (`9e-4` to `1.6e-3`), except
where `2R` exceeds the far-field cutoff of 0.25.

### Sealed curved wall, tetrahedra (shell)

`T = 0.1`. Least-squares orders over all levels of the three seeds (linear,
cosine, suppressed) or of seed 1 (uniform, gaussian); values at `N = 28`
(mean over seeds); the seed-to-seed spread at `N = 28` is 1-5%.

| kernel | support L2 | junction current | graph `V1D` at the terminal | far field |
|---|---|---|---|---|
| uniform | 1.44 +- 0.09 (1.4e-2) | 1.86 +- 0.11 (3.9e-3) | 1.21 +- 0.19 | 1.59 +- 0.07 |
| linear | 1.95 +- 0.06 (3.7e-3) | 2.15 +- 0.11 (1.8e-3) | 1.61 +- 0.12 | 1.53 +- 0.02 |
| gaussian | 1.66 +- 0.11 (7.3e-3) | 1.65 +- 0.09 (3.7e-3) | 1.12 +- 0.18 | 1.51 +- 0.06 |
| cosine | 1.85 +- 0.05 (6.1e-3) | 1.83 +- 0.06 (3.0e-3) | 1.29 +- 0.06 | 1.52 +- 0.02 |
| linear, junction suppressed | 1.76 +- 0.04 (2.2e-4) | - | 1.64 +- 0.08 | 1.50 +- 0.03 |

With the support cut by the curved wall the linear kernel has the smallest
errors and the highest orders; the cosine kernel loses its box order. The
graph terminal voltage and the far field converge at the same orders with or
without the junction. These shell orders are for `T = 0.1`; the tissue
operator's order on these meshes depends on the window. With the junction
suppressed, seed 1 and the per-level `dt` of the
[sphericalShell tetConvergence](../../../../monodomainPseudoECG/sphericalShell/setup/studies/tetConvergence/README.md)
study, this case reproduces that study's whole-domain `Vm` L2 errors to
1.4e-7 relative: pairwise orders 1.66, 1.64, 1.87 at `T = 0.2` and 2.38,
1.98, 1.95 at `T = 2`, against 1.39, 1.49, 1.75 at `T = 0.1`. The
mesh-integrated normalisation exceeds the exact one by 2.0e-2 (`N = 10`) to
2.8e-3 (`N = 28`), about second order: the flat facets place part of the
mesh inside the true cavity sphere. On 4 ranks the `N = 14` run reproduces
the serial one to about 1e-10.

The continuum normalisation is taken over the true tissue part of the ball,
but the cell integrals cover the faceted mesh, so the reference deposit
itself integrates to `I (1 + delta)` with `delta = sum omega_K/Wc - 1`
(`meshIntegratedNormalisationError` in `_map.csv`): 2.0e-2, 1.1e-2, 5.2e-3,
2.8e-3 for `N = 10` to 28 (linear, mean of three seeds), the same size as the
relative junction current error at each level (2.3e-2 to 2.3e-3). On the
endocardium the junction order therefore cannot be separated from the
second-order convergence of the faceted wall; the interior terminals below
have no facet in their supports. Per terminal (mean of three seeds, fit over
all four levels): linear current 2.20 (terminal 0) and 1.89 +- 0.34
(terminal 1), graph terminal voltage 1.66 and 1.33; cosine current 1.88 and
1.77.

### Terminals inside the wall, tetrahedra

`summarize_interior.py`: `-graphY 0.5590169943749474`, terminals at
`r = 0.75`, full-ball supports, `R = 0.2`, `T = 0.1`. Means over the seeds
present (three for linear, one for the other kernels) and both terminals;
`h` is the shell verifier's mean cell size (0.059, 0.043, 0.030, 0.022);
least-squares orders over the four levels, values at `N = 28`.

| case | support cells | support L2 | junction current | graph `V1D` at the terminal | deposit L1 |
|---|---|---|---|---|---|
| unidirectional, linear | 144-3270 | 2.18 +- 0.06 (3.7e-3) | 2.26 +- 0.08 (3.5e-3) | 1.92 +- 0.05 | 2.07 +- 0.06 |
| unidirectional, cosine | 146-3306 | 1.99 +- 0.04 (5.4e-3) | 2.00 +- 0.05 (5.5e-3) | 1.93 | 1.98 +- 0.06 |
| unidirectional, gaussian | 146-3306 | 1.93 +- 0.03 (6.7e-3) | 1.90 +- 0.06 (6.8e-3) | 1.93 | 1.70 +- 0.05 |
| unidirectional, uniform | 146-3306 | 1.59 +- 0.03 (1.2e-2) | 2.25 +- 0.23 (4.4e-3) | 1.93 | 1.16 +- 0.04 |
| bidirectional, linear | 144-3270 | 2.17 +- 0.06 (3.9e-3) | 2.34 +- 0.09 (2.9e-3) | 1.75 +- 0.04 (8.4e-4) | 2.07 +- 0.06 |
| endocardial, unidirectional, linear (seed 1) | 100-2000 | 1.90 +- 0.25 | 1.92 +- 0.58 (pairwise 4.39, -0.08, 2.53) | 1.93 | 1.85 +- 0.30 |

The exact current is about 1.057; the seed-to-seed spread of the linear
current error at `N = 28` is 7.7%. Under unidirectional transfer the graph
error is the standalone graph error for every kernel. The far field and the
whole-domain error converge at 1.56 and 1.99, the shell operator's orders in
this window (1.54, 1.54 with the junction suppressed). The bidirectional
charge balance is 9.4e-14 or better.

Suppression (linear, unidirectional, error with the junction over error with
it suppressed): within `R` 18.2, 13.9, 10.8, 10.3; whole domain 2.46, 1.96,
1.65, 1.56; far field 1.03, 1.01, 1.01, 1.01 (`N = 10` to 28).

Time step at `N = 14`, junction active, `dt = 4.5e-3` halved four times:
orders of the successive differences of the junction current error 1.00,
1.00, 1.00; of the gathered tissue voltage 0.99, 1.00, 1.00. With the
junction suppressed the gathered voltage gives 0.98, 0.98, 0.99.

### Overlapping supports and a branch point

`-graphY 0.49 -graphY2 0.61 -radius 0.15`: four terminals, each 0.12 from its
neighbour and inside its support; the two chains meet at a node of degree
three. Means over the four terminals (and the three seeds for
unidirectional), support cells 62-1379, fit over `N = 10-28`:

| case | support L2 | junction current | graph `V1D` at the terminal | charge balance |
|---|---|---|---|---|
| unidirectional, linear | 2.07 +- 0.05 (1.4e-2) | 2.14 +- 0.02 (8.2e-3) | 1.92 +- 0.05 | - |
| bidirectional, linear (seed 1) | 2.08 +- 0.12 (1.6e-2) | 2.17 +- 0.06 (7.6e-3) | 1.66 +- 0.03 (2.6e-3) | 5.9e-14 |

### Unidirectional transfer on the box

Linear kernel, regime (a), `N = 16-96`, terminal 0, fit over `N >= 32`:
junction current 1.93 +- 0.12 (3.4e-3 at `N = 96`), support L2 2.10 +- 0.13
(9.6e-3), far field 2.00 (2.1e-5). The bidirectional ladder gives 1.93
(3.3e-3), 2.09 (9.8e-3) and 2.00 (2.1e-5): the tissue-side junction does not
depend on the direction of transfer. The graph terminal error is the
standalone graph error (4.7e-7 at `N = 96`).
