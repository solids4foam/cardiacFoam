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

Shell (`monodomain1D3DShell`), bidirectional, `R = 0.2`, `R_pvj = 1`: gmsh
element size `1/N`, `N = 10, 14, 20, 28`, kernels linear, cosine, uniform,
gaussian, and linear with `R_pvj = 1e6` on seed 1; linear, cosine and the
suppressed control also on seeds 2 and 3; `N = 14` once more on 4 ranks.

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
for c in 1 1.5 2 3 4 6; do run $box scan_cosine_c$c -N 40 $cont -kernel cosine -radius $(awk -v c=$c 'BEGIN{printf "%.17g", c/40}'); done
for k in uniform linear cosine
do
    for c in 1.5 2 3 4 6 8; do run $box scan64_${k}_c$c -N 64 $cont -kernel $k -radius $(awk -v c=$c 'BEGIN{printf "%.17g", c/64}'); done
done
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

Results of the last run of this study (27 October 2026). Orders are
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

The junction quantities follow the map: erratic and about first order for
the uniform kernel, near second order for linear and gaussian, second order
for cosine. The far field converges at second order with the same constant
for every kernel.

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
every bidirectional run (box and shell, serial and 4 ranks).

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

`R = h` is below the stability limit of the explicit tissue-side junction,
`dt q/(R^3 R_pvj chi Cm) < 2` with `q` = 0.48, 0.76, 1.29, 1.12 for the
uniform, linear, gaussian and cosine kernels on a half-ball; with `dt ~ h^2`
the minimum stable `R/h` grows as `h^(-1/3)`. At `R = h` and `N = 24`, where
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

### Sealed curved wall, tetrahedra (shell)

Least-squares orders over all levels of the three seeds (linear, cosine,
suppressed) or of seed 1 (uniform, gaussian); values at `N = 28` (mean over
seeds); the seed-to-seed spread at `N = 28` is 1-5%.

| kernel | support L2 | junction current | graph `V1D` at the terminal | far field |
|---|---|---|---|---|
| uniform | 1.44 +- 0.09 (1.4e-2) | 1.86 +- 0.11 (3.9e-3) | 1.21 +- 0.19 | 1.59 +- 0.07 |
| linear | 1.95 +- 0.06 (3.7e-3) | 2.15 +- 0.11 (1.8e-3) | 1.61 +- 0.12 | 1.53 +- 0.02 |
| gaussian | 1.66 +- 0.11 (7.3e-3) | 1.65 +- 0.09 (3.7e-3) | 1.12 +- 0.18 | 1.51 +- 0.06 |
| cosine | 1.85 +- 0.05 (6.1e-3) | 1.83 +- 0.06 (3.0e-3) | 1.29 +- 0.06 | 1.52 +- 0.02 |
| linear, junction suppressed | 1.76 +- 0.04 (2.2e-4) | - | 1.64 +- 0.08 | 1.50 +- 0.03 |

With the support cut by the curved wall the linear kernel has the smallest
errors and the highest orders; the cosine kernel loses its box order. The
graph terminal voltage and the far field converge at the sealed-shell
tetrahedral operator's own orders, with or without the junction. The
mesh-integrated normalisation exceeds the exact one by 2.0e-2 (`N = 10`) to
2.8e-3 (`N = 28`), about second order: the flat facets place part of the
mesh inside the true cavity sphere. On 4 ranks the `N = 14` run reproduces
the serial one to about 1e-10.
