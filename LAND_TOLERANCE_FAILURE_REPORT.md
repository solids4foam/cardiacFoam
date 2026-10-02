# Land spring case: tolerance check

## Question

Does tightening the solid nonlinear tolerance allow the uncapped scalar
`LandNiederer` spring-supported slab case to complete?

## Case and run

- Case: `tutorials/electromechanicsProtocols/springSupportedSlab`
- Isolated run directory: `/tmp/cardiac-spring-active-Land-scalar-rTol005-clean-20261002`
- Active tension model: scalar `LandNiederer` (no `maxLambdaRate` cap)
- Spring stiffness: `kEnds = 1e7` Pa/m
- Time step: `deltaT = 1e-5` s
- Requested end time: `0.02` s
- Solid nonlinear relative tolerance: `rTol = 0.005`
- Run mode: 6 MPI ranks; clean mesh generation and fresh decomposition

## Result

The solver started and advanced to approximately `0.0181` s. The solid solver
then reached its 1000-iteration limit and MPI aborted. The recorded run status
is `1`, so this is a failed run. Tightening `rTol` to `0.005` did not make the
case complete. Previous tolerance variations also failed near the same time.

The last sampled active tension was approximately `1.09` kPa. This is a
plausible order of magnitude for the trace at that point, but it does not
establish that the uncapped coupled simulation is numerically or
physiologically validated.

A separate 250 ms near-isometric Land run using the optional stretch-rate cap
completed. That result is a different configuration and does not demonstrate
that this uncapped spring case is stable.

## Evidence

- Solver log: `/tmp/cardiac-spring-active-Land-scalar-rTol005-clean-20261002/log.cardiacFoam`
- Run status: `/tmp/cardiac-spring-active-Land-scalar-rTol005-clean-20261002/run.status`
- Active-tension probes: `/tmp/cardiac-spring-active-Land-scalar-rTol005-clean-20261002/postProcessing/Taprobes/solid/0/Ta`

The solver log shows the maximum-iteration failure; `run.status` contains `1`.

## Interpretation

The observed failure is in the coupled solid solve, not a startup, mesh,
decomposition, or active-tension registration failure. Tolerance refinement
alone has not resolved it.

## Reproduction on a second machine

The failure reproduces on macOS (OpenFOAM v2412), with cardiacFoam `origin/main`
(`34ce49af`) built in full mode against solids4foam `electromechanical-guccione`
(`ed7f957f`, with `solidRobin`). Runs to 30 ms, 6 ranks unless stated:

| Run | Active tension | Stretch-rate cap | `deltaT` [s] | `rTol` | Result |
|---|---|---|---|---|---|
| A | `LandNiederer` | none | 1e-5 | 0.005 | fails at 18.21 ms |
| B | `LandNiederer` | none | 1e-5 | 0.02 | fails at 18.23 ms |
| E | `LandNiederer` | none | 5e-6 | 0.02 | fails at 17.64 ms |
| D | `LandNiederer` | ±20 s^-1 | 1e-5 | 0.02 | completes |
| C | `LandNiedererTWorld` (tutorial default) | ±20 s^-1 (built in) | 1e-5 | 0.02 | completes (serial) |

The failure is therefore independent of the machine, the solids4foam version
and `rTol`.

## Cause

`LandNiederer` computes the fibre stretch rate per point as
`(lambda - prevLambda)/deltaT`, from the previous step's solid solution, and
feeds it into the cross-bridge distortion terms (`A dlambda/dt`). Ta is
computed before the solid solve and held fixed during its outer iterations.
The loop lambda -> rate -> Ta -> lambda is therefore explicit:

- The solid initial residual stays near 1e-3 throughout the uncapped run,
  about 100x the `LandNiedererTWorld` run. From about 18.0 ms it alternates
  from step to step with growing amplitude (1.52e-3, 1.55e-3, 1.53e-3,
  1.57e-3, ... 1.92e-3) until a step needs more than 1000 outer iterations.
  The solid needs 11 outer iterations per step until then.
- Halving `deltaT` makes the case fail earlier, as expected for a rate that
  divides by `deltaT`.
- The failure occurs at low tension (Ta about 1.1 kPa near the stimulated end),
  so it is not a tension-magnitude problem.

`LandNiedererTWorld` computes the rate the same way but clamps it to ±20 s^-1,
which is why the tutorial default is stable.

Supplying the rate from the solver (e.g. from the solid velocity `U`) would not
remove the loop: it is the same one-step-lagged quantity. Removing the loop
requires evaluating Ta with the new-time rate, i.e. iterating the active
tension and the solid within each step.

## Resolution

`maxLambdaRate` now defaults to 20 s^-1 in `LandNiederer` and
`LandNiedererBatched`, matching the fixed safeguard of `LandNiedererTWorld`
and `LandNiedererTWorldBatched`. Set `maxLambdaRate GREAT` to run the uncapped
model.

- The cap does not change the solution before the instability: up to 17.5 ms
  the capped and uncapped runs agree (Ta to 1e-3 kPa, end displacement to
  0.003 um) and have the same residual history.
- With the new default, and no `maxLambdaRate` entry in the case, scalar
  `LandNiederer` and `LandNiedererBatched` both complete the 30 ms case.

## Related: springSupportedSlab regression

The tutorial regression uses `LandNiedererTWorld` and is not affected by the
cap. It does depend on the solids4foam version. The reference was generated
with solids4foam `d28c6527` plus `solidRobin`. With `ed7f957f`:

| Check | Result |
|---|---|
| Vm at 20 and 50 ms | unchanged |
| Spring law at 50 and 150 ms | exact |
| Ta at the slab centre, 150 ms | 18.12 vs 17.34 kPa (+4.5 %) |
| End displacement, 150 ms | 0.784 vs 0.857 mm (-8.5 %) |

The mechanical response changed between the two solids4foam versions, which
include `electroMechanicalLaw` changes (#393, #414). The regression reference
must be regenerated from the solids4foam commit that cardiacFoam pins, once
that commit contains `solidRobin`.
