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
alone has not resolved it. The next diagnostic should compare the failed
solid residual and displacement iteration history with the spring load,
stretch-rate input, and active-tension evolution around `18.1` ms. The
rate-capped run should remain identified as a separate stabilization
sensitivity, not silently substituted for the uncapped model.
