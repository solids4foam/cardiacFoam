# PIMPLE / 3D-1D coupling design analysis

Date: 2026-07-12

## Scope

This note analyzes the accepted CF-024 defect without changing solver code. It
separates observed behavior from design options and records the external
numerical basis needed before choosing a repair.

## Observed current behavior

`pimpleStaggeredElectrophysicsAdvanceScheme::advance()` calls
`pimpleControl::loop()` around:

1. secondary/Purkinje coupling preparation;
2. the complete conduction-domain advance;
3. primary/myocardium coupling preparation;
4. the complete myocardium advance.

The same `pimpleControl` is then looped again inside the monodomain and
bidomain implicit diffusion overloads. Official OpenFOAM v2512 source shows
that every `loop()` call increments the shared `corr_` counter and resets it
after `nCorrPIMPLE_` is exhausted. The inner loop therefore consumes outer
correctors.

The defect is broader than nested counter ownership:

- `monodomain1DSolver::advance()` integrates its ionic ODEs for the full `dt`
  and performs a backward-Euler cable step. Calling it for every corrector
  advances mutable time-level state repeatedly.
- `myocardiumDomain::advance()` integrates the myocardial ionic ODEs for the
  full `dt` on every corrector.
- `pvjMapper::depositCoupling()` adds the volumetric current with `+=`, while
  `myocardiumDomain::prepareTimeStep()` resets the source field only once
  before the outer loop. Corrector currents therefore accumulate instead of
  replacing the previous iterate.
- The PIMPLE scheme has no `hasPotentialDomain()` branch and does not execute
  the split reaction -> global potential -> local diffusion sequence used by
  `staggeredElectrophysicsAdvanceScheme` for bath/global-`phiE` workflows.
- OpenFOAM `residualControl` monitors solver-field residuals; the present scheme
  defines no convergence residual for the actual interface unknown, the PMJ
  current.

## External numerical basis

Vergara et al., *A coupled 3D-1D numerical monodomain solver for cardiac
electrical activation in the myocardium with detailed Purkinje network*, JCP
308 (2016), DOI `10.1016/j.jcp.2015.12.016`, formulate the directly relevant
two-way Purkinje/myocardium problem. Their Algorithm 1 performs a fixed-point
iteration on the PMJ current `gamma` within each physical timestep:

1. keep the previous-time solution fixed;
2. solve the discretized myocardium problem using `gamma(k)`;
3. solve the discretized Purkinje problem using `gamma(k)`;
4. compute `gamma(k+1)` from the newly solved endpoint/tissue potentials;
5. stop when the change in `gamma` is below tolerance.

The paper proves a sufficient contraction condition for this iteration. Its
essential time-discretization property is that coupling iterations re-solve the
same `n -> n+1` discrete problems; they do not repeatedly advance the ODE/PDE
states by another `dt`.

Cardiac reaction-diffusion literature likewise treats ionic reaction and
diffusion as successive subproblems over one physical timestep. Repeating a
full reaction step is not an outer-corrector operation unless the state is
restored to the timestep checkpoint before every nonlinear iterate.

Sources:

- OpenFOAM v2512 `pimpleControl::loop()` source:
  <https://api.openfoam.com/2512/pimpleControl_8C_source.html>
- Vergara et al. open manuscript:
  <https://eprints.whiterose.ac.uk/id/eprint/92788/1/PF-monodomain-jcp-revised-II.pdf>
- Ogiermann et al., reaction-diffusion operator splitting:
  <https://doi.org/10.1002/cnm.3670>

## Design options

### A. Remove only the inner diffusion loop

Rejected as incomplete. It fixes corrector-counter consumption but still
repeats the conduction and myocardial full-timestep ODE advances and still
accumulates PVJ source terms.

### B. Run reaction/conduction once, iterate only tissue diffusion

This avoids repeated time advancement but is not a strongly coupled 3D-1D
solve. The Purkinje state is frozen, so retrograde feedback is not converged.
It is better described as a staggered/semi-implicit method, which the repository
already provides.

### C. True fixed-point coupling on PMJ current

This is the strongest evidence-based option. Introduce a coupling iteration
whose explicit interface unknown and residual are the terminal PMJ currents.
Each iteration solves both discrete-time domains from the same timestep-start
checkpoint with the latest current, replaces (does not accumulate) interface
sources, updates the PMJ current, and checks absolute/relative convergence.

This should be a coupling-specific controller rather than treating PIMPLE
field residuals as the interface convergence criterion. OpenFOAM linear solvers
and `pimpleControl` may still be used inside a domain solve where appropriate,
but one object must have one loop owner.

### D. Monolithic 3D-1D solve

Potentially strongest mathematically, but far larger in scope and inconsistent
with the repository's modular domain architecture. Not justified as the first
repair.

## Recommended staged path

1. Immediately mark the existing PIMPLE scheme experimental/unsafe, or reject
   it at runtime, until a validated replacement exists.
2. Define checkpoint/restore or non-committing trial-solve semantics for both
   myocardium and conduction domains.
3. Define `gamma` replacement semantics so each iterate starts from external
   stimulus plus the current coupling iterate, never accumulated prior
   iterates.
4. Implement a dedicated fixed-point controller with `maxIterations`, absolute
   and relative PMJ-current tolerances, and optional under-relaxation. Add
   Aitken acceleration only after the basic contraction behavior is measured.
5. Commit the converged trial states once, then update activation time and ECG
   once at the physical timestep boundary.
6. Treat bath/global-`phiE` as a separate coupled participant and specify its
   placement explicitly before enabling this scheme for bath bidomain cases.

## Required validation before implementation approval

- Exact call counts for reaction, diffusion, coupling evaluation, source
  replacement, activation update, and ECG update.
- `nOuterCorrectors = 1` equivalence with the existing staggered method for a
  configuration where both represent the same algorithm.
- PMJ-current residual histories and convergence under 1, 2, and more
  iterations.
- A two-way 1D-3D regression covering orthodromic and retrograde propagation.
- Time-step refinement showing convergence to a stable reference.
- Serial/decomposed equivalence.
- Monodomain first; bidomain and bath/global-`phiE` only after their exact
  iteration placement is specified.
- OpenFOAM v2312-v2512 compilation and numerical checks.

## Open questions for maintainer discussion

- Should the fixed-point update be Jacobi (both domains use `gamma(k)`) as in
  the cited algorithm, or block Gauss-Seidel (the second domain uses the newest
  first-domain state)? The latter may converge faster but needs its own
  stability evidence.
- Which state is checkpointed: all ionic states plus Vm and activation history,
  or can ionic reaction be computed once and held explicit as in the cited
  semi-implicit discretization?
- Is the target coupling residual absolute current, relative current, voltage
  mismatch, or a combined norm? The first implementation should expose and log
  the chosen physical units.
- Should bath/global-`phiE` participate inside every interface iteration or only
  after converged 1D-3D activation/currents?
