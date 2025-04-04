# Solver settings for cardiac benchmark code

The cardiac benchmark paper describes most of the parameters used in this code. However, each group still has some freedom to choose the way they solve the equations. Here we describe the choices made for this code.

## Discretization Method
We use the finite element framework FEniCS (legacy) version 2019.1.0 to solve the PDEs.

## Degree of discretization Method
We chose to use $\mathbb{P}_2$ finite elements for the displacement field.

## Quadrature rules and degree
Gaussian quadrature of degree 4

### PDE solver
Use use the SuperLU_DIST sparse solver which is part of the PETSc suite with the following convergence criteria for the Newton iterations
- relative tolerance: 1e-5
- absolute tolerance: 1e-5
- maximum_iterations: 50

and the following stopping criteria for the Krylov solver

- relative tolerance: 1e-10
- absolute tolerance: 1e-10
- maximum_iterations: 1000

## Time integration

For time integration we use the Generalized $\alpha$-method {cite}`erlicher2002analysis` with $\alpha_m = 0.2$ and $\alpha_f = 0.4$ and a time step $\Delta t = 0.001$
