# Cardiac Mechanics Benchmark Problems

## Overview

These cases test the solids4foam solvers against the cardiac mechanics benchmark problems presented in [Land et al. (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707707).

There are three benchmark problems:
 - 1. Deformation of a beam;
 - 2. Inflation of an idealised ventricle;
 - 3. Inflation and active contraction of an idealised ventricle.

All three cases use the Guccione incompressible hyperelastic law, a special case of the Fung orthotropic incompressible hyperelastic law; for example, see the [FEBio benchmarks](https://febio.org/knowledgebase/case-studies/structural-mechanics/cardiac-mechanics-benchmark-problems).


## Problem 1: Deformation of a beam

### Problem definition

As described in [Land et al. (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707707), this case is summarised as:
- Geometry: the undeformed geometry is the region x ∈ [0, 10], y ∈ [0, 1], z ∈ [0, 1] mm. Constitutive parameters: transversely isotropic, $$C = 2$$ kPa, $$b_f = 8$$, $$b_t = 2$$, $$b_{fs} = 4$$.
- Fibre direction: constant along the long axis, i.e. (1, 0, 0).
- Dirichlet boundary conditions: the left face (x = 0) is fixed in all directions.
- Pressure boundary conditions: 0.004 kPa is applied to the bottom face (z = 0).

### Expected results

The Z coordinate/position of a point (10, 0.5, 1) mm (centre of the top edge on the beam's free end) should be between 4.14 mm and 4.20 mm.


## Problem 2: Inflation of an idealised ventricle

See [Land et al. (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707707).


## Problem 3: Inflation and active contraction of an idealised ventricle

See [Land et al. (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4707707).
