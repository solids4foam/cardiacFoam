# Case description
This is the cardiac tissue electrophysiology simulation benchmark from Niederer et al. (2011) Verification of cardiac tissue electrophysiology simulators using an N-version benchmark, Phil. Trans. R. Soc. A (2011) 369, 4331–4351 doi:10.1098/rsta.2011.0139.

The case consists of a cuboid domain, which is externally stimulated in one corner, and the time taken for the depolarisation wave to reach the far corner is measured. Niederer et al. (2011) present the results from eleven separate groups, some of which use finite element implementations and some of which use finite difference codes; none use the finite volume method; as such, this would appear to be the first application of the finite volume method to this problem.

# Geometry
20 × 7 × 3 mm cuboid.

# Mesh
The benchmark proposes three cell widths:
 - 0.5 mm
 - 0.2 mm
 - 0.1 mm

For the current case, these can be changed in `system/blockMeshDict`.

# Time-step
The benchmark proposes three time-step sizes:
 - 0.05 ms
 - 0.01 ms
 - 0.005 ms

For the current case, these can be changed in `system/controlDict`.

# Governing equations

The electrophysiology is governed by a scalar reaction-diffusion partial differential equation derived from the monodomain assumption. The primary unknown is the scalar transmembrane voltage, representing the voltage difference between the inside and outside of the cardiac cells. The reaction term is governed by an ionic model, which expresses a local current in terms of the transmembrane voltage and internal "gating" variables; these ionic models are governed by multiple ordinary differential equations. The Niederer et al. (2011) benchmark stimulates the use of the ten Tusscher & Panfilov ionic model; however, as a starting point, the current case uses the simpler, minimal ionic model by Bueno-Orovio, governed by three ordinary differential equations. The Bueno-Orovio parameters have been chosen to mimic the more complex ten Tusscher & Panfilov model.

# Fields

The initial conditions and boundary conditions for the primary unknown - transmembrane voltage `Vm` - are found in `0/Vm`.

The electrophysiological and ionic model parameters are given in `constant/electroActivationProperties`; in addition, numerical settings for solving the ionic model equations are found in the same file.

# Expected results
Activation times are plotted on the plane with normal (-1.05e-05 1.4E-4 3.E-5) and origin (0 0 0.007).

# How to run the case
To run the case, execute the included `Allrun` script, e.g.

    > ./Allrun

The `Allrun` script creates the mesh using `blockMesh` and then runs the `electroActivationFoam` solver. Alternatively, the mesh and solver can be run manually, e.g.

    > blockMesh
    > electroActivationFoam

For example, after updating the mesh density in `system/blockMeshDict`, the mesh can be regenerated with `blockMesh`. The results and mesh are cleared using the included `Allclean` script:

    > ./Allclean
