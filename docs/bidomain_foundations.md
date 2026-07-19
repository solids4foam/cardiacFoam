# Mathematical and Numerical Foundations of the Cardiac Bidomain Equations: Unstructured Finite Volume Schemes, Convergence Dynamics, and Interface Boundary Layers

## Mathematical Formulation and Homogenization of the Bidomain Model

The electrical activity of cardiac tissue is characterized by highly cooperative, multiscale physiological dynamics. To simulate the propagation of cardiac action potentials at the tissue and organ scales, the bidomain model serves as the standard continuum framework. This model is derived from the microscopic cellular structure of the myocardium through mathematical homogenization, specifically volume averaging or unfolding operators. This homogenization processes the discrete cellular arrangement—comprising individual cylindrical myocytes connected by high-resistance gap junctions—into a homogenized macroscopic continuum where both the intracellular and extracellular spaces exist co-locally at every point in the cardiac domain.

Mathematically, the intracellular potential $\phi_i$ and the extracellular potential $\phi_e$ are defined across the entire cardiac domain, and their difference yields the transmembrane potential $V_m$:

$$V_m = \phi_i - \phi_e$$ [cite: 3, 4, 8]

By asserting that electrical charge is conserved and modeling the cell membrane as a parallel combination of a capacitor and non-linear ion channels, the bidomain equations are formulated as a coupled system of partial differential equations (PDEs):

$$\chi \left( C_m \frac{\partial V_m}{\partial t} + I_{\text{ion}}(u, V_m) \right) - \nabla \cdot \left(\mathbf{M}_i \nabla \phi_i\right) = I_i$$ [cite: 8, 11]

$$\chi \left( C_m \frac{\partial V_m}{\partial t} + I_{\text{ion}}(u, V_m) \right) + \nabla \cdot \left(\mathbf{M}_e \nabla \phi_e\right) = -I_e$$ [cite: 8, 11]

Within this system, $\chi$ represents the membrane surface-to-volume ratio, $C_m$ denotes the specific membrane capacitance, and $I_i$ and $I_e$ represent external current stimuli applied to the intracellular and extracellular spaces, respectively. The total transmembrane ionic current density $I_{\text{ion}}$ depends on the transmembrane voltage and a set of internal state variables $u$ governed by ordinary differential equations (ODEs) describing ion channel gates and intracellular concentrations.

The anisotropic conduction properties of the tissue are represented by the intracellular and extracellular conductivity tensors, $\mathbf{M}_i$ and $\mathbf{M}_e$, which are constructed from the local alignment of cardiac muscle fibers. These tensors are symmetric and coercive, meaning their spectrum is strictly bounded away from zero:

$$\lambda_{-} \vert\xi\vert^2 \le \mathbf{M}(\mathbf{x}) \xi \cdot \xi \le \lambda_{+} \vert\xi\vert^2 \quad \forall \xi \in \mathbb{R}^d, \quad \text{with} \quad 0 < \lambda_{-} \le \lambda_{+} < +\infty$$ [cite: 14, 15]

A central complexity of the bidomain model is that these conductivity tensors have unequal anisotropy ratios, meaning that $\mathbf{M}_i$ is not proportional to $\mathbf{M}_e$. Under the highly restrictive and physically unrealistic assumption of equal anisotropy ratios ($\mathbf{M}_i = k \mathbf{M}_e$), the system simplifies to the monodomain model, which is represented by a single reaction-diffusion equation for $V_m$. However, the presence of unequal anisotropy ratios requires solving the fully coupled bidomain system. This inequality introduces numerical challenges, as any external electrical field applied to the tissue induces complex virtual electrode polarization (VEP) patterns deep within the bulk of the myocardium.

## Spatial Discretization Stencils: Structured vs. Unstructured Meshes

Selecting the spatial mesh topology is a critical factor in determining the physical accuracy and computational efficiency of cardiac simulations. Historically, structured Cartesian grids have been widely used due to the simplicity of their uniform stencils and their compatibility with standard finite difference methods. Despite their ease of implementation, structured grids introduce severe limitations when resolving the complex, curved surfaces of the human ventricles.

A primary issue with structured grids is the generation of a jagged, staircase representation of curved cardiac boundaries. In simulations of clinical interventions, such as high-voltage defibrillation shocks, these artificial geometric steps act as unphysical conductive discontinuities. These jagged steps produce spurious, localized electrical potential gradients and artificial VEP hotspots along the boundaries, which can compromise the integrity of wave initiation and propagation predictions.

To overcome these geometric limitations, unstructured meshes composed of triangles, tetrahedra, or triangular prisms are used. These unstructured elements conform smoothly to curved anatomical boundaries, eliminating boundary artifacts and allowing the computational domain to represent patient-specific geometries accurately. Additionally, the myocardium is organized into layers of fibers that rotate through the wall thickness. Unstructured meshes allow the local coordinate systems of individual elements to align with these rotating fiber orientations, providing a more accurate representation of the anisotropic conductivity tensors.

On a microstructural level, cardiac tissue is not a continuous syncytium but is divided by collagenous cleavage planes and laminar sheets of fibers. These cleavage planes act as barriers to intracellular current while allowing extracellular current to flow freely, significantly affecting macroscopic conduction velocity and wave front safety. Standard finite difference and finite element methods struggle to capture these microstructural details without excessive grid refinement. In contrast, unstructured finite volume methods can represent cleavage planes as zero-volume interfaces between adjacent control volumes. This allows the resistance of these intercellular clefts to be varied explicitly, enabling microstructural simulations that are computationally efficient to construct and solve.

To construct conforming unstructured grids for complex geometries, specific geometric strategies are applied to the mesh elements. For example, in implementations using hexahedral finite volumes, quadrilateral sub-faces are split into triangles. This sub-triangulation enhances the accuracy of the numerical integration over complex, curved surfaces and preserves the fidelity of the rotating fiber fields within the numerical flux calculations.

| Attribute | Structured (Cartesian) Grids | Unstructured Conforming Meshes |
| :--- | :--- | :--- |
| **Boundary Fidelity** | Jagged, staircase representation of curved surfaces | Conforming, smooth representation of anatomy |
| **Anisotropy Representation** | Approximated on uniform grid directions; prone to grid-alignment bias | Aligns naturally with rotating fiber vectors |
| **Microstructural Modeling** | Incapable of representing cleavage planes without extremely fine resolution | Cleavage planes represented as zero-volume interfaces |
| **Artifact Generation** | Generates false virtual electrode polarizations at boundary steps | Mimics the physical potential distribution accurately |
| **Matrix Properties** | Uniform band structure; highly optimized memory access | Irregular sparsity patterns; requires indirect indexing |

## Cell-Centered Finite Volume Schemes: Consistent Flux Approximations

### The K-Orthogonality Constraint and TPFA Limitations

In cell-centered finite volume methods, the computational domain $\Omega$ is partitioned into a set of non-overlapping control volumes $K$. Integrating the anisotropic diffusion equation over each cell $K$ and applying the divergence theorem yields the sum of normal fluxes across the faces $\sigma$ that form the boundary $\partial K$:

$$\int_K \nabla \cdot (\mathbf{M} \nabla u) \, dx = \sum_{\sigma \in \partial K} \int_{\sigma} (\mathbf{M} \nabla u) \cdot \mathbf{n}_{K,\sigma} \, d\gamma$$ [cite: 15, 21, 22]

Here, $\mathbf{n}_{K,\sigma}$ is the outward unit normal vector to the face $\sigma$.

The classical Two-Point Flux Approximation (TPFA) assumes that the normal diffusive flux across a shared face $\sigma$ between adjacent cells $K$ and $L$ depends only on the values at the cell centers $u_K$ and $u_L$:

$$F_{K,\sigma} = T_{\sigma} (u_L - u_K)$$

Where $T_{\sigma}$ is the face transmissibility coefficient. This simple formulation is consistent and converges to the exact solution only if the mesh satisfies the K-orthogonality constraint. This constraint requires that the vector connecting the cell centers $x_K$ and $x_L$ is parallel to the conormal direction $\mathbf{M} \mathbf{n}_{K,\sigma}$ across the interface.

When simulating cardiac tissue on unstructured meshes, the K-orthogonality condition is violated due to geometric element distortion and the strong anisotropy of the conductivity tensors $\mathbf{M}_i$ and $\mathbf{M}_e$. Under these conditions, the TPFA scheme exhibits an $O(1)$ consistency error that does not vanish with mesh refinement. This inconsistency can lead to severe grid-alignment errors, artificial conduction blocks, and inaccurate conduction velocities.

### Multi-Point Flux Approximation (MPFA) Formulations

To recover consistency on general unstructured meshes, Multi-Point Flux Approximation (MPFA) schemes express the flux across a face $\sigma$ as a linear combination of the potentials from a larger set of neighboring cell centers. The O-method and the L-method are the primary MPFA variants used to handle non-orthogonal grids and anisotropic media:

- **The MPFA O-method** constructs local interaction regions around each mesh vertex. By introducing auxiliary degrees of freedom on the faces and enforcing continuity of both the potential and the normal flux across the sub-faces within each interaction region, the face potentials are eliminated. This results in a multi-point flux formula where the stencil for a single face involves all cells sharing the vertex. While the O-method is highly consistent, it can violate the discrete minimum and maximum principles on distorted grids or under highly anisotropic conductivities, leading to spurious oscillations and unphysical potentials.
- **The MPFA L-method** is designed to address the monotonicity issues of the O-method, particularly in highly anisotropic media. The L-method minimizes the size of the flux stencil by selecting the closest neighboring cells that align with the dominant principal direction of anisotropy. This approach reduces the cell stencil from nine points to seven points on quadrilateral meshes, which expands the monotonicity range and helps suppress non-physical oscillations.

$$\begin{aligned} \text{MPFA O-method Stencil (Quadrilateral):} \quad & 6\text{-point face flux stencil}, \quad 9\text{-point cell stencil} \\ \text{MPFA L-method Stencil (Quadrilateral):} \quad & 4\text{-point face flux stencil}, \quad 7\text{-point cell stencil} \end{aligned}$$ [cite: 26]

These multi-point flux methods are mathematically proven to be first-order convergent in both the potential and the flux fields on general distorted meshes.

### Parallel Solver Performance: MPFA vs. TPFA

While MPFA schemes are more complex and require larger stencils than TPFA, their integration into high-performance parallel solvers is highly feasible. Computational studies using decoupled parallel Schwarz preconditioned solvers show that the performance of MPFA on unstructured meshes is comparable to that of TPFA. The additional computational overhead of assembling and solving the multi-point stencil is offset by the improved physical accuracy and the reduction of mesh-induced conduction blocks.

## Discrete Duality Finite Volume (DDFV) Framework for the Bidomain Model

The Discrete Duality Finite Volume (DDFV) method provides a robust framework for discretizing diffusion operators on highly distorted or non-conforming meshes without requiring K-orthogonality. This method achieves strong consistency and stability by solving the governing equations simultaneously on multiple, coupled meshes.

### Mesh Structure and Operator Construction

The DDFV spatial discretization is based on three associated meshes:
1. A primal mesh $\mathcal{T}$ consisting of almost arbitrary polygons.
2. A dual mesh $\mathcal{T}^*$ formed by connecting the centers of the primal cells, with vertices located at the primal cell centers.
3. A third mesh $\mathfrak{D}$ consisting of "diamond cells" $\mathcal{D}_{\sigma}$ constructed around the edges of the primal and dual meshes.

```
       Primal Cell Center (u_K)
                o (x_K)
               / \
              /   \
  Dual Vertex/     \ Dual Vertex
    (u_F1)  o-------o (u_F2)
    (y_F1)  \  \D_s / (y_F2)
             \     /
              \   /
               \ /
                o (x_L)
       Primal Cell Center (u_L)
```

Degrees of freedom are located at both the primal cell centers $x_K \in \mathcal{T}$ and the dual mesh vertices $y_{F} \in \mathcal{T}^*$. For each primal edge $\sigma$ shared by cells $K$ and $L$ with endpoints $y_{F_1}$ and $y_{F_2}$, the corresponding diamond cell $\mathcal{D}_{\sigma}$ is the quadrilateral defined by the vertices $(x_K, y_{F_1}, x_L, y_{F_2})$ ordered clockwise.

The discrete gradient of a scalar field $u_h$ on the diamond cell $\mathcal{D}_{\sigma}$ is defined as:

$$\nabla_{\mathcal{D}} u_h = \frac{1}{2 \vert\mathcal{D}_{\sigma}\vert} \left[ \left(u_L - u_K\right) \vert\sigma\vert \mathbf{n}_{K,L} + \left(u_{F_2} - u_{F_1}\right) \vert d_{K,L}\vert \mathbf{n}_{F_1,F_2} \right]$$ [cite: 31]

Here, $\vert\mathcal{D}_{\sigma}\vert$ is the area of the diamond cell, $\vert\sigma\vert$ is the length of the primal edge, and $\vert d_{K,L}\vert$ represents the distance between the primal centers $x_K$ and $x_L$. The vectors $\mathbf{n}_{K,L}$ and $\mathbf{n}_{F_1,F_2}$ are the unit normals to the dual and primal edges, respectively.

This gradient formulation is exact for linear functions. The discrete divergence operator for a vector field $\mathbf{p}_h$ is constructed as the formal discrete adjoint of the negative gradient operator. This construction ensures that a discrete Green's formula (integration by parts) is satisfied exactly at the discrete level:

$$\sum_{K \in \mathcal{T}} \vert K\vert (\text{div}_{\mathcal{T}} \mathbf{p}_h)_K u_K + \sum_{D \in \mathfrak{D}} \vert D\vert \mathbf{p}_D \cdot (\nabla_{\mathcal{D}} u_h)_D = \sum_{\sigma \in \partial\Omega} \vert\sigma\vert (\mathbf{p}_{\sigma} \cdot \mathbf{n}) u_{\sigma}$$ [cite: 29, 31]

### Equivalence to Non-Conforming Finite Elements

The DDFV method is closely linked to finite element formulations. Specifically, this finite volume scheme is mathematically equivalent to a non-conforming finite element method where the basis functions are defined on the diamond cells $\mathfrak{D}$ of the third mesh. Under general geometric conditions on these diamond cells, this equivalence allows the application of classic finite element error analysis techniques to the finite volume formulation.

### Convergence Rates and Stabilization

The DDFV scheme has been proven to converge to weak solutions of the cardiac bidomain equations on highly distorted and non-conforming meshes. Under regular mesh conditions, the convergence rates exhibit specific asymptotic behavior in different norms:
- **$H^1$ Norm Convergence**: First-order accuracy ($O(h)$) is achieved for both the potentials and the discrete gradients, which is essential for capturing the steep spatial activation wavefronts of cardiac depolarization.
- **$L^2$ Norm Convergence**: While theoretical analysis guarantees first-order convergence, numerical experiments consistently demonstrate second-order accuracy ($O(h^2)$) for the potentials on general, unstructured grids.
- **Superconvergence**: On homothetically refined grids, superconvergence is observed, where the error in the gradients reduces at a rate higher than $O(h)$.

To maintain stability on highly distorted grids, DDFV schemes often incorporate a penalization term in the discrete equations. This stabilization penalizes the difference between the primal and dual representations of the potentials, forcing them to converge to the same weak limit and ensuring the overall energy stability of the system.

## Control Volume Finite Element Methods (CVFEM) and Vertex-Centered Formulations

The Control Volume Finite Element Method (CVFEM), also known as the Finite Volume Element Method (FVEM), combines the geometric flexibility of finite element methods with the conservative properties of finite volumes.

### The Dual Grid Structure

In a vertex-centered CVFEM, the computational domain is first discretized using a primal triangular or tetrahedral finite element mesh $\mathcal{T}_h$. A secondary dual mesh is then constructed, typically a Donald dual mesh or median dual mesh, by connecting the barycenters of the primal elements to the midpoints of their edges. The control volumes $K$ are associated with the vertices of the primal mesh, and the degrees of freedom are located at these vertices, as in standard conforming finite element methods.

```
          Primal Vertex (u_1)
                 o
                /|\
               / | \
              /  |  \
             /  (B)  \  Primal Element
            /  /   \  \
           /  /     \  \
          o--o-------o--o
      Primal Vertex    Primal Vertex
         (u_2)            (u_3)
```

*(B) denotes the element barycenter. The dashed lines inside the element represent boundaries of the dual control volume surrounding each vertex.*

To solve the equations, the convective and diffusive terms are integrated over each control volume. This integration requires an interpolation function that describes how the variables vary within each finite element. Typically, a linear function is used for the diffusive terms, while an exponential function of the local Peclet number is used for the convective terms to maintain stability.

### The Mass Matrix Formulation

A key advantage of the CVFEM framework for bidomain simulations is its formulation of the mass matrix. In cardiac propagation modeling, standard finite difference, cell-centered finite volume, and mass-lumped finite element methods can introduce significant errors in conduction velocity as the grid spacing increases. This grid-sensitivity requires fine mesh resolutions (typically $\sim 200\,\mu\text{m}$) to avoid numerical artifact blocks.

In CVFEM, the mass matrix is assembled by integrating the linear interpolants exactly over each control volume. This consistent mass matrix formulation reduces discretization error at larger grid spacings, allowing accurate simulations on coarser meshes. This helps mitigate the sensitivity of conduction velocity to spatial discretization, though some tuning of the conductivity tensors remains common practice to match physiological targets.

### DMP Preservation with Godunov Fluxes

Preserving the Discrete Maximum Principle (DMP) is essential to prevent unphysical oscillations in the transmembrane potential, particularly during the steep upstroke of the action potential. In standard linear finite element methods, the DMP is only guaranteed under strict geometric constraints, such as requiring that all angles of a triangular mesh are acute. On highly distorted unstructured grids, these conditions are often violated.

To ensure the DMP without geometric restrictions, positive non-linear CVFE schemes have been developed. These schemes discretize the anisotropic diffusion terms on the dual mesh using fluxes provided by the conforming finite element reconstruction on the primal mesh. The remaining terms, such as the non-linear ionic reactions, are discretized using a non-classical upwind or Godunov-type flux approximation on the dual control volumes. This upwind treatment ensures that the off-diagonal entries of the resulting system matrix remain non-positive, satisfying the conditions of an M-matrix and preventing unphysical oscillations regardless of mesh distortion.

## The Tissue-Bath Interface and Boundary Layer Modeling

### Mathematical Conditions at the Boundary

When modeling cardiac electrical activity, the heart $\Omega_H$ is typically surrounded by a surrounding conductive volume conductor $\Omega_B$ representing the blood in the cavities or the body tissues. This setup is known as the tissue-bath interface problem. Since cells do not extend past the epicardium or endocardium, the intracellular current cannot flow out of the heart, requiring a no-flux boundary condition at the interface $\partial\Omega_H$:

$$\mathbf{n} \cdot \left(\mathbf{M}_i \nabla \phi_i\right) = 0 \quad \text{on} \quad \partial\Omega_H$$ [cite: 7, 8]

Substituting $\phi_i = V_m + \phi_e$, this boundary condition is expressed in terms of the solved potentials as:

$$\mathbf{n} \cdot \left(\mathbf{M}_i \nabla V_m + \mathbf{M}_i \nabla \phi_e\right) = 0 \quad \text{on} \quad \partial\Omega_H$$ [cite: 8]

Conversely, the extracellular space is physically continuous with the surrounding bath. This continuity requires that the extracellular potential $\phi_e$ and the normal current flux match the bath potential $\phi_b$ and its associated flux across the interface:

$$\phi_e = \phi_b \quad \text{on} \quad \partial\Omega_H$$ [cite: 40]

$$\mathbf{n} \cdot \left(\mathbf{M}_e \nabla \phi_e\right) = \mathbf{n} \cdot \left(\sigma_b \nabla \phi_b\right) \quad \text{on} \quad \partial\Omega_H$$ [cite: 40]

Here, $\sigma_b$ represents the isotropic conductivity of the bath, and $\mathbf{n}$ is the outward unit normal vector pointing from the heart into the bath. Within the bath domain, the potential is governed by Laplace's equation:

$$\nabla^2 \phi_b = 0 \quad \text{in} \quad \Omega_B$$ [cite: 40]

### Analytical Resolution of the Boundary Layer

The tissue-bath interface is characterized by a narrow boundary layer where the potentials vary rapidly. This behavior is governed by the tissue length constant $\lambda$, which is typically on the scale of $\sim 1$ mm. Resolving this boundary layer numerically requires fine spatial grids near the interface.

To address this challenge, Patel and Roth developed an analytical decoupling approach. Using the linear transformations for $\phi_m$ and $\psi$, they decoupled the bidomain equations and introduced analytical auxiliary functions to resolve the rapid spatial decay near the boundary. For a tissue of semi-infinite extent, they expressed the potentials using an exponential decay function:

$$\phi_m = \phi_m^0 + A e^{-n/\lambda}$$ [cite: 9]

$$\psi = \psi^0 + B e^{-n/\lambda}$$ [cite: 9]

Here, $n$ represents the distance normal to the interface, and $\phi_m^0$ and $\psi^0$ are the background potentials. For domains of finite thickness, the simple exponential decay is replaced with hyperbolic functions (such as $\cosh(n/\lambda)$ or $\sinh(n/\lambda)$) or modified Bessel functions. These functions exhibit the same rapid decay as the exponential function but satisfy the physical boundary conditions at the opposite faces of the finite domain. Incorporating these auxiliary functions into the numerical discretization allows the boundary layer to be resolved accurately without requiring excessive grid refinement at the interface.

### Discrete Implementation: Conforming FVM vs. Decoupled Solvers

The numerical implementation of these interface conditions depends on the coupling of the solver and the mesh structure.

In conforming unstructured finite volume methods, the control volumes are constructed to align with the tissue-bath interface $\partial\Omega_H$. The interface conditions are enforced directly during the assembly of the numerical fluxes across the boundary faces. Because flux conservation is built into the finite volume formulation, these conditions are satisfied naturally without requiring ghost nodes or artificial boundary approximations.

In decoupled solvers, where the equations for $V_m$ and $\phi_e$ are solved sequentially, the intracellular boundary condition $\mathbf{n} \cdot (\mathbf{M}_i \nabla V_m + \mathbf{M}_i \nabla \phi_e) = 0$ must be split. This splitting can introduce time-lag errors and numerical instability unless small time steps are used.

Solving the bidomain equations as a fully coupled system avoids these splitting errors. In a coupled solver, the boundary conditions are integrated directly into the global linear system. Although this results in a larger system of algebraic equations at each time step, fully coupled implicit solvers converge in fewer iterations. This coupled approach can provide overall computational speedups of 50% to 80% compared to decoupled schemes, while enforcing the boundary conditions more rigorously at the interface.

| Boundary Formulation | Mesh Requirements | Interface Representation | Computational Cost | Handling of Decoupled Boundaries |
| :--- | :--- | :--- | :--- | :--- |
| **Conforming FVM** | Grid faces must align with the interface $\partial\Omega_H$ [cite: 13] | Sharp, explicit representation of the boundary | Moderate to high; requires conforming mesh generation | Natural; boundary conditions are satisfied directly in flux evaluations |
| **Smoothed Boundary (Diffuse Interface)** | No boundary tracking; uses uniform voxel meshes | Spatially diffuse interface of finite thickness | Lower mesh cost; requires higher local resolution | Handled implicitly via phase-field scaling functions |
| **Analytical Decoupling (Patel-Roth)** | Compatible with both structured and unstructured grids | Resolved analytically using auxiliary functions | Lower; reduces the need for fine grids at the boundary | Approximates the boundary layer decay analytically |
