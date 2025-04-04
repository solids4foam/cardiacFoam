
# Computing the cavity volume
One of the quantities reported for the benchmark problem is the cavity volume.
The simplest way to compute this quantity is by an integral over the cavity, i.e.
by

$$
\begin{align}
V &= \int_\Omega \mathrm{d}x = \int_{\Omega_0} J \mathrm{d}X,
\end{align}
$$ (vol0)

where $J=\det(F)$ and $F$ is the deformation gradient, or,
equivalently, the Jacobian of the coordinate transformation from $X$
to $x$. This is the approach described in the benchmark description.
The $F$ and $J$ fielmathrm{d}s inside the cavity are computed from a displacement field
constructed by harmonic lifting.

Our approach for computing the volume does not require a mesh of the cavity,
but instead applies the divergence theorem to turn the volume integral into
a surface integral over the endocardial surface. The divergence theorem states
that

$$
\int_\Omega \nabla\cdot f \, \mathrm{d}x = \int_{\partial\Omega} f\cdot \mathbf{n} \, \mathrm{d}s,
$$

holmathrm{d}s for any smooth function $f$ on a volume $\Omega$, with $\mathbf{n}$ being the
outward surface normal of $\partial\Omega$. If we introduce an
arbitrary vector-valued function $\mathbf{g}$ satisfying $\nabla\cdot \mathbf{g} = 1$, we have

$$
V =\int_{\Omega} \mathrm{d}x = \int_{\Omega}\nabla\cdot \mathbf{g} \, \mathrm{d}x
= \int_{\partial\Omega}\mathbf{g}\cdot\mathbf{n} \, \mathrm{d}s .
$$

There are many simple choices for the function $\mathbf{g}$, for instance $\mathbf{g}=(x_1,0,0)$ or $\mathbf{g}=\mathbf{x}/3$.

One limitation of this approach is that the surface integral is performed over the
entire boundary $\partial\Omega$, i.e., a closed surface, while the cavity volumes
of interest are usually open at the base. However, if we know the characteristics
of the geometry it is often possible to choose the function $g$ so that the
contribution from the missing part of the surface is zero, and in this case the
volume calculation is accurate. For the benchmark problem, the open
part of the surface is a plane with surface normal $(1,0,0)$. The plane moves
in the $x_1$ direction, but the surface normal stays constant through the
deformation. We can choose, for instance, $\mathbf{g}=(0,x_2,0)$ to ensure that the
contribution from this part of the surface is always zero. We then have

$$
\begin{align}
V &= \int_{\partial\Omega_{\rm endo}}(0,x_2,0) \cdot\mathbf{n} \, \mathrm{d}s ,
\end{align}
$$ (vol1)

where $\partial\Omega_{\rm endo}$ is the deformed configuration of the (open)
endocardial surface.


We want to convert the surface integral in {eq}`vol1` to an
integral in the reference configuration. We have

$$
\begin{align}
\mathbf{u}(\mathbf{X},t) &= \mathbf{x}(\mathbf{X},t)-\mathbf{X},
\end{align}
$$ (displ)

where $\mathbf{X}$ is the position vector in the undeformed (material) configuration.
From Nanson's formula we have that

$$
\begin{align}
\mathbf{n} \, \mathrm{d}s &= J\mathbf{F}^{-T}\mathbf{\eta} \, \mathrm{d}S,
\label{nanson}
\end{align}
$$ (nanson)

where $\mathbf{\eta}$ and $\mathrm{d}S$ are the surface normal and area of the undeformed surface.
Using {eq}`displ`-{eq}`nanson`, we can convert {eq}`vol1` to

$$
\begin{align}
V &=  \int_{\partial\Omega_{\rm endo,\, 0}}(u_2+X_2)\cdot
J\mathbf{F}^{-T}\mathbf{\eta} \, \mathrm{d}S ,
\end{align}
$$

where $\partial\Omega_{\rm endo,\, 0}$ is the undeformed configuration
of the endocardial surface.
