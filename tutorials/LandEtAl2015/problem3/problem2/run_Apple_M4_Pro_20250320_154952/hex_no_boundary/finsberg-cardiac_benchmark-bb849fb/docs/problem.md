# Problem formulation

We want to solve a dynamic version of the cardiac mechanics equations,
including the inertia term that is usually neglected. The boundary
conditions are as follows:

-   a standard pressure term on the endocardium ($\Gamma_{\rm endo}$);

-   normal forces on the epicardium ($\Gamma_{\rm epi}$) are a
    combination of a linear spring and linear friction;

-   zero tangent force on the epicardium; and

-   the traction (tangent + normal) on the base ($\Gamma_{\rm top}$) is
    the sum of a linear spring and linear friction.

In the benchmark description, the equations are written as

$$
\begin{aligned}
  \rho \ddot{u} - \nabla\cdot(J\sigma F^{-T}) &= 0, \mbox{ in } \Omega ,\\
  \sigma J F^{-T}N &= pJF^{-T}N, \mbox{ on } \Gamma_{\rm endo}, \\
  \sigma JF^{-T}N\cdot N + \alpha_{\rm epi}u\cdot N + \beta_{\rm epi}\dot{u}\cdot N &= 0, \mbox{ on }  \Gamma_{\rm epi}, \\
  \sigma JF^{-T}N\times N &=0, \mbox{ on }  \Gamma_{\rm epi}, \\
  \sigma JF^{-T}N + \alpha_{\rm top}u + \beta_{\rm top}\dot{u} &= 0, \mbox{ on } \Gamma_{\rm top},
\end{aligned}
$$

where we have used the notation $\sigma$ for the Cauchy
stress to be more in line with standard notation. In terms of the first Piola-Kirchhoff stress, given by $P = J\sigma F^{-T}$, the problem reads


$$
\begin{aligned}
  \rho \ddot{u} - \nabla\cdot P &= 0, \mbox{ in } \Omega ,\\
\end{aligned}
$$ (dyn_eq0)

$$
\begin{aligned}
  PN &= pJF^{-T}N, \mbox{ on } \Gamma_{\rm endo}, \\
  PN\cdot N + \alpha_{\rm epi}u\cdot N + \beta_{\rm epi}\dot{u}\cdot N &= 0, \mbox{ on }  \Gamma_{\rm epi}, \\
  PN\times N &=0, \mbox{ on }  \Gamma_{\rm epi}, \\
  PN + \alpha_{\rm top}u + \beta_{\rm top}\dot{u} &= 0, \mbox{ on } \Gamma_{\rm top}.
\end{aligned}
$$ (dyn_bc3)

If we introduce the notation $a= \ddot{u}, v=\dot{u}$
for the acceleration and velocity, respectively, a weak form
of {eq}`dyn_eq0` - {eq}`dyn_bc3` can be written as

$$
\begin{aligned}
\int_{\Omega} \rho a \cdot w \, \mathop{}\!\mathrm{d}{}X+ \int_{\Omega} P:Grad(w) \, \mathop{}\!\mathrm{d}{}X-\int_{\Gamma_{\rm endo}} p I JF^{-T}N \cdot w \, \mathop{}\!\mathrm{d}{}S\\
+\int_{\Gamma_{\rm epi}} \big(\alpha_{\rm epi} u \cdot N + \beta_{\rm epi} v \cdot N \big) w \cdot N \, \mathop{}\!\mathrm{d}{}S\\
+\int_{\Gamma_{\rm top}} \alpha_{\rm top} u \cdot w + \beta_{\rm top} v \cdot w \, \mathop{}\!\mathrm{d}{}S= 0 \quad \forall w \in H^1(\Omega).
\end{aligned}
$$ (weak1)

In order to integrate {eq}`weak1` in time, we need to express $a$ and $v$ in
terms of the displacement $u$. This can be done by numerous methods, a few of
which will be discussed below.
