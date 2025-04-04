# Time integration schemes

## The Newmark $\beta$ scheme

A popular method for time integration is the Newmark $\beta$ scheme,
which gives $a,v$ as

$$
v_{i+1} = v_i + (1-\gamma) \Delta t ~ a_i + \gamma \Delta t ~ a_{i+1}
$$ (N_v)

$$
a_{i+1} = \frac{u_{i+1} - (u_i + \Delta t ~ v_i + (0.5 - \beta) \Delta t^2 ~ a_i)}{\beta \Delta t^2}
$$ (N_a)

Inserting these expressions into {eq}`weak1` and assuming $u_i,v_i,a_i$ known, we obtain the time discrete weak form

$$
\begin{aligned}
\int_{\Omega} \rho a_{i+1} \cdot w \, \mathrm{d}X + \int_{\Omega} P(u_{i+1},v_{i+1}):Grad(w) \, \mathrm{d}X -\int_{\Gamma_{\rm endo}} p I JF_{i+1}^{-T}N \cdot w \, \mathrm{d}S  \\
+\int_{\Gamma_{\rm epi}} \big(\alpha_{\rm epi} u_{i+1} \cdot N + \beta_{\rm epi} v_{i+1} \cdot N \big) w \cdot N \, \mathrm{d}S  \\
+\int_{\Gamma_{\rm top}} \alpha_{\rm top} u_{i+1} \cdot w + \beta_{\rm top} v_{i+1} \cdot w \, \mathrm{d}S = 0 \quad \forall w \in H^1(\Omega),
\end{aligned}
$$ (weak2)

where $F_{i+1} = I + \nabla u_{i+1}$ and $v_i,a_i$ are given by
{eq}`N_v` - {eq}`N_a`.

```{note}
We could insert these expressions into {eq}`weak2` to eliminate $v_i,a_i$
and yield a non-linear weak form with $u_{i+1}$ as the only unknown. However, the
resulting form becomes fairly complex, and we avoid this task by implementing {eq}`weak2` and
{eq}`N_v` - {eq}`N_a`. directly as UFL forms in the code, and leave the algebra to UFL.
```

A common choice is $\gamma=1/2, \beta=1/4$, which yields
the *average constant acceleration (middle point rule)* method, while $\gamma=1/2, \beta=1/6$
yields the *linear acceleration method* where the acceleration is linearly varying between
$t$ and $t+\Delta t$. Both methods are unconditionally stable and second order accurate, but
they introduce no numerical damping and can suffer from spurious oscillations. Choosing $\gamma > 1/2$
introduces energy dissipation which can avoid such oscillations, but reduces the
accuracy to first order. The method is unconditionally stable for $1/2 \leq \gamma \leq 2\beta$.


## The generalized $\alpha$ method

Several methods have been derived with the purpose of avoiding the non-physical
oscillations sometimes seen in the Newmark method for $\gamma =1/2$, while retaining
the second order convergence in $\Delta t$. One class of such methods is called
the generalized $\alpha$ or G-$\alpha$ methods, which introduce additional parameters
$\alpha_f$ and $\alpha_m$. The methods use the Newmark approximations in
{eq}`N_v` - {eq}`N_a` to approximate $v$ and $a$, but evaluate the terms of the
weak form at times $t_{i+1}-\Delta t\alpha_m$ and  $t_{i+1}-\Delta t\alpha_f$.
Specifically, the inertia term is evaluated at $t_{i+1}-\Delta t\alpha_m$, and the
other terms at $t_{i+1}-\Delta t\alpha_f$. The weak form becomes

$$
\begin{aligned}
\int_{\Omega} \rho a_{i+1-\alpha_m} \cdot w \, \mathrm{d}X + \int_{\Omega} P_{i+1-\alpha_f}:Grad(w) \, \mathrm{d}X -\int_{\Gamma_{\rm endo}} p I JF_{i+1-\alpha_f}^{-T}N \cdot w \, dS  \\
+\int_{\Gamma_{\rm epi}} \big(\alpha_{\rm epi} u_{i+1-\alpha_f} \cdot N + \beta_{\rm epi} v_{i+1-\alpha_f} \cdot N \big) w \cdot N \, dS  \\
+\int_{\Gamma_{\rm top}} \alpha_{\rm top} u_{i+1-\alpha_f} \cdot w + \beta_{\rm top} v_{i+1-\alpha_f} \cdot w \, dS = 0 \quad \forall w \in H^1(\Omega),
\end{aligned}
$$ (weak3)

with

$$
\begin{align*}
  u_{i+1-\alpha_f} &= (1-\alpha_f)u_{i+1}-\alpha_f u_i, \\
  v_{i+1-\alpha_f} &= (1-\alpha_f)v_{i+1}-\alpha_f v_i, \\
  a_{i+1-\alpha_m} &= (1-\alpha_m)a_{i+1}-\alpha_m a_i,
\end{align*}
$$

$v_{i+1},a_{i+1}$ given by {eq}`N_v` - {eq}`N_a`, and

$$
\begin{align*}
F_{i+1-\alpha_f} &= I + \nabla u_{i+1-\alpha_f}, \\
P_{i+1-\alpha_f} &= P(u_{i+1-\alpha_f}, v_{i+1-\alpha_f}).
\end{align*}
$$

The only difference between {eq}`weak3` and {eq}`weak2` is the time point at which the
terms of the equation are evaluated.

Different choices of the four parameters $\alpha_m, \alpha_f, \beta, \gamma$ yield
methods with different accuracy and stability properties. Tables 1--3 in {cite}`erlicher2002analysis`
provides an overview of parameter choices for methods in the literature,
as well as conditions for stability and convergence.
We have used the choice $\alpha_m =0.2, \alpha_f=0.4$, and

$$
\begin{align*}
  \gamma &= 1/2 + \alpha_f-\alpha_m ,\\
  \beta &= \frac{(\gamma + 1/2)^2}{4} .
\end{align*}
$$

For this choice the solver converges through the time interval of interest, and
the convergence is second order.


### Alternative formulations of the G-$\alpha$ method.


- In addition to the choice of the four parameters, different choices
  can be made in how the quantities $F_{i+1-\alpha_f}$ and $P_{i+1-\alpha_f}$
  are approximated. Alternative choices include

  $$
  \begin{align*}
   F_{i+1-\alpha_f} = (1-\alpha_f) (I + \nabla u_{i+1}) + \alpha_f(I + \nabla u_{i}), \\
   P_{i+1-\alpha_f} = (1-\alpha_f) P(u_{i+1}, v_{i+1}) + \alpha_f P (I + \nabla u_{i}, v_{i}) ,
 \end{align*}
 $$

 and other options, see, for instance, {cite}`erlicher2002analysis` and references therein.
 These alternatives have not yet been explored.

- So far, we have implemented the weak form {eq}`weak3` and the relations
 {eq}`N_v` - {eq}`N_a` directly, and relied on UFL to handle the algebra to turn it
 into a problem with $u_{i+1}$ as the only unknown. This seems to work well, but
 there are alternatives. For instance, Eq. (14) in {cite}`erlicher2002analysis` is
 a relatively compact equation for $u_{i+1}$ based on the G-$\alpha$ method. It should
 be possible to derive a weak form similar to this equation for our problem, which could be
 solved and compared with the UFL approach.



### The Newmark $\beta$ scheme
We consider first a general dynamic equation on the form

$$
M\ddot{u} + C\dot{u} + f^{\rm int}(u) = f_{\rm ext},
$$

where $M$ is the mass matrix, $C$ is the damping matrix, and $f^{\rm int}, f^{\rm ext}$
are internal and external forces. Note that the mechanics equation above can also
be written on this general form, since the stress tensor $P$ includes a viscous
component which is a linear function of the velocity $\dot{u}$.

The Newmark $\beta$ method is then commonly written as

$$
\begin{align}
v_{i+1} &= v_i + (1-\gamma) \Delta t ~ a_i + \gamma \Delta t ~ a_{i+1}, \\
a_{i+1} &= \frac{u_{i+1} - (u_i + \Delta t ~ v_i + (0.5 - \beta) \Delta t^2 ~ a_i)}{\beta \Delta t^2}, \\
Ma_{i+1} &+ C v_{i+1} + f^{\rm int}(u_{i+1}) = f^{\rm ext}_{i+1} .
\end{align}
$$
