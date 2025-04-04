r"""
Implementation of the mechanics problem

For time integration we employ the generalized :math:`\alpha`-method [1]_.

    .. [1] Silvano Erlicher, Luca Bonaventura, Oreste Bursi.
        The analysis of the Generalized-alpha method for
        non-linear dynamic problems. Computational Mechanics,
        Springer Verlag, 2002, 28, pp.83-104, doi:10.1007/s00466-001-0273-z
"""
import abc
import typing

import dolfin

try:
    import ufl_legacy as ufl
except ImportError:
    import ufl

from .geometry import HeartGeometry
from .material import HolzapfelOgden
from .solver import NonlinearProblem
from .solver import NonlinearSolver

T = typing.TypeVar("T", dolfin.Function, dolfin.Vector)


def interpolate(x0: T, x1: T, alpha: float):
    r"""Interpolate beteween :math:`x_0` and :math:`x_1`
    to find `math:`x_{1-\alpha}`

    Parameters
    ----------
    x0 : T
        First point
    x1 : T
        Second point
    alpha : float
        Amount of interpolate

    Returns
    -------
    T
        `math:`x_{1-\alpha}`
    """
    return alpha * x0 + (1 - alpha) * x1


class Problem(abc.ABC):
    def __init__(
        self,
        geometry: HeartGeometry,
        material: HolzapfelOgden,
        parameters: typing.Optional[
            typing.Dict[str, typing.Union[dolfin.Constant, str]]
        ] = None,
        solver_parameters=None,
    ) -> None:
        """Constructor

        Parameters
        ----------
        geometry : HeartGeometry
            The geometry
        material : HolzapfelOgden
            The material
        parameters : typing.Dict[str, Union[dolfin.Constant, str], optional
            Problem parameters, by default None. See
            `Problem.default_parameters`
        """
        self.geometry = geometry
        self.material = material

        parameters = parameters or {}

        self.parameters = type(self).default_parameters()
        self.parameters.update(parameters)

        self.solver_parameters = NonlinearSolver.default_solver_parameters()
        if solver_parameters is not None:
            self.solver_parameters.update(**solver_parameters)
        self._init_spaces()
        self._init_forms()

    @staticmethod
    @abc.abstractmethod
    def default_parameters():
        ...

    def _init_spaces(self):
        """Initialize function spaces"""
        mesh = self.geometry.mesh

        family, degree = self.parameters["function_space"].split("_")

        element = dolfin.VectorElement(family, mesh.ufl_cell(), int(degree))
        self.u_space = dolfin.FunctionSpace(mesh, element)
        self.u = dolfin.Function(self.u_space)
        self.u_test = dolfin.TestFunction(self.u_space)
        self.du = dolfin.TrialFunction(self.u_space)

        self.u_old = dolfin.Function(self.u_space)
        self.v_old = dolfin.Function(self.u_space)
        self.a_old = dolfin.Function(self.u_space)

    def _acceleration_form(self, a: dolfin.Function, w: dolfin.TestFunction):
        return ufl.inner(self.parameters["rho"] * a, w) * ufl.dx

    def _first_piola(self, F: ufl.Coefficient, v: dolfin.Function):
        F_dot = ufl.grad(v)
        l = F_dot * ufl.inv(F)  # Holzapfel eq: 2.139
        d = 0.5 * (l + l.T)  # Holzapfel 2.146
        E_dot = ufl.variable(F.T * d * F)  # Holzapfel 2.163

        return ufl.diff(self.material.strain_energy(F), F) + F * ufl.diff(
            self.material.W_visco(E_dot),
            E_dot,
        )

    def _form(self, u: dolfin.Function, v: dolfin.Function, w: dolfin.TestFunction):
        F = ufl.variable(ufl.grad(u) + ufl.Identity(3))
        P = self._first_piola(F, v)
        epi = ufl.dot(self.parameters["alpha_epi"] * u, self.N) + ufl.dot(
            self.parameters["beta_epi"] * v,
            self.N,
        )
        top = self.parameters["alpha_top"] * u + self.parameters["beta_top"] * v

        return (
            ufl.inner(P, ufl.grad(w)) * ufl.dx
            + self._pressure_term(F, w)
            + ufl.inner(epi * w, self.N) * self.ds(self.epi)
            + ufl.inner(top, w) * self.ds(self.top)
        )

    @abc.abstractmethod
    def _pressure_term(self, F, w):
        ...

    def v(
        self,
        a: T,
        v_old: T,
        a_old: T,
    ) -> T:
        r"""
        Velocity computed using the generalized
        :math:`alpha`-method

        .. math::
            v_{i+1} = v_i + (1-\gamma) \Delta t a_i + \gamma \Delta t a_{i+1}

        Parameters
        ----------
        a : T
            Current acceleration
        v_old : T
            Previous velocity
        a_old: T
            Previous acceleration
        Returns
        -------
        T
            The current velocity
        """
        dt = self.parameters["dt"]
        return v_old + (1 - self._gamma) * dt * a_old + self._gamma * dt * a

    def a(
        self,
        u,
        u_old,
        v_old,
        a_old,
    ) -> typing.Union[dolfin.Function, dolfin.Vector]:
        r"""
        Acceleration computed using the generalized
        :math:`alpha`-method

        .. math::
            a_{i+1} = \frac{u_{i+1} - (u_i + \Delta t v_i + (0.5 - \beta) \Delta t^2 a_i)}{\beta \Delta t^2}

        Parameters
        ----------
        u : T
            Current displacement
        u_old : T
            Previous displacement
        v_old : T
            Previous velocity
        a_old: T
            Previous acceleration
        Returns
        -------
        T
            The current acceleration
        """
        dt = self.parameters["dt"]
        dt2 = dt**2
        beta = self._beta
        return (u - (u_old + dt * v_old + (0.5 - beta) * dt2 * a_old)) / (beta * dt2)

    def _update_fields(self) -> None:
        """Update old values of displacement, velocity
        and acceleration
        """
        a = self.a(
            u=self.u.vector(),
            u_old=self.u_old.vector(),
            v_old=self.v_old.vector(),
            a_old=self.a_old.vector(),
        )
        v = self.v(a=a, v_old=self.v_old.vector(), a_old=self.a_old.vector())

        self.a_old.vector()[:] = a
        self.v_old.vector()[:] = v
        self.u_old.vector()[:] = self.u.vector()

    @property
    def ds(self):
        """Surface measure"""
        return ufl.ds(domain=self.geometry.mesh, subdomain_data=self.geometry.ffun)

    @property
    def epi(self):
        """Marker for the epicardium"""
        return self.geometry.markers["EPI"][0]

    @property
    def top(self):
        """Marker for the top or base"""
        return self.geometry.markers["BASE"][0]

    @property
    def N(self):
        """Facet Noraml"""
        return ufl.FacetNormal(self.geometry.mesh)

    def _init_forms(self) -> None:
        """Initialize ufl forms"""
        w = self.u_test

        # Markers
        if self.geometry.markers is None:
            raise RuntimeError("Missing markers in geometry")

        alpha_m = self.parameters["alpha_m"]
        alpha_f = self.parameters["alpha_f"]

        a_new = self.a(u=self.u, u_old=self.u_old, v_old=self.v_old, a_old=self.a_old)
        v_new = self.v(a=a_new, v_old=self.v_old, a_old=self.a_old)

        virtual_work = self._acceleration_form(
            interpolate(self.a_old, a_new, alpha_m),
            w,
        ) + self._form(
            interpolate(self.u_old, self.u, alpha_f),
            interpolate(self.v_old, v_new, alpha_f),
            w,
        )
        jacobian = ufl.derivative(
            virtual_work,
            self.u,
            self.du,
        )

        self._problem = NonlinearProblem(J=jacobian, F=virtual_work, bcs=[])
        self.solver = NonlinearSolver(
            self._problem,
            self.u,
            parameters=self.solver_parameters,
        )

    def von_Mises(self) -> ufl.Coefficient:
        r"""Compute the von Mises stress tensor :math`\sigma_v`, with

        .. math::

            \sigma_v^2 = \frac{1}{2} \left(
                (\mathrm{T}_{11} - \mathrm{T}_{22})^2 +
                (\mathrm{T}_{22} - \mathrm{T}_{33})^2 +
                (\mathrm{T}_{33} - \mathrm{T}_{11})^2 +
            \right) - 3 \left(
                \mathrm{T}_{12} + \mathrm{T}_{23} + \mathrm{T}_{31}
            \right)

        Returns
        -------
        ufl.Coefficient
            The von Mises stress tensor
        """
        u = self.u
        a = self.a(u=self.u, u_old=self.u_old, v_old=self.v_old, a_old=self.a_old)
        v = self.v(a=a, v_old=self.v_old, a_old=self.a_old)

        F = ufl.variable(ufl.grad(u) + ufl.Identity(3))
        J = ufl.det(F)
        P = self._first_piola(F, v)

        # Cauchy
        T = pow(J, -1.0) * P * F.T
        von_Mises_squared = 0.5 * (
            (T[0, 0] - T[1, 1]) ** 2
            + (T[1, 1] - T[2, 2]) ** 2
            + (T[2, 2] - T[0, 0]) ** 2
        ) + 3 * (T[0, 1] + T[1, 2] + T[2, 0])

        return ufl.sqrt(abs(von_Mises_squared))

    @property
    def _gamma(self) -> dolfin.Constant:
        """Parameter in the generalized alpha-method"""
        return dolfin.Constant(
            0.5 + self.parameters["alpha_f"] - self.parameters["alpha_m"],
        )

    @property
    def _beta(self) -> dolfin.Constant:
        """Parameter in the generalized alpha-method"""
        return dolfin.Constant((self._gamma + 0.5) ** 2 / 4.0)

    def solve(self) -> bool:
        """Solve the system"""
        _, conv = self.solver.solve()

        if not conv:
            self.u.assign(self.u_old)
            self._init_forms()
            return False

        self._update_fields()
        return True


class LVProblem(Problem):
    @property
    def endo(self):
        """Marker for the endocardium"""
        return self.geometry.markers["ENDO"][0]

    def _pressure_term(self, F, w):
        return ufl.inner(
            self.parameters["p"] * ufl.det(F) * ufl.inv(F).T * self.N,
            w,
        ) * self.ds(
            self.endo,
        )

    @staticmethod
    def default_parameters() -> typing.Dict[str, dolfin.Constant]:
        return dict(
            alpha_top=dolfin.Constant(1e5),
            alpha_epi=dolfin.Constant(1e8),
            beta_top=dolfin.Constant(5e3),
            beta_epi=dolfin.Constant(5e3),
            p=dolfin.Constant(0.0),
            rho=dolfin.Constant(1e3),
            dt=dolfin.Constant(1e-3),
            alpha_m=dolfin.Constant(0.2),
            alpha_f=dolfin.Constant(0.4),
            function_space="P_2",
        )


class BiVProblem(Problem):
    @property
    def endo_lv(self):
        """Marker for the endocardium"""
        return self.geometry.markers["ENDO_LV"][0]

    @property
    def endo_rv(self):
        """Marker for the endocardium"""
        return self.geometry.markers["ENDO_RV"][0]

    def _pressure_term(self, F, w):
        return ufl.inner(
            self.parameters["plv"] * ufl.det(F) * ufl.inv(F).T * self.N,
            w,
        ) * self.ds(self.endo_lv) + ufl.inner(
            self.parameters["prv"] * ufl.det(F) * ufl.inv(F).T * self.N,
            w,
        ) * self.ds(
            self.endo_rv,
        )

    @staticmethod
    def default_parameters() -> typing.Dict[str, dolfin.Constant]:
        return dict(
            alpha_top=dolfin.Constant(1e6),
            alpha_epi=dolfin.Constant(1e8),
            beta_top=dolfin.Constant(5e3),
            beta_epi=dolfin.Constant(5e3),
            plv=dolfin.Constant(0.0),
            prv=dolfin.Constant(0.0),
            rho=dolfin.Constant(1e3),
            dt=dolfin.Constant(1e-3),
            alpha_m=dolfin.Constant(0.2),
            alpha_f=dolfin.Constant(0.4),
            function_space="P_2",
        )
