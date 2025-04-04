import math
import pprint
from typing import Dict
from typing import Optional
from typing import Tuple

import numpy as np
import scipy.integrate


def default_parameters_benchmark1() -> Dict[str, float]:
    r"""Default parameters for the pressure model in benchmark 1

    Returns
    -------
    Dict[str, float]
        Default parameters

    Notes
    -----
    The default parameters are

    .. math::
        t_{\mathrm{sys} - \mathrm{pre}} &= 0.17 \\
        t_{\mathrm{dias} - \mathrm{pre}} &= 0.484 \\
        \gamma &= 0.005 \\
        a_{\mathrm{max}} &= 5.0 \\
        a_{\mathrm{min}} &= -30.0 \\
        \alpha_{\mathrm{pre}} &= 5.0 \\
        \alpha_{\mathrm{mid}} &= 1.0 \\
        \sigma_{\mathrm{pre}} &= 7000.0 \\
        \sigma_{\mathrm{mid}} &= 16000.0 \\
    """
    return dict(
        t_sys_pre=0.17,
        t_dias_pre=0.484,
        gamma=0.005,
        a_max=5.0,
        a_min=-30.0,
        alpha_pre=5.0,
        alpha_mid=1.0,
        sigma_pre=7000.0,
        sigma_mid=16000.0,
    )


def default_lv_parameters_benchmark2() -> Dict[str, float]:
    r"""Default parameters for the LV pressure model in benchmark 2

    Returns
    -------
    Dict[str, float]
        Default parameters

    Notes
    -----
    The default parameters are

    .. math::
        t_{\mathrm{sys} - \mathrm{pre}} &= 0.17 \\
        t_{\mathrm{dias} - \mathrm{pre}} &= 0.484 \\
        \gamma &= 0.005 \\
        a_{\mathrm{max}} &= 5.0 \\
        a_{\mathrm{min}} &= -30.0 \\
        \alpha_{\mathrm{pre}} &= 5.0 \\
        \alpha_{\mathrm{mid}} &= 15.0 \\
        \sigma_{\mathrm{pre}} &= 12000.0 \\
        \sigma_{\mathrm{mid}} &= 16000.0 \\
    """
    return dict(
        t_sys_pre=0.17,
        t_dias_pre=0.484,
        gamma=0.005,
        a_max=5.0,
        a_min=-30.0,
        alpha_pre=5.0,
        alpha_mid=15.0,
        sigma_pre=12000.0,
        sigma_mid=16000.0,
    )


def default_rv_parameters_benchmark2():
    r"""Default parameters for the RV pressure model in benchmark 2

    Returns
    -------
    Dict[str, float]
        Default parameters

    Notes
    -----
    The default parameters are

    .. math::
        t_{\mathrm{sys} - \mathrm{pre}} &= 0.17 \\
        t_{\mathrm{dias} - \mathrm{pre}} &= 0.484 \\
        \gamma &= 0.005 \\
        a_{\mathrm{max}} &= 5.0 \\
        a_{\mathrm{min}} &= -30.0 \\
        \alpha_{\mathrm{pre}} &= 5.0 \\
        \alpha_{\mathrm{mid}} &= 10.0 \\
        \sigma_{\mathrm{pre}} &= 3000.0 \\
        \sigma_{\mathrm{mid}} &= 4000.0 \\
    """
    return dict(
        t_sys_pre=0.17,
        t_dias_pre=0.484,
        gamma=0.005,
        a_max=5.0,
        a_min=-30.0,
        alpha_pre=1.0,
        alpha_mid=10.0,
        sigma_pre=3000.0,
        sigma_mid=4000.0,
    )


def pressure_function(
    t_span: Tuple[float, float],
    parameters: Dict[str, float],
    t_eval: Optional[np.ndarray] = None,
) -> np.ndarray:
    r"""Time-dependent pressure derived from the Bestel model [3]_.

    Parameters
    ----------
    t_span : Tuple[float, float]
        A tuple representing start and end of time
    parameters : Dict[str, float]
        Parameters used in the model, see :func:`default_parameters`
    t_eval : Optional[np.ndarray], optional
        Time points to evaluate the solution, by default None.
        If not provided, the default points from `scipy.integrate.solve_ivp`
        will be used

    Returns
    -------
    np.ndarray
        An array of pressure points

    Notes
    -----
    We consider a time-dependent pressure derived from the Bestel model.
    The solution :math:`p = p(t)` is characterized as solution to the evolution equation

    .. math::
        \dot{p}(t) = -|b(t)|p(t) + \sigma_{\mathrm{mid}}|b(t)|_+
        + \sigma_{\mathrm{pre}}|g_{\mathrm{pre}}(t)|

    being b(\cdot) the activation function described below:

    .. math::
        b(t) =& a_{\mathrm{pre}}(t) + \alpha_{\mathrm{pre}}g_{\mathrm{pre}}(t)
        + \alpha_{\mathrm{mid}} \\
        a_{\mathrm{pre}}(t) :=& \alpha_{\mathrm{max}} \cdot f_{\mathrm{pre}}(t)
        + \alpha_{\mathrm{min}} \cdot (1 - f_{\mathrm{pre}}(t)) \\
        f_{\mathrm{pre}}(t) =& S^+(t - t_{\mathrm{sys}-\mathrm{pre}}) \cdot
         S^-(t  t_{\mathrm{dias} - \mathrm{pre}}) \\
        g_{\mathrm{pre}}(t) =& S^-(t - t_{\mathrm{dias} - \mathrm{pre}})

    with :math:`S^{\pm}` given by

    .. math::
        S^{\pm}(\Delta t) = \frac{1}{2}(1 \pm \mathrm{tanh}(\frac{\Delta t}{\gamma}))


    """
    print(f"Solving pressure model with parameters: {pprint.pformat(parameters)}")

    f = (
        lambda t: 0.25
        * (1 + math.tanh((t - parameters["t_sys_pre"]) / parameters["gamma"]))
        * (1 - math.tanh((t - parameters["t_dias_pre"]) / parameters["gamma"]))
    )
    a = lambda t: parameters["a_max"] * f(t) + parameters["a_min"] * (1 - f(t))

    f_pre = lambda t: 0.5 * (
        1 - math.tanh((t - parameters["t_dias_pre"]) / parameters["gamma"])
    )
    b = lambda t: a(t) + parameters["alpha_pre"] * f_pre(t) + parameters["alpha_mid"]

    def rhs(t, p):
        return (
            -abs(b(t)) * p
            + parameters["sigma_mid"] * max(b(t), 0)
            + parameters["sigma_pre"] * max(f_pre(t), 0)
        )

    res = scipy.integrate.solve_ivp(
        rhs,
        t_span,
        [0.0],
        t_eval=t_eval,
        method="Radau",
    )
    return res.y.squeeze()
