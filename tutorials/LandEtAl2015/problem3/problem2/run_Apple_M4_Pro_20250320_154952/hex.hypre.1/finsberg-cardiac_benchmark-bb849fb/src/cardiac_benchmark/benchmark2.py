import dolfin
import logging
import numpy as np
from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Union

from . import activation_model
from . import postprocess
from . import pressure_model
from .geometry import BiVGeometry
from .material import HolzapfelOgden
from .problem import BiVProblem
from .utils import _update_parameters

logger = logging.getLogger(__name__)

dolfin.parameters["form_compiler"]["quadrature_degree"] = 6
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["optimize"] = True


def solve(
    problem: BiVProblem,
    tau: dolfin.Constant,
    activation: np.ndarray,
    lv_pressure: np.ndarray,
    rv_pressure: np.ndarray,
    plv: dolfin.Constant,
    prv: dolfin.Constant,
    time: np.ndarray,
    collector: postprocess.DataCollector,
    store_freq: int = 1,
) -> None:
    """Solve the problem for benchmark 2

    Parameters
    ----------
    problem : BiVProblem
        The problem
    tau : dolfin.Constant
        Constant in the model representing the activation
    activation : np.ndarray
        An array of activation points
    lv_pressure : np.ndarray
        An array of LV pressure points
    rv_pressure : np.ndarray
        An array of RV pressure points
    plv : dolfin.Constant
        Constant in the model representing the LV pressure
    prv : dolfin.Constant
        Constant in the model representing the RV pressure
    time : np.ndarray
        Time stamps
    collector : postprocess.DataCollector
        Datacollector used to store the results
    store_freq : int, optional
        Frequency of how often to store the results, by default 1

    Raises
    ------
    RuntimeError
        If the solver does not converge
    """
    for i, (t, a, plv_, prv_) in enumerate(
        zip(time, activation, lv_pressure, rv_pressure),
    ):
        logger.info(
            f"{i}: Solving for time {t:.3f} with tau = {a:.3f}, lvp = {plv_:.3f} and rvp = {prv_:.3f}",
        )

        tau.assign(a)
        plv.assign(plv_)
        prv.assign(prv_)
        converged = problem.solve()
        if not converged:
            raise RuntimeError

        if i % store_freq == 0:
            dolfin.info("Store solution")
            collector.store(t)


def default_parameters():
    """Default parameters for Benchmark 2"""
    return dict(
        problem_parameters=BiVProblem.default_parameters(),
        material_parameters=HolzapfelOgden.default_parameters(),
        lv_pressure_parameters=pressure_model.default_lv_parameters_benchmark2(),
        rv_pressure_parameters=pressure_model.default_rv_parameters_benchmark2(),
        activation_parameters=activation_model.default_parameters(),
        geometry_path="biv_geometry.h5",
        outpath="results_benchmark2.h5",
        T=1.0,
    )


def run(
    geometry_path: Union[str, Path] = "biv_geometry.h5",
    activation_parameters: Optional[Dict[str, float]] = None,
    lv_pressure_parameters: Optional[Dict[str, float]] = None,
    rv_pressure_parameters: Optional[Dict[str, float]] = None,
    material_parameters: Optional[Dict[str, Union[float, dolfin.Constant]]] = None,
    problem_parameters: Optional[Dict[str, Union[float, dolfin.Constant]]] = None,
    outpath: Union[str, Path] = "results_benchmark2.h5",
    T: float = 1.0,
):
    """Run benchmark 2

    Parameters
    ----------
    geometry_path : Union[str, Path], optional
        Path to the geometry, by default "biv_geometry.h5"
    activation_parameters : Optional[Dict[str, float]], optional
        Parameters for the activation model, by default None
    lv_pressure_parameters : Optional[Dict[str, float]], optional
        Parameters for the pressure model in the LV, by default None
    rv_pressure_parameters : Optional[Dict[str, float]], optional
        Parameters for the pressure model in the RV, by default None
    material_parameters : Optional[Dict[str, Union[float, dolfin.Constant]]], optional
        Parameters for the material model, by default None
    problem_parameters : Optional[Dict[str, Union[float, dolfin.Constant]]], optional
        Parameters for the problem, by default None
    outpath : Union[str, Path], optional
        Path to where to save the results, by default "results.h5"
    T : float, optional
        End time of simulation, by default 1.0

    Raises
    ------
    OSError
        If output file is not an HDF5 file
    """
    outdir = Path(outpath).parent
    outdir.mkdir(parents=True, exist_ok=True)

    geo = BiVGeometry.from_file(geometry_path)

    problem_parameters = _update_parameters(
        BiVProblem.default_parameters(),
        problem_parameters,
    )
    material_parameters = _update_parameters(
        HolzapfelOgden.default_parameters(),
        material_parameters,
    )
    lv_pressure_parameters = _update_parameters(
        pressure_model.default_lv_parameters_benchmark2(),
        lv_pressure_parameters,
    )
    rv_pressure_parameters = _update_parameters(
        pressure_model.default_rv_parameters_benchmark2(),
        rv_pressure_parameters,
    )
    activation_parameters = _update_parameters(
        activation_model.default_parameters(),
        activation_parameters,
    )

    dt = float(problem_parameters["dt"])
    tau = dolfin.Constant(0.0)
    time = np.arange(dt, T, dt)

    t_eval = time - float(problem_parameters["alpha_f"]) * dt
    lv_pressure = pressure_model.pressure_function(
        (0, T),
        t_eval=t_eval,
        parameters=lv_pressure_parameters,
    )
    rv_pressure = pressure_model.pressure_function(
        (0, T),
        t_eval=t_eval,
        parameters=rv_pressure_parameters,
    )
    activation = activation_model.activation_function(
        (0, T),
        t_eval=t_eval,
        parameters=activation_parameters,
    )

    if dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
        np.save(
            outdir / "pressure_model.npy",
            {
                "time": time,
                "activation": activation,
                "lv_pressure": lv_pressure,
                "rv_pressure": rv_pressure,
                "lv_pressure_parameters": lv_pressure_parameters,
                "rv_pressure_parameters": rv_pressure_parameters,
                "activation_parameters": activation_parameters,
            },
        )

        postprocess.plot_activation_pressure_function(
            t=time,
            activation=activation,
            lv_pressure=lv_pressure,
            rv_pressure=rv_pressure,
            outdir=outdir,
        )

    plv = dolfin.Constant(0.0)
    prv = dolfin.Constant(0.0)
    problem_parameters["plv"] = plv
    problem_parameters["prv"] = prv

    material = HolzapfelOgden(
        f0=geo.f0,
        s0=geo.s0,
        tau=tau,
        parameters=material_parameters,
    )
    problem = BiVProblem(
        geometry=geo,
        material=material,
        parameters=problem_parameters,
    )

    problem.solve()

    result_filepath = Path(outpath)
    if result_filepath.suffix != ".h5":
        msg = f"Expected output path to be to type HDF5 with suffix .h5, got {result_filepath.suffix}"
        raise OSError(msg)
    result_filepath.parent.mkdir(exist_ok=True)

    collector = postprocess.DataCollector(
        result_filepath,
        problem=problem,
        pressure_parameters={
            "lv": lv_pressure_parameters,
            "rv": rv_pressure_parameters,
        },
        activation_parameters=activation_parameters,
    )

    solve(
        problem=problem,
        tau=tau,
        activation=activation,
        lv_pressure=lv_pressure,
        rv_pressure=rv_pressure,
        plv=plv,
        prv=prv,
        time=time,
        collector=collector,
        store_freq=1,
    )
