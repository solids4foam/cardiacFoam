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
from .geometry import LVGeometry
from .material import HolzapfelOgden
from .problem import LVProblem
from .utils import _update_parameters

HERE = Path(__file__).absolute().parent
logger = logging.getLogger(__name__)

dolfin.parameters["form_compiler"]["quadrature_degree"] = 6
dolfin.parameters["form_compiler"]["cpp_optimize"] = True
dolfin.parameters["form_compiler"]["representation"] = "uflacs"
dolfin.parameters["form_compiler"]["optimize"] = True


def solve(
    problem: LVProblem,
    tau: dolfin.Constant,
    activation: np.ndarray,
    pressure: np.ndarray,
    p: dolfin.Constant,
    time: np.ndarray,
    collector: postprocess.DataCollector,
    store_freq: int = 1,
) -> None:
    """Solve the problem for benchmark 1

    Parameters
    ----------
    problem : LVProblem
        The problem
    tau : dolfin.Constant
        Constant in the model representing the activation
    activation : np.ndarray
        An array of activation points
    pressure : np.ndarray
        Constant in the model representing the pressure
    p : dolfin.Constant
        An array of pressure points
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
    for i, (t, a, p_) in enumerate(zip(time, activation, pressure)):
        logger.info(f"{i}: Solving for time {t:.3f} with tau = {a} and pressure = {p_}")

        tau.assign(a)
        p.assign(p_)
        converged = problem.solve()
        if not converged:
            raise RuntimeError

        if i % store_freq == 0:
            logger.info("Store solution")
            collector.store(t)


def get_geometry(
    path: Path,
    mesh_parameters: Optional[Dict[str, float]] = None,
    fiber_parameters: Optional[Dict[str, Union[float, str]]] = None,
) -> LVGeometry:
    """Get the LV geometry from a path. If the file does not
    exist, generate a new geometry with the given parameters
    and save it to the path.

    Parameters
    ----------
    path : Path
        The path to the geometry
    mesh_parameters : Optional[Dict[str, float]], optional
        Parameter for the mesh, by default None
    fiber_parameters : Optional[Dict[str, Union[float, str]]], optional
        Parameter for the fibers, by default None

    Returns
    -------
    LVGeometry
        _description_
    """
    if not path.is_file():
        mesh_parameters = _update_parameters(
            LVGeometry.default_mesh_parameters(),
            mesh_parameters,
        )
        fiber_parameters = _update_parameters(
            LVGeometry.default_fiber_parameters(),
            fiber_parameters,
        )
        geo = LVGeometry.from_parameters(
            fiber_params=fiber_parameters,
            mesh_params=mesh_parameters,
        )
        geo.save(path)
    return LVGeometry.from_file(path)


def default_parameters():
    """Default parameters for Benchmark 1"""
    return dict(
        problem_parameters=LVProblem.default_parameters(),
        activation_parameters=activation_model.default_parameters(),
        pressure_parameters=pressure_model.default_parameters_benchmark1(),
        material_parameters=HolzapfelOgden.default_parameters(),
        mesh_parameters=LVGeometry.default_mesh_parameters(),
        fiber_parameters=LVGeometry.default_fiber_parameters(),
        outpath="results.h5",
        geometry_path="geometry.h5",
    )


def run(
    problem_parameters: Optional[Dict[str, Union[float, dolfin.Constant]]] = None,
    activation_parameters: Optional[Dict[str, float]] = None,
    pressure_parameters: Optional[Dict[str, float]] = None,
    material_parameters: Optional[Dict[str, Union[float, dolfin.Constant]]] = None,
    mesh_parameters: Optional[Dict[str, float]] = None,
    fiber_parameters: Optional[Dict[str, Union[float, str]]] = None,
    zero_pressure: bool = False,
    zero_activation: bool = False,
    outpath: Union[str, Path] = "results.h5",
    geometry_path: Union[str, Path] = "geometry.h5",
) -> None:
    """Run benchmark 1

    Parameters
    ----------
    problem_parameters : Optional[Dict[str, Union[float, dolfin.Constant]]], optional
        Parameters for the problem, by default None
    activation_parameters : Optional[Dict[str, float]], optional
        Parameters for the activation model, by default None
    pressure_parameters : Optional[Dict[str, float]], optional
        Parameters for the pressure model, by default None
    material_parameters : Optional[Dict[str, Union[float, dolfin.Constant]]], optional
        Parameters for the material model, by default None
    mesh_parameters : Optional[Dict[str, float]], optional
        Parameters for the mesh, by default None
    fiber_parameters : Optional[Dict[str, Union[float, str]]], optional
        Parameters for the fibers, by default None
    zero_pressure : bool, optional
        If True, set the pressure to zero, by default False
    zero_activation : bool, optional
        If True set the activation to zero, by default False
    outpath : Union[str, Path], optional
        Path to where to save the results, by default "results.h5"
    geometry_path : Union[str, Path], optional
        Path to the geometry, by default "geometry.h5"

    Raises
    ------
    OSError
        If output file is not an HDF5 file
    """
    outdir = Path(outpath).parent
    outdir.mkdir(parents=True, exist_ok=True)

    problem_parameters = _update_parameters(
        LVProblem.default_parameters(),
        problem_parameters,
    )
    pressure_parameters = _update_parameters(
        pressure_model.default_parameters_benchmark1(),
        pressure_parameters,
    )
    activation_parameters = _update_parameters(
        activation_model.default_parameters(),
        activation_parameters,
    )
    material_parameters = _update_parameters(
        HolzapfelOgden.default_parameters(),
        material_parameters,
    )
    mesh_parameters = _update_parameters(
        LVGeometry.default_mesh_parameters(),
        mesh_parameters,
    )
    fiber_parameters = _update_parameters(
        LVGeometry.default_fiber_parameters(),
        fiber_parameters,
    )

    dt = float(problem_parameters["dt"])
    tau = dolfin.Constant(0.0)
    time = np.arange(dt, 1, dt)

    t_eval = time - float(problem_parameters["alpha_f"]) * dt
    pressure = pressure_model.pressure_function(
        (0, 1),
        t_eval=t_eval,
        parameters=pressure_parameters,
    )
    if zero_pressure:
        # We set the pressure to zero
        pressure[:] = 0.0

    activation = activation_model.activation_function(
        (0, 1),
        t_eval=t_eval,
        parameters=activation_parameters,
    )
    if zero_activation:
        activation[:] = 0.0

    if dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
        np.save(
            outdir / "pressure_model.npy",
            {
                "time": time,
                "activation": activation,
                "pressure": pressure,
                "pressure_parameters": pressure_parameters,
                "activation_parameters": activation_parameters,
            },
        )

        postprocess.plot_activation_pressure_function(
            t=time,
            activation=activation,
            lv_pressure=pressure,
            outdir=outdir,
        )

    p = dolfin.Constant(0.0)
    problem_parameters["p"] = p

    geo = get_geometry(
        path=Path(geometry_path),
        fiber_parameters=fiber_parameters,
        mesh_parameters=mesh_parameters,
    )
    material = HolzapfelOgden(
        f0=geo.f0,
        s0=geo.s0,
        tau=tau,
        parameters=material_parameters,
    )
    problem = LVProblem(
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
        pressure_parameters={"lv": pressure_parameters},
        activation_parameters=activation_parameters,
    )

    solve(
        problem=problem,
        tau=tau,
        activation=activation,
        pressure=pressure,
        p=p,
        time=time,
        collector=collector,
        store_freq=1,
    )
