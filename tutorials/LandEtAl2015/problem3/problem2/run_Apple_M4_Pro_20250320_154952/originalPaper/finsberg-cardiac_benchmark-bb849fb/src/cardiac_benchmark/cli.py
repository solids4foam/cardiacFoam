"""Console script for cardiac-benchmark."""
import datetime
import json
import logging
import pprint
from enum import Enum
from pathlib import Path
from typing import Optional

import dolfin
import typer

from .geometry import BiVGeometry
from .geometry import LVGeometry
from .postprocess import DataLoader
from .utils import ConstantEncoder


app = typer.Typer()
logger = logging.getLogger(__name__)


class Resolution(str, Enum):
    fine = "fine"
    coarse = "coarse"


class Step0Case(str, Enum):
    a = "a"
    b = "b"


class Step2Case(str, Enum):
    a = "a"
    b = "b"
    c = "c"


def setup_logging(loglevel=logging.INFO):
    logging.basicConfig(
        level=loglevel,
        format=(
            "[%(asctime)s](proc %(process)d) - %(levelname)s - "
            "%(module)s:%(funcName)s:%(lineno)d %(message)s"
        ),
    )
    dolfin.set_log_level(logging.WARNING)
    for module in ["matplotlib", "h5py", "FFC", "UFL"]:
        logger = logging.getLogger(module)
        logger.setLevel(logging.WARNING)


def version_callback(show_version: bool):
    """Prints version information."""
    if show_version:
        from . import __version__, __program_name__

        typer.echo(f"{__program_name__} {__version__}")
        raise typer.Exit()


def license_callback(show_license: bool):
    """Prints license information."""
    if show_license:
        from . import __license__

        typer.echo(f"{__license__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Show version",
    ),
    license: bool = typer.Option(
        None,
        "--license",
        callback=license_callback,
        is_eager=True,
        help="Show license",
    ),
):
    # Do other global stuff, handle other global options here
    return


def benchmark1_step_0_1(
    step: int,
    case: str = "",
    geometry_path: Path = typer.Option(),
    outdir: Optional[Path] = typer.Option(None),
    run_benchmark: bool = True,
    run_postprocess: bool = True,
    alpha_m: float = 0.2,
    alpha_f: float = 0.4,
    zero_pressure: bool = False,
    zero_activation: bool = False,
    function_space: str = "P_2",
    loglevel: int = logging.INFO,
    a: float = 59.0,
    a_f: float = 18472.0,
    a_fs: float = 216.0,
    a_s: float = 2481.0,
    sigma_0: float = 150e3,
) -> int:
    setup_logging(loglevel=loglevel)
    if outdir is not None:
        outdir = Path(outdir).absolute()
    else:
        outdir = Path.cwd() / "results"

    if geometry_path is None:
        geometry_path = Path.cwd() / "geometry.h5"

    outdir.mkdir(exist_ok=True, parents=True)
    outpath = outdir / "result.h5"

    from . import benchmark1

    params = benchmark1.default_parameters()

    params["problem_parameters"]["alpha_m"] = alpha_m
    params["problem_parameters"]["alpha_f"] = alpha_f
    params["problem_parameters"]["function_space"] = function_space

    params["material_parameters"]["a"] = a
    params["material_parameters"]["a_f"] = a_f
    params["material_parameters"]["a_fs"] = a_fs
    params["material_parameters"]["a_s"] = a_s

    params["activation_parameters"]["sigma_0"] = sigma_0

    params["zero_pressure"] = zero_pressure
    params["zero_activation"] = zero_activation
    params["outpath"] = outpath.as_posix()
    params["geometry_path"] = geometry_path.as_posix()

    parameters = params.copy()
    parameters["benchmark"] = 1
    parameters["step"] = step
    parameters["case"] = case
    parameters["outdir"] = outdir.as_posix()
    parameters["timestamp"] = datetime.datetime.now().isoformat()

    logger.info(
        f"Running benchmark 1 step {step}, "
        f"case {case} with parameters {pprint.pformat(parameters)}",
    )
    if dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
        logger.info(f"Output will be saved to {outdir}")
        (outdir / "parameters.json").write_text(
            json.dumps(parameters, cls=ConstantEncoder),
        )

    if run_benchmark:
        benchmark1.run(**params)  # type: ignore

    loader = DataLoader(outpath)
    if run_postprocess:
        loader.postprocess_all(folder=outdir)

    return 0


@app.command(help="Run benchmark 1 - step 0")
def benchmark1_step0(
    case: Step0Case,
    geometry_path: Path = typer.Option(),
    outdir: Optional[Path] = typer.Option(None),
    run_benchmark: bool = True,
    run_postprocess: bool = True,
    alpha_m: float = 0.2,
    alpha_f: float = 0.4,
    function_space: str = "P_2",
    loglevel: int = logging.INFO,
) -> int:
    if outdir is None:
        outdir = Path("results_benchmark1_step0_caseA")

    if Step0Case[case] == Step0Case.a:
        zero_pressure = True
        zero_activation = False
    elif Step0Case[case] == Step0Case.b:
        zero_activation = True
        zero_pressure = False
    else:
        raise ValueError(f"Invalid case {Step0Case[case].name} for step 0")

    return benchmark1_step_0_1(
        step=0,
        case=Step0Case[case].name,
        outdir=outdir,
        run_benchmark=run_benchmark,
        run_postprocess=run_postprocess,
        alpha_m=alpha_m,
        alpha_f=alpha_f,
        zero_pressure=zero_pressure,
        zero_activation=zero_activation,
        geometry_path=geometry_path,
        function_space=function_space,
        loglevel=loglevel,
    )


@app.command(help="Run benchmark 1 - step 1")
def benchmark1_step1(
    geometry_path: Path = typer.Option(),
    outdir: Optional[Path] = typer.Option(None),
    run_benchmark: bool = True,
    run_postprocess: bool = True,
    alpha_m: float = 0.2,
    alpha_f: float = 0.4,
    function_space: str = "P_2",
    loglevel: int = logging.INFO,
) -> int:
    if outdir is None:
        outdir = Path("results_benchmark1_step1")
    benchmark1_step_0_1(
        step=1,
        outdir=outdir,
        run_benchmark=run_benchmark,
        run_postprocess=run_postprocess,
        alpha_m=alpha_m,
        alpha_f=alpha_f,
        zero_activation=False,
        zero_pressure=False,
        geometry_path=geometry_path,
        function_space=function_space,
        loglevel=loglevel,
    )

    return 0


@app.command(help="Run benchmark 1 - step 2")
def benchmark1_step2(
    case: Step2Case,
    geometry_path: Path = typer.Option(),
    outdir: Optional[Path] = typer.Option(None),
    run_benchmark: bool = True,
    run_postprocess: bool = True,
    alpha_m: float = 0.2,
    alpha_f: float = 0.4,
    function_space: str = "P_2",
    loglevel: int = logging.INFO,
) -> int:
    if outdir is None:
        outdir = Path(f"results_benchmark1_step2/case{Step2Case[case].name}")

    if Step2Case[case] == Step2Case.a:
        a = 177.0
        a_f = 55416.0
        a_fs = 648.0
        a_s = 7443.0
        sigma_0 = 2e5
    elif Step2Case[case] == Step2Case.b:
        a = 295.0
        a_f = 92360.0
        a_fs = 1080.0
        a_s = 12405.0
        sigma_0 = 1e5
    elif Step2Case[case] == Step2Case.c:
        a = 19.0
        a_f = 6157.0
        a_fs = 72.0
        a_s = 827.0
        sigma_0 = 2e5
    else:
        raise ValueError(f"Invalid case {Step2Case[case].name} for step 2")

    return benchmark1_step_0_1(
        step=2,
        case=case,
        outdir=outdir,
        run_benchmark=run_benchmark,
        run_postprocess=run_postprocess,
        alpha_m=alpha_m,
        alpha_f=alpha_f,
        geometry_path=geometry_path,
        function_space=function_space,
        loglevel=loglevel,
        a=a,
        a_f=a_f,
        a_fs=a_fs,
        a_s=a_s,
        sigma_0=sigma_0,
    )


@app.command(help="Run benchmark 2")
def benchmark2(
    geometry_path: Path = typer.Option(),
    outdir: Optional[Path] = typer.Option(None),
    run_benchmark: bool = True,
    run_postprocess: bool = True,
    alpha_m: float = 0.2,
    alpha_f: float = 0.4,
    function_space: str = "P_2",
    loglevel: int = logging.INFO,
    T: float = 1.0,
) -> int:
    setup_logging(loglevel=loglevel)
    if outdir is not None:
        outdir = Path(outdir).absolute()
    else:
        outdir = Path.cwd() / "results_benchmark2"

    outdir.mkdir(exist_ok=True, parents=True)
    outpath = outdir / "result.h5"

    from . import benchmark2

    params = benchmark2.default_parameters()
    params["outpath"] = outpath
    params["problem_parameters"]["alpha_m"] = alpha_m
    params["problem_parameters"]["alpha_f"] = alpha_f
    params["geometry_path"] = geometry_path.as_posix()
    params["T"] = T
    params["problem_parameters"]["function_space"] = function_space

    parameters = params.copy()
    parameters["benchmark"] = 2
    parameters["outdir"] = outdir.as_posix()
    parameters["timestamp"] = datetime.datetime.now().isoformat()

    logger.info(f"Running benchmakr 2 with parameters {pprint.pformat(parameters)}")
    if dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
        logger.info(f"Output will be saved to {outdir}")
        (outdir / "parameters.json").write_text(
            json.dumps(parameters, cls=ConstantEncoder),
        )

    if run_benchmark:
        benchmark2.run(**params)

    loader = DataLoader(outpath)
    if run_postprocess:
        loader.postprocess_all(folder=outdir)

    return 0


class FiberSpaces(str, Enum):
    P1 = "P1"
    P2 = "P2"


@app.command(help="Download and extract data for benchmark 1")
def download_data_benchmark1(
    outdir: Optional[Path] = typer.Option(None),
    fiber_space: FiberSpaces = FiberSpaces.P2,
):
    setup_logging(loglevel=logging.INFO)
    import urllib.request

    mesh_data = {
        "mesh_xdmf": "https://drive.google.com/uc?id=1xpEQB0EgKZCgdC2pG2C7cUIXkH1LvShe&export=download",
        "mesh_h5": "https://drive.google.com/uc?id=1xjLRJbOBViq8Lq3Ktk9STkjzcFtF6rTP&export=download",
    }
    fiber_data = {
        FiberSpaces.P1: {
            "fiber": "https://drive.google.com/uc?id=15_aXtjjeYC9T8uNbByhSBvHTbkfXkUWc&export=download",
            "sheet": "https://drive.google.com/uc?id=1tr2niWHHvkLMrL7hs0H_VorUb6WYhqSh&export=download",
            "sheet_normal": "https://drive.google.com/uc?id=1fF7vBsfAbEDd_bfElFc4QLlndRkl53rN&export=download",
        },
        FiberSpaces.P2: {
            "fiber": "https://drive.google.com/uc?id=1TRUSkndBDQIajA5-Ai3P_Po5rewTsepD&export=download",
            "sheet": "https://drive.google.com/uc?id=1mejFt0YGHSGAth-NvW7OuOZRMzbPFUY_&export=download",
            "sheet_normal": "https://drive.google.com/uc?id=1NdMJPrNbf3K7qm6iOF5sMSeoioUe5GOe&export=download",
        },
    }

    if outdir is None:
        outdir = Path.cwd() / "data_benchmark1"
    outdir.mkdir(exist_ok=True, parents=True)
    fiberdir = outdir / "fibers" / FiberSpaces[fiber_space].name.lower()
    fiberdir.mkdir(exist_ok=True, parents=True)
    logger.info(f"Download data for benchmark 1 to {outdir}")

    logger.info("Download ellipsoid_0.005.h5")
    urllib.request.urlretrieve(
        mesh_data["mesh_h5"],
        outdir / "ellipsoid_0.005.h5",
    )
    logger.info("Download ellipsoid_0.005.h5")
    urllib.request.urlretrieve(
        mesh_data["mesh_xdmf"],
        outdir / "ellipsoid_0.005.xdmf",
    )
    logger.info("Download ellipsoid_fiber.h5")
    urllib.request.urlretrieve(
        fiber_data[FiberSpaces[fiber_space]]["fiber"],
        fiberdir / "ellipsoid_fiber.h5",
    )
    logger.info("Download ellipsoid_sheet.h5")
    urllib.request.urlretrieve(
        fiber_data[FiberSpaces[fiber_space]]["sheet"],
        fiberdir / "ellipsoid_sheet.h5",
    )
    logger.info("Download ellipsoid_sheet_normal.h5")
    urllib.request.urlretrieve(
        fiber_data[FiberSpaces[fiber_space]]["sheet_normal"],
        fiberdir / "ellipsoid_sheet_normal.h5",
    )


@app.command(help="Download and extract data for benchmark 2")
def download_data_benchmark2(
    resolution: Resolution,
    outdir: Optional[Path] = typer.Option(None),
):
    setup_logging(loglevel=logging.INFO)
    import urllib.request

    data = {
        Resolution.coarse: {
            "mesh_xdmf": "https://drive.google.com/uc?id=140Qe5hNkYuCHI6Zoe8NkL7xvB3gleAt5&export=download",
            "mesh_h5": "https://drive.google.com/uc?id=16FsaTNfexETHxHS_y4U0rFACOnVI1j_j&export=download",
            "fiber": "https://drive.google.com/uc?id=1m8fajRggnVsimKdX_FlgqcXbiHIwLB-M&export=download",
            "sheet": "https://drive.google.com/uc?id=1NGCtHPCnAd-KWL9fuA8JC4EiWC1JBWd8&export=download",
            "sheet_normal": "https://drive.google.com/uc?id=15XfHDPkG4kCZMH-LC2mBORASD_LTKDb_&export=download",
            # "mesh_xdmf": "https://drive.google.com/uc?id=1ijbYNyEFxmRAgdztN2Aw9Ll_VsWz7BSE&export=download",
            # "mesh_h5": "https://drive.google.com/uc?id=1ZOl2NWeQbQVhKcLFlC81ZTEZ-EY_Dwrj&export=download",
            # "fiber": "https://drive.google.com/uc?id=1ZIlNPxLxdCVdMqTlZzoWdliBqKN-qNWd&export=download",
            # "sheet": "https://drive.google.com/uc?id=1EmYFua-_dFU3_OY_0PJVmnOS2QFyYIf4&export=download",
            # "sheet_normal": "https://drive.google.com/uc?id=1rejGIziVkAJFssLFF38oVMaWh9a4F94I&export=download",
        },
        Resolution.fine: {
            "mesh_xdmf": "https://drive.google.com/uc?id=1W5Ui7sWU5tj_dEMMlKBc_ftmGA5Zc_HJ&export=download",
            "mesh_h5": "https://drive.google.com/uc?id=1dyPVle3-PlSMQx9IcDWT1EprsoM2ACjV&export=download",
            "fiber": "https://drive.google.com/uc?id=18zV-5nUV9ZMFt9sZSn3BFmbuNgvy6TV6&export=download",
            "sheet": "https://drive.google.com/uc?id=1WF1AqlsQjn_ypVfWjZOhdt_o-HnDbxUu&export=download",
            "sheet_normal": "https://drive.google.com/uc?id=1Zxr_NaM9Ee1rpHVJWYMnm0AJjNHDx_0L&export=download",
            # "mesh_xdmf": "https://drive.google.com/uc?id=1CTG17jPe6y6aZMnqmSM4n-OHgc9S3QWN&export=download",
            # "mesh_h5": "https://drive.google.com/uc?id=16f9LXN3jB-nEDtvhbIA1O0QZB6cCGd78&export=download",
            # "fiber": "https://drive.google.com/uc?id=1WzTMCfwy1Sxjnp7cmPT5WW6vBQ15aET-&export=download",
            # "sheet": "https://drive.google.com/uc?id=13ESPrIFBvy2Qi8FflcZtefmy2obvpuiz&export=download",
            # "sheet_normal": "https://drive.google.com/uc?id=1gTORU_BarcNKJAiFc9WjSLwGAIomJghN&export=download",
        },
    }

    if outdir is None:
        outdir = Path.cwd() / f"data_benchmark2_{resolution.name}"
    outdir.mkdir(exist_ok=True, parents=True)
    fiberdir = outdir / "fibers"
    fiberdir.mkdir(exist_ok=True)
    logger.info(f"Download {resolution.name} data for benchmark to to {outdir}")

    logger.info("Download bi_ventricular.h5")
    urllib.request.urlretrieve(
        data[resolution]["mesh_h5"],
        outdir / "bi_ventricular.h5",
    )
    logger.info("Download bi_ventricular.xdmf")
    urllib.request.urlretrieve(
        data[resolution]["mesh_xdmf"],
        outdir / "bi_ventricular.xdmf",
    )
    logger.info("Download bi_ventricular_fiber.h5")
    urllib.request.urlretrieve(
        data[resolution]["fiber"],
        fiberdir / "bi_ventricular_fiber.h5",
    )
    logger.info("Download bi_ventricular_sheet.h5")
    urllib.request.urlretrieve(
        data[resolution]["sheet"],
        fiberdir / "bi_ventricular_sheet.h5",
    )
    logger.info("Download bi_ventricular_sheet_normal.h5")
    urllib.request.urlretrieve(
        data[resolution]["sheet_normal"],
        fiberdir / "bi_ventricular_sheet_normal.h5",
    )


@app.command(help="Convert downloaded data for benchmark 1")
def convert_data_benchmark1(
    data_folder: Path,
    outpath: Optional[Path] = typer.Option(None),
    outdir: Optional[Path] = typer.Option(None),
    to_paraview: bool = True,
    fiber_space: FiberSpaces = FiberSpaces.P2,
):
    setup_logging()
    fiberdir = data_folder / "fibers" / FiberSpaces[fiber_space].name.lower()
    mesh_file = data_folder / "ellipsoid_0.005.xdmf"
    assert mesh_file.is_file(), f"Missing {mesh_file}"
    fiber_file = fiberdir / "ellipsoid_fiber.h5"
    assert fiber_file.is_file(), f"Missing {fiber_file}"
    sheet_file = fiberdir / "ellipsoid_sheet.h5"
    assert sheet_file.is_file(), f"Missing {sheet_file}"
    sheet_normal_file = fiberdir / "ellipsoid_sheet_normal.h5"
    assert sheet_normal_file.is_file(), f"Missing {sheet_normal_file}"

    if outdir is None:
        outdir = data_folder

    outdir.mkdir(exist_ok=True, parents=True)
    logger.info(f"Load geometry from {data_folder}")
    geo = LVGeometry.from_files(
        mesh_file=mesh_file,
        fiber_file=fiber_file,
        sheet_file=sheet_file,
        sheet_normal_file=sheet_normal_file,
    )

    if outpath is not None:
        geo.save(outpath)

    if not to_paraview:
        return

    mesh_path = (outdir / "mesh.xdmf").as_posix()
    logger.info(f"Save mesh to {mesh_path}")
    with dolfin.XDMFFile(mesh_path) as xdmf:
        xdmf.write(geo.mesh)

    ffun_path = (outdir / "ffun.xdmf").as_posix()
    logger.info(f"Save ffun to {ffun_path}")
    with dolfin.XDMFFile(ffun_path) as xdmf:
        xdmf.write(geo.ffun)

    microstructure_path = (outdir / "microstructure.xdmf").as_posix()
    logger.info(f"Save microstructure to {microstructure_path}")
    with dolfin.XDMFFile(microstructure_path) as xdmf:
        xdmf.write_checkpoint(geo.f0, "fiber", 0.0, dolfin.XDMFFile.Encoding.HDF5, True)
        xdmf.write_checkpoint(geo.s0, "sheet", 0.0, dolfin.XDMFFile.Encoding.HDF5, True)
        xdmf.write_checkpoint(
            geo.n0,
            "sheet_normal",
            0.0,
            dolfin.XDMFFile.Encoding.HDF5,
            True,
        )


@app.command(help="Convert downloaded data for benchmark 2")
def convert_data_benchmark2(
    data_folder: Path,
    outpath: Optional[Path] = typer.Option(None),
    outdir: Optional[Path] = typer.Option(None),
    to_paraview: bool = True,
):
    setup_logging()
    mesh_file = data_folder / "bi_ventricular.xdmf"
    assert mesh_file.is_file(), f"Missing {mesh_file}"
    fiber_file = data_folder / "fibers/bi_ventricular_fiber.h5"
    assert fiber_file.is_file(), f"Missing {fiber_file}"
    sheet_file = data_folder / "fibers/bi_ventricular_sheet.h5"
    assert sheet_file.is_file(), f"Missing {sheet_file}"
    sheet_normal_file = data_folder / "fibers/bi_ventricular_sheet_normal.h5"
    assert sheet_normal_file.is_file(), f"Missing {sheet_normal_file}"

    if outdir is None:
        outdir = data_folder

    outdir.mkdir(exist_ok=True, parents=True)
    logger.info(f"Load geometry from {data_folder}")
    geo = BiVGeometry.from_files(
        mesh_file=mesh_file,
        fiber_file=fiber_file,
        sheet_file=sheet_file,
        sheet_normal_file=sheet_normal_file,
    )

    if outpath is not None:
        geo.save(outpath)

    if not to_paraview:
        return

    mesh_path = (outdir / "mesh.xdmf").as_posix()
    logger.info(f"Save mesh to {mesh_path}")
    with dolfin.XDMFFile(mesh_path) as xdmf:
        xdmf.write(geo.mesh)

    ffun_path = (outdir / "ffun.xdmf").as_posix()
    logger.info(f"Save ffun to {ffun_path}")
    with dolfin.XDMFFile(ffun_path) as xdmf:
        xdmf.write(geo.ffun)

    microstructure_path = (outdir / "microstructure.xdmf").as_posix()
    logger.info(f"Save microstructure to {microstructure_path}")
    with dolfin.XDMFFile(microstructure_path) as xdmf:
        xdmf.write_checkpoint(geo.f0, "fiber", 0.0, dolfin.XDMFFile.Encoding.HDF5, True)
        xdmf.write_checkpoint(geo.s0, "sheet", 0.0, dolfin.XDMFFile.Encoding.HDF5, True)
        xdmf.write_checkpoint(
            geo.n0,
            "sheet_normal",
            0.0,
            dolfin.XDMFFile.Encoding.HDF5,
            True,
        )


@app.command(help="Create geometry")
def create_geometry(
    path: Path,
    alpha_endo: float = -60.0,
    alpha_epi: float = 60.0,
    function_space: str = "P_2",
    mesh_size_factor: float = 1.0,
):
    geo = LVGeometry.from_parameters(
        fiber_params={
            "alpha_endo": alpha_endo,
            "alpha_epi": alpha_epi,
            "function_space": function_space,
        },
        mesh_params={
            "mesh_size_factor": mesh_size_factor,
        },
    )
    dolfin.File("mesh.pvd") << geo.mesh
    geo.save(path)
