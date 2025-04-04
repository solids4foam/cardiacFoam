import weakref
from pathlib import Path
from typing import Dict
from typing import List
from typing import NamedTuple
from typing import Optional
from typing import Tuple
from typing import Union

import dolfin
import h5py
import matplotlib.pyplot as plt
import numpy as np

try:
    import ufl_legacy as ufl
except ImportError:
    import ufl
from dolfin import FiniteElement  # noqa: F401
from dolfin import tetrahedron  # noqa: F401
from dolfin import VectorElement  # noqa: F401

from .geometry import BiVGeometry
from .geometry import load_geometry
from .geometry import save_geometry
from .material import HolzapfelOgden
from .problem import BiVProblem  # noqa: F401
from .problem import LVProblem  # noqa: F401
from .problem import Problem


class SavedProblem(NamedTuple):
    problem: Problem
    time_stamps_str: np.ndarray
    time_stamps: np.ndarray
    signature: ufl.finiteelement.FiniteElementBase


def close_h5file(h5file):
    if h5file is not None:
        h5file.close()


def close_h5pyfile(h5pyfile):
    if h5pyfile is not None:
        h5pyfile.__exit__()


def save_problem(
    fname,
    problem: Problem,
    pressure_parameters: Optional[Dict[str, Dict[str, float]]],
    activation_parameters: Optional[Dict[str, float]],
):
    path = Path(fname)
    if path.is_file():
        path.unlink()
    save_geometry(fname, problem.geometry)

    comm = dolfin.MPI.comm_world
    if problem.geometry.mesh is not None:
        comm = problem.geometry.mesh.mpi_comm()

    if dolfin.MPI.rank(comm) == 0:
        with h5py.File(fname, "a") as h5file:
            group = h5file.create_group("problem_parameters")
            for k, v in problem.parameters.items():
                if isinstance(v, str):
                    continue  # Skip the function space
                group.create_dataset(k, data=float(v))
            group.attrs.create("cls", type(problem).__name__)

            group = h5file.create_group("material_parameters")
            for k, v in problem.material.parameters.items():
                group.create_dataset(k, data=float(v))

            if pressure_parameters is not None:
                group = h5file.create_group("pressure_parameters")

                for key in ["lv", "rv"]:
                    if key not in pressure_parameters:
                        continue
                    subgroup = group.create_group(key)
                    for k, v in pressure_parameters[key].items():
                        subgroup.create_dataset(k, data=float(v))

            if activation_parameters is not None:
                group = h5file.create_group("activation_parameters")
                for k, v in activation_parameters.items():
                    group.create_dataset(k, data=float(v))

            h5file.create_group("plv")
            if isinstance(problem, BiVProblem):
                h5file.create_group("prv")
            h5file.create_group("tau")


def load_problem(fname) -> SavedProblem:
    path = Path(fname)
    if not path.is_file():
        raise FileNotFoundError(f"File {path} does not exist")

    required_keys = ["mesh", "u", "problem_parameters", "material_parameters"]

    material_parameters = {}  # type: ignore
    problem_parameters = {}  # type: ignore
    signature = None  # type: ignore
    time_stamps_str = None  # type: ignore
    time_stamps = None  # type: ignore
    cls = None  # type: ignore
    if dolfin.MPI.rank(dolfin.MPI.comm_world) == 0:
        with h5py.File(path, "r") as h5file:
            # Check that we have the required keys
            for key in required_keys:
                if key not in h5file:
                    raise ValueError(f"Missing {key} in results file. Cannot load data")

            time_stamps_str = sorted(h5file["u"].keys(), key=lambda x: float(x))
            time_stamps = np.array(time_stamps_str, dtype=float)

            if len(time_stamps_str) == 0:
                raise ValueError(
                    "No timestamps found in results file. Cannot load data",
                )

            signature = h5file["u"][time_stamps_str[0]].attrs["signature"]
            signature = eval(signature)

            for k, v in h5file["material_parameters"].items():
                material_parameters[k] = v[...].tolist()

            for k, v in h5file["problem_parameters"].items():
                problem_parameters[k] = v[...].tolist()
            cls = eval(h5file["problem_parameters"].attrs["cls"])

    material_parameters = dolfin.MPI.comm_world.bcast(material_parameters, root=0)
    signature = dolfin.MPI.comm_world.bcast(signature, root=0)
    time_stamps = dolfin.MPI.comm_world.bcast(time_stamps, root=0)
    time_stamps_str = dolfin.MPI.comm_world.bcast(time_stamps_str, root=0)
    problem_parameters = dolfin.MPI.comm_world.bcast(problem_parameters, root=0)
    problem_parameters["function_space"] = f"{signature.family()}_{signature.degree()}"
    cls = dolfin.MPI.comm_world.bcast(cls, root=0)
    geometry = load_geometry(path)

    tau = dolfin.Constant(0.0)
    material = HolzapfelOgden(
        f0=geometry.f0,
        s0=geometry.s0,
        tau=tau,
        parameters=material_parameters,
    )

    problem = cls(
        geometry=geometry,
        material=material,
        parameters=problem_parameters,
    )

    return SavedProblem(
        problem=problem,
        time_stamps_str=time_stamps_str,
        time_stamps=time_stamps,
        signature=signature,
    )


class DataCollector:
    """Class for saving the simulation results

    Parameters
    ----------
    path : str
        HDF5 file to where to save the results
    problem : :class:`problem.Problem`
        The Mechanics problem
    pressure_parameters : Optional[Dict[str, Dict[str, float]]], optional
        Parameters used in the pressure model, by default None
    activation_parameters : Optional[Dict[str, float]], optional
        Parameters used in the activation model, by default None
    """

    def __init__(
        self,
        path,
        problem: Problem,
        pressure_parameters: Optional[Dict[str, Dict[str, float]]] = None,
        activation_parameters: Optional[Dict[str, float]] = None,
    ) -> None:
        self._path = Path(path)
        self._comm = dolfin.MPI.comm_world

        dolfin.MPI.barrier(self._comm)
        if self._comm.rank == 0:
            if self._path.is_file():
                # Delete file if is allready exist
                self._path.unlink()
        dolfin.MPI.barrier(self._comm)

        self._path = path
        self.u = problem.u_old
        self.v = problem.v_old
        self.a = problem.a_old
        if isinstance(problem, LVProblem):
            self.plv = problem.parameters["p"]
            self._problem_is_biv = False
        else:
            self.plv = problem.parameters["plv"]
            self.prv = problem.parameters["prv"]
            self._problem_is_biv = True

        self.tau = problem.material.tau

        if problem.geometry.mesh is not None:
            self._comm = problem.geometry.mesh.mpi_comm()

        save_problem(path, problem, pressure_parameters, activation_parameters)

    @property
    def path(self) -> str:
        """Path to the results as a string"""
        return self._path.as_posix()

    def store(self, t: float) -> None:
        """Save data at given time step

        Parameters
        ----------
        t : float
            Time stamp
        """
        with dolfin.HDF5File(self._comm, self.path, "a") as h5file:
            h5file.write(self.u, f"/u/{t:.4f}")
            h5file.write(self.v, f"/v/{t:.4f}")
            h5file.write(self.a, f"/a/{t:.4f}")

        if dolfin.MPI.rank(self._comm) == 0:
            with h5py.File(self.path, "a") as h5file:
                h5file["plv"].create_dataset(f"{t:.4f}", data=float(self.plv))
                if self._problem_is_biv:
                    h5file["prv"].create_dataset(f"{t:.4f}", data=float(self.prv))

                h5file["tau"].create_dataset(f"{t:.4f}", data=float(self.tau))


class DataLoader:
    def __init__(
        self,
        result_file: Union[str, Path],
    ) -> None:
        self._h5file_py = None
        self._h5file = None
        self._path = Path(result_file)
        if not self._path.is_file():
            raise FileNotFoundError(f"File {result_file} does not exist")

        self._problem = load_problem(self._path)

        if self.geometry.mesh is None:
            raise RuntimeError("Cannot load mesh")

        mesh = self.geometry.mesh
        self.comm = mesh.mpi_comm()
        x_max = self.geometry.mesh.coordinates().max(0)[0]
        self._shift = dolfin.Constant((x_max, 0, 0))

        self._h5file = dolfin.HDF5File(
            mesh.mpi_comm(),
            self.path,
            "r",
        )

        if dolfin.MPI.rank(self.comm) == 0:
            self._h5file_py = h5py.File(self.path, "r")

        if self.geometry.markers is None:
            raise RuntimeError("Cannot load markers")

        self.ds = ufl.Measure(
            "exterior_facet",
            domain=mesh,
            subdomain_data=self.geometry.ffun,
        )

        V = dolfin.FunctionSpace(mesh, self.signature)
        self.u = dolfin.Function(V, name="displacement")
        self.v = dolfin.Function(V, name="velocity")
        self.a = dolfin.Function(V, name="acceleration")

        self.stress_space = dolfin.FunctionSpace(self.geometry.mesh, "DG", 1)
        self._finalizer_h5file = weakref.finalize(self, close_h5file, self._h5file)
        self._finalizer_h5pyfile = weakref.finalize(
            self,
            close_h5pyfile,
            self._h5file_py,
        )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    def cleanup(self):
        self._finalizer_h5file()
        self._finalizer_h5pyfile()

    @property
    def problem(self):
        return self._problem.problem

    @property
    def signature(self):
        return self._problem.signature

    @property
    def time_stamps(self):
        return self._problem.time_stamps

    @property
    def time_stamps_str(self):
        return self._problem.time_stamps_str

    @property
    def geometry(self):
        return self.problem.geometry

    @property
    def geo_is_biv(self):
        return isinstance(self.geometry, BiVGeometry)

    def _check_t(self, t):
        if not isinstance(t, str):
            t = f"{t:.4f}"

        if t not in self.time_stamps_str:
            raise KeyError(f"Invalid time stamp {t}")

    def get_u(self, t):
        self._check_t(t)
        self._h5file.read(self.u, f"u/{t}/")
        return self.u

    def get_u_p_tau(self, t):
        self._check_t(t)
        self._h5file.read(self.u, f"u/{t}/")
        plv = []
        prv = []
        tau = []
        if dolfin.MPI.rank(self.comm) == 0:
            plv = self._h5file_py["plv"][t][...].tolist()
            if self.geo_is_biv:
                prv = self._h5file_py["prv"][t][...].tolist()

            tau = self._h5file_py["tau"][t][...].tolist()
        plv = self.comm.bcast(plv, root=0)
        prv = self.comm.bcast(prv, root=0)
        tau = self.comm.bcast(tau, root=0)
        return (self.u, plv, prv, tau)

    def get(self, t):
        self._check_t(t)

        self._h5file.read(self.u, f"u/{t}/")
        self._h5file.read(self.v, f"v/{t}/")
        self._h5file.read(self.a, f"a/{t}/")
        return self.u, self.v, self.a

    @property
    def path(self) -> str:
        return self._path.as_posix()

    def von_Mises_stress_at_point(
        self,
        point: Tuple[float, float, float],
    ) -> List[float]:
        print(f"Compute von Mises stress at point {point}")
        stress = []
        for t in self.time_stamps_str:
            print(f"Time {t}", end="\r")
            stress.append(self.von_Mises_stress_at_timepoint(t)(point))
        return stress

    def von_Mises_stress_at_timepoint(self, t):
        u, plv, prv, tau = self.get_u_p_tau(t)
        self.problem.u.assign(u)
        if self.geo_is_biv:
            self.problem.parameters["plv"] = plv
            self.problem.parameters["prv"] = prv
        else:
            self.problem.parameters["p"] = plv
        self.problem.parameters["tau"] = tau
        return dolfin.project(self.problem.von_Mises(), self.stress_space)

    def to_xdmf(self, path):
        print("Write displacement to xdmf")
        xdmf = dolfin.XDMFFile(self.geometry.mesh.mpi_comm(), path)
        for t in self.time_stamps:
            print(f"Time {t}", end="\r")
            u, v, a = self.get(t)

            xdmf.write_checkpoint(u, "displacement", t, append=True)
            xdmf.write_checkpoint(v, "velocity", t, append=True)
            xdmf.write_checkpoint(a, "acceleration", t, append=True)
        xdmf.close()

    def deformation_at_point(self, p):
        print(f"Compute deformation at point {p}")
        us = []
        for t in self.time_stamps:
            print(f"Deformation at time {t}", end="\r")
            us.append(self._deformation_at_time_point(self.get(t)(p)))

        return np.array(us)

    def _volume_form(self, u):
        X = ufl.SpatialCoordinate(self.geometry.mesh) - self._shift
        N = ufl.FacetNormal(self.geometry.mesh)
        F = ufl.grad(u) + ufl.Identity(3)

        u1 = ufl.as_vector([0.0, u[1], 0.0])
        X1 = ufl.as_vector([0.0, X[1], 0.0])
        return (-1.0 / 1.0) * ufl.dot(X1 + u1, ufl.cofac(F) * N)

    def _volume_at_timepoint(self, u, marker):
        return dolfin.assemble(self._volume_form(u) * self.ds(marker))

    def cavity_volume(self, marker):
        print("Compute cavity volume...")
        vols = []
        for t in self.time_stamps:
            print(f"Volume at time {t}", end="\r")
            u = self.get_u(t)
            vols.append(self._volume_at_timepoint(u, marker=marker))
        return np.array(vols)

    def postprocess_all(self, folder: Optional[Union[str, Path]] = None):
        if folder is None:
            folder = self._path.with_suffix("")
        outfolder = Path(folder)
        outfolder.mkdir(parents=True, exist_ok=True)

        np.save(outfolder / "time_stamps.npy", self.time_stamps)

        comm = self.geometry.mesh.mpi_comm()

        u_path = outfolder / "displacement.xdmf"
        u_xdmf = dolfin.XDMFFile(comm, u_path.as_posix())

        von_Mises_path = outfolder / "von_Mises_stress.xdmf"
        von_mises_xdmf = dolfin.XDMFFile(comm, von_Mises_path.as_posix())
        s = dolfin.Function(self.stress_space)
        p0 = (0.025, 0.03, 0)
        p1 = (0, 0.03, 0)
        # Only for benchmark 2
        p2 = (0.025, 0, 0.072)

        lv_vols = []
        rv_vols = []
        if self.geo_is_biv:
            lv_marker = self.geometry.markers["ENDO_LV"][0]
            rv_marker = self.geometry.markers["ENDO_RV"][0]
        else:
            lv_marker = self.geometry.markers["ENDO"][0]
            rv_marker = None

        up0 = []
        up1 = []
        up2: List[float] = []
        von_mises_p0 = []
        von_mises_p1 = []
        von_mises_p2: List[float] = []
        for t in self.time_stamps_str:
            print(f"Time {t}", end="\r")
            u = self.get_u(t)

            s.assign(self.von_Mises_stress_at_timepoint(t))
            von_mises_p0.append(s(p0))
            von_mises_p1.append(s(p1))

            u_xdmf.write(u, float(t))
            von_mises_xdmf.write(s, float(t))

            up0.append(u(p0))
            up1.append(u(p1))

            if self.geo_is_biv:
                up2.append(u(p2))
                von_mises_p2.append(s(p2))

            lv_vols.append(self._volume_at_timepoint(u, marker=lv_marker))
            if rv_marker is not None:
                rv_vols.append(self._volume_at_timepoint(u, marker=rv_marker))
        u_xdmf.close()
        von_mises_xdmf.close()

        plot_componentwise_displacement(
            up0=np.array(up0),
            up1=np.array(up1),
            up2=np.array(up2),
            time_stamps=self.time_stamps,
            fname=outfolder / "componentwise_displacement.png",
        )
        plot_von_Mises_stress(
            sp0=np.array(von_mises_p0),
            sp1=np.array(von_mises_p1),
            sp2=np.array(von_mises_p2),
            time_stamps=self.time_stamps,
            fname=outfolder / "von_Mises_stress.png",
        )
        plot_volume(
            lv_volumes=lv_vols,
            rv_volumes=rv_vols,
            time_stamps=self.time_stamps,
            fname=outfolder / "volume.png",
        )


class DisplacementData(NamedTuple):
    up0: np.ndarray
    up1: np.ndarray
    time: np.ndarray


def load_true_displacement_data(disp_path: Union[str, Path]) -> DisplacementData:
    disp_path = Path(disp_path)
    if not disp_path.is_file():
        raise FileNotFoundError(f"Cannot find file {disp_path}")
    results_disp = np.load(disp_path, allow_pickle=True)
    time = results_disp["times.npy"]
    results_disp_dict = results_disp["disp_dict.npy"].item()
    up0 = np.array(list(results_disp_dict["top"].values())).T
    up1 = np.array(list(results_disp_dict["middle"].values())).T

    return DisplacementData(up0=up0, up1=up1, time=time)


def load_true_volume_data(vol_path: Union[str, Path]) -> np.ndarray:
    vol_path = Path(vol_path)
    if not vol_path.is_file():
        raise FileNotFoundError(f"Cannot find file {vol_path}")
    results = np.load(vol_path, allow_pickle=True)
    return results["vol_lst.npy"]


def plot_componentwise_displacement(
    up0: np.ndarray,
    up1: np.ndarray,
    time_stamps: np.ndarray,
    up2: Optional[np.ndarray] = None,
    fname="componentwise_displacement.png",
):
    if up2 is None:
        up2 = np.array([])
    fname = Path(fname).with_suffix(".png")

    N = 3 if len(up2) > 0 else 2

    basefname = fname.with_suffix("").as_posix()
    np.save(Path(basefname + "_up0").with_suffix(".npy"), up0)
    np.save(Path(basefname + "_up1").with_suffix(".npy"), up1)
    if N == 3:
        np.save(Path(basefname + "_up2").with_suffix(".npy"), up2)

    fig, ax = plt.subplots(N, 1, sharex=True)
    ax[0].plot(time_stamps, up0[:, 0], label="x")
    ax[0].plot(time_stamps, up0[:, 1], label="y")
    ax[0].plot(time_stamps, up0[:, 2], label="z")
    ax[0].set_ylabel("$u(p_0)$[m]")

    ax[1].plot(time_stamps, up1[:, 0], label="x")
    ax[1].plot(time_stamps, up1[:, 1], label="y")
    ax[1].plot(time_stamps, up1[:, 2], label="z")
    ax[1].set_ylabel("$u(p_1)$[m]")
    ax[1].set_xlabel("Time [s]")

    if N == 3:
        ax[2].plot(time_stamps, up2[:, 0], label="x")
        ax[2].plot(time_stamps, up2[:, 1], label="y")
        ax[2].plot(time_stamps, up2[:, 2], label="z")
        ax[2].set_ylabel("$u(p_2)$[m]")
        ax[2].set_xlabel("Time [s]")

    for axi in ax:
        axi.legend()
        axi.grid()
    fig.savefig(fname, dpi=300)


def plot_von_Mises_stress(
    sp0: np.ndarray,
    sp1: np.ndarray,
    time_stamps: np.ndarray,
    sp2: Optional[np.ndarray] = None,
    fname="von_Mises_stress.png",
):
    if sp2 is None:
        sp2 = np.array([])

    fname = Path(fname).with_suffix(".png")

    basefname = fname.with_suffix("").as_posix()
    np.save(Path(basefname + "_sp0").with_suffix(".npy"), sp0)
    np.save(Path(basefname + "_sp1").with_suffix(".npy"), sp1)

    fig, ax = plt.subplots()
    ax.plot(time_stamps, sp0, label=r"$\sigma_v(p_0)$")
    ax.plot(time_stamps, sp1, label=r"$\sigma_v(p_1)$")
    if len(sp2) > 0:
        np.save(Path(basefname + "_sp2").with_suffix(".npy"), sp2)
        ax.plot(time_stamps, sp2, label=r"$\sigma_v(p_2)$")
    ax.set_ylabel("von Mises stress [Pa]")
    ax.set_xlabel("Time [s]")
    ax.legend()
    ax.grid()
    fig.savefig(fname, dpi=300)


def plot_volume(lv_volumes, time_stamps, rv_volumes=None, fname="volume.png"):
    if rv_volumes is None:
        rv_volumes = []
    fname = Path(fname).with_suffix(".png")
    basefname = fname.with_suffix("").as_posix()

    volume_data = {"lv": lv_volumes, "time": time_stamps}

    fig, ax = plt.subplots()
    ax.plot(time_stamps, lv_volumes, label="LV volume")
    if len(rv_volumes) > 0:
        volume_data["rv"] = rv_volumes
        ax.plot(time_stamps, rv_volumes, label="RV volume")
        ax.legend()
    ax.set_ylabel("Volume [m^3]")
    ax.set_xlabel("Time [s]")
    ax.grid()

    ax.set_title("Volume throug time")
    fig.savefig(fname, dpi=300)

    np.save(Path(basefname).with_suffix(".npy"), volume_data, allow_pickle=True)


def plot_volume_comparison(
    volumes,
    time_stamps,
    volumes_cmp=None,
    time_stamps_cmp=None,
    fname="volume_comparison.png",
):
    fname = Path(fname).with_suffix(".png")
    basefname = fname.with_suffix("").as_posix()
    np.save(Path(basefname).with_suffix(".npy"), volumes)

    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(time_stamps, volumes, label="volume")
    ax[0].plot(time_stamps_cmp, volumes_cmp, label="volume (true)")
    ax[1].plot(
        time_stamps,
        volumes - np.interp(time_stamps, time_stamps_cmp, volumes_cmp),
        label="volume - volume (true)",
    )
    ax[0].set_ylabel("Volume [$m^3$]")
    ax[1].set_xlabel("Time [s]")
    for axi in ax:
        axi.grid()
        axi.legend()
    ax[0].set_title("Volume through time")
    fig.tight_layout()
    fig.savefig(fname, dpi=300)


def plot_activation_pressure_function(
    outdir,
    t,
    activation,
    lv_pressure,
    rv_pressure=None,
):
    fig, ax = plt.subplots()
    ax.plot(t, activation)
    ax.set_title("Activation fuction \u03C4(t)")
    ax.set_ylabel("Pressure [Pa]")
    ax.set_xlabel("Time [s]")
    fig.savefig(Path(outdir) / "activation_function.png")

    fig, ax = plt.subplots()
    ax.plot(t, lv_pressure, label="LV")
    if rv_pressure is not None:
        ax.plot(t, rv_pressure, label="RV")
        ax.legend()
    ax.set_title("Pressure fuction p(t)")
    ax.set_ylabel("Pressure [Pa]")
    ax.set_xlabel("Time [s]")
    fig.savefig(Path(outdir) / "pressure_function.png")
