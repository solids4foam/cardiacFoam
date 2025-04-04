import abc
import logging
import math
import warnings
from collections import namedtuple
from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Tuple
from typing import Union

import dolfin

try:
    import gmsh
except ImportError:
    warnings.warn("gmsh not installed - mesh generation not possible")
    hash_gmsh = False
else:
    hash_gmsh = True
import h5py
import meshio
from dolfin import FiniteElement  # noqa: F401
from dolfin import tetrahedron  # noqa: F401
from dolfin import VectorElement  # noqa: F401

from .microstructure import create_microstructure

GmshGeometry = namedtuple(
    "GmshGeometry",
    ["mesh", "cfun", "ffun", "efun", "vfun", "markers"],
)


logger = logging.getLogger(__name__)


def create_mesh(mesh, cell_type):
    # From http://jsdokken.com/converted_files/tutorial_pygmsh.html
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out_mesh = meshio.Mesh(
        points=mesh.points,
        cells={cell_type: cells},
        cell_data={"name_to_read": [cell_data]},
    )
    return out_mesh


def read_meshfunction(fname, obj):
    with dolfin.XDMFFile(Path(fname).as_posix()) as f:
        f.read(obj, "name_to_read")


def gmsh2dolfin(msh_file):
    logger.debug(f"Convert file {msh_file} to dolfin")
    msh = meshio.gmsh.read(msh_file)

    vertex_mesh = create_mesh(msh, "vertex")
    line_mesh = create_mesh(msh, "line")
    triangle_mesh = create_mesh(msh, "triangle")
    tetra_mesh = create_mesh(msh, "tetra")

    vertex_mesh_name = Path("vertex_mesh.xdmf")
    meshio.write(vertex_mesh_name, vertex_mesh)

    line_mesh_name = Path("line_mesh.xdmf")
    meshio.write(line_mesh_name, line_mesh)

    triangle_mesh_name = Path("triangle_mesh.xdmf")
    meshio.write(triangle_mesh_name, triangle_mesh)

    tetra_mesh_name = Path("mesh.xdmf")
    meshio.write(
        tetra_mesh_name,
        tetra_mesh,
    )

    mesh = dolfin.Mesh()

    with dolfin.XDMFFile(tetra_mesh_name.as_posix()) as infile:
        infile.read(mesh)

    cfun = dolfin.MeshFunction("size_t", mesh, 3)
    read_meshfunction(tetra_mesh_name, cfun)
    tetra_mesh_name.unlink()
    tetra_mesh_name.with_suffix(".h5").unlink()

    ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
    read_meshfunction(triangle_mesh_name, ffun_val)
    ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)
    for value in ffun_val.values():
        mesh.domains().set_marker(value, 2)
    ffun.array()[ffun.array() == max(ffun.array())] = 0
    triangle_mesh_name.unlink()
    triangle_mesh_name.with_suffix(".h5").unlink()

    efun_val = dolfin.MeshValueCollection("size_t", mesh, 1)
    read_meshfunction(line_mesh_name, efun_val)
    efun = dolfin.MeshFunction("size_t", mesh, efun_val)
    efun.array()[efun.array() == max(efun.array())] = 0
    line_mesh_name.unlink()
    line_mesh_name.with_suffix(".h5").unlink()

    vfun_val = dolfin.MeshValueCollection("size_t", mesh, 0)
    read_meshfunction(vertex_mesh_name, vfun_val)
    vfun = dolfin.MeshFunction("size_t", mesh, vfun_val)
    vfun.array()[vfun.array() == max(vfun.array())] = 0
    vertex_mesh_name.unlink()
    vertex_mesh_name.with_suffix(".h5").unlink()

    markers = msh.field_data

    return GmshGeometry(
        mesh=mesh,
        vfun=vfun,
        efun=efun,
        ffun=ffun,
        cfun=cfun,
        markers=markers,
    )


def create_benchmark_ellipsoid_mesh_gmsh(
    mesh_name,
    r_short_endo=0.025,
    r_short_epi=0.035,
    r_long_endo=0.09,
    r_long_epi=0.097,
    psize_ref=0.005,
    mu_apex_endo=-math.pi,
    mu_base_endo=-math.acos(5 / 17),
    mu_apex_epi=-math.pi,
    mu_base_epi=-math.acos(5 / 20),
    mesh_size_factor=1.0,
):
    if not hash_gmsh:
        raise ImportError("gmsh is not installed - unable to create mesh")
    gmsh.initialize()

    gmsh.option.setNumber("Geometry.CopyMeshingMethod", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    # gmsh.option.setNumber("Mesh.Algorithm3D", 7)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    # gmsh.option.setNumber("Mesh.Smoothing", 100)
    gmsh.option.setNumber("Mesh.MeshSizeFactor", mesh_size_factor)

    def ellipsoid_point(mu, theta, r_long, r_short, psize):
        return gmsh.model.geo.addPoint(
            r_long * math.cos(mu),
            r_short * math.sin(mu) * math.cos(theta),
            r_short * math.sin(mu) * math.sin(theta),
            psize,
        )

    center = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)

    apex_endo = ellipsoid_point(
        mu=mu_apex_endo,
        theta=0.0,
        r_short=r_short_endo,
        r_long=r_long_endo,
        psize=psize_ref / 2.0,
    )

    base_endo = ellipsoid_point(
        mu=mu_base_endo,
        theta=0.0,
        r_short=r_short_endo,
        r_long=r_long_endo,
        psize=psize_ref,
    )

    apex_epi = ellipsoid_point(
        mu=mu_apex_epi,
        theta=0.0,
        r_short=r_short_epi,
        r_long=r_long_epi,
        psize=psize_ref / 2.0,
    )

    base_epi = ellipsoid_point(
        mu=mu_base_epi,
        theta=0.0,
        r_short=r_short_epi,
        r_long=r_long_epi,
        psize=psize_ref,
    )

    apex = gmsh.model.geo.addLine(apex_endo, apex_epi)
    base = gmsh.model.geo.addLine(base_endo, base_epi)
    endo = gmsh.model.geo.add_ellipse_arc(apex_endo, center, apex_endo, base_endo)
    epi = gmsh.model.geo.add_ellipse_arc(apex_epi, center, apex_epi, base_epi)

    ll1 = gmsh.model.geo.addCurveLoop([apex, epi, -base, -endo])

    s1 = gmsh.model.geo.addPlaneSurface([ll1])

    sendoringlist = []
    sepiringlist = []
    sendolist = []
    sepilist = []
    sbaselist = []
    vlist = []

    out = [(2, s1)]
    for _ in range(4):
        out = gmsh.model.geo.revolve(
            [out[0]],
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            math.pi / 2,
        )

        sendolist.append(out[4][1])
        sepilist.append(out[2][1])
        sbaselist.append(out[3][1])
        vlist.append(out[1][1])

        gmsh.model.geo.synchronize()
        bnd = gmsh.model.getBoundary([out[0]])

        sendoringlist.append(bnd[1][1])
        sepiringlist.append(bnd[3][1])

    phys_apex_endo = gmsh.model.addPhysicalGroup(0, [apex_endo])
    gmsh.model.setPhysicalName(0, phys_apex_endo, "ENDOPT")

    phys_apex_epi = gmsh.model.addPhysicalGroup(0, [apex_epi])
    gmsh.model.setPhysicalName(0, phys_apex_epi, "EPIPT")

    phys_epiring = gmsh.model.addPhysicalGroup(1, sepiringlist)
    gmsh.model.setPhysicalName(1, phys_epiring, "EPIRING")

    phys_endoring = gmsh.model.addPhysicalGroup(1, sendoringlist)
    gmsh.model.setPhysicalName(1, phys_endoring, "ENDORING")

    phys_base = gmsh.model.addPhysicalGroup(2, sbaselist)
    gmsh.model.setPhysicalName(2, phys_base, "BASE")

    phys_endo = gmsh.model.addPhysicalGroup(2, sendolist)
    gmsh.model.setPhysicalName(2, phys_endo, "ENDO")

    phys_epi = gmsh.model.addPhysicalGroup(2, sepilist)
    gmsh.model.setPhysicalName(2, phys_epi, "EPI")

    phys_myo = gmsh.model.addPhysicalGroup(3, vlist)
    gmsh.model.setPhysicalName(3, phys_myo, "MYOCARDIUM")

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write(Path(mesh_name).as_posix())

    gmsh.finalize()


def save_geometry(fname, geo: "HeartGeometry") -> None:
    logger.info(f"Save geometry to {fname}")
    fname = Path(fname)
    if fname.is_file():
        logger.info(f"File {fname} already exists. Remove existing file")
        fname.unlink()

    mpi_comm = dolfin.MPI.comm_world
    if geo.mesh is not None:
        mpi_comm = geo.mesh.mpi_comm()
    with dolfin.HDF5File(mpi_comm, fname.as_posix(), "w") as h5file:
        for attr in ["mesh", "cfun", "ffun", "efun", "vfun", "f0", "s0", "n0"]:
            obj = getattr(geo, attr, None)
            if obj is None:
                continue

            h5file.write(obj, attr)
        if geo.markers is not None:
            for name, (marker, dim) in geo.markers.items():
                h5file.attributes("mesh")[name] = f"{marker}_{dim}"


def load_geometry(fname) -> "HeartGeometry":
    fname = Path(fname)
    if not fname.is_file():
        raise FileNotFoundError(f"File {fname} does not exist")
    logger.info(f"Load geometry from {fname}")
    mesh = dolfin.Mesh()

    element = None  # type: ignore
    markers: Dict[str, Tuple[int, int]] = {}

    if mesh.mpi_comm().rank == 0:
        with h5py.File(fname, "r") as h5file:
            element = eval(h5file["f0"].attrs["signature"])
            for k, v in h5file["mesh"].attrs.items():
                v_split = v.decode().split("_")
                markers[k] = (int(v_split[0]), int(v_split[1]))

    element = mesh.mpi_comm().bcast(element, root=0)
    markers = mesh.mpi_comm().bcast(markers, root=0)

    with dolfin.HDF5File(mesh.mpi_comm(), fname.as_posix(), "r") as h5file:
        h5file.read(mesh, "mesh", False)
        ffun = dolfin.MeshFunction("size_t", mesh, 2)
        h5file.read(ffun, "ffun")

        element._quad_scheme = "default"
        V = dolfin.FunctionSpace(mesh, element)
        f0 = dolfin.Function(V)
        h5file.read(f0, "f0")
        s0 = dolfin.Function(V)
        h5file.read(s0, "s0")
        n0 = dolfin.Function(V)
        h5file.read(n0, "n0")

    cls = BiVGeometry if "ENDO_RV" in markers else LVGeometry
    return cls(
        mesh=mesh,
        ffun=ffun,
        markers=markers,
        f0=f0,
        s0=s0,
        n0=n0,
    )


class HeartGeometry(abc.ABC):
    def __init__(
        self,
        mesh: Optional[dolfin.Mesh] = None,
        cfun: Optional[dolfin.MeshFunction] = None,
        ffun: Optional[dolfin.MeshFunction] = None,
        efun: Optional[dolfin.MeshFunction] = None,
        vfun: Optional[dolfin.MeshFunction] = None,
        markers: Optional[Dict[str, Tuple[int, int]]] = None,
        f0: Optional[dolfin.Function] = None,
        s0: Optional[dolfin.Function] = None,
        n0: Optional[dolfin.Function] = None,
    ):
        self.mesh = mesh
        self.cfun = cfun
        self.ffun = ffun
        self.efun = efun
        self.vfun = vfun
        self.f0 = f0
        self.s0 = s0
        self.n0 = n0
        self.markers = markers

    @classmethod
    def from_file(cls, fname):
        return load_geometry(fname)

    @staticmethod
    @abc.abstractmethod
    def load_mesh(mesh_file):
        ...

    @classmethod
    def from_files(
        cls,
        mesh_file: Union[str, Path],
        fiber_file: Union[str, Path, None] = None,
        sheet_file: Union[str, Path, None] = None,
        sheet_normal_file: Union[str, Path, None] = None,
    ) -> "HeartGeometry":
        mesh_file = Path(mesh_file)
        if not mesh_file.is_file():
            raise FileNotFoundError(f"File {mesh_file} does not exist")
        mesh, ffun, markers = cls.load_mesh(mesh_file)

        complete_microstructure = True
        if fiber_file is None:
            complete_microstructure = False
            logger.warning("Missing file for fibers. ")
        else:
            fiber_file = Path(fiber_file)

        if sheet_file is None:
            complete_microstructure = False
            logger.warning("Missing file for sheets. ")
        else:
            sheet_file = Path(sheet_file)

        if sheet_normal_file is None:
            complete_microstructure = False
            logger.warning("Missing file for sheet_normals. ")
        else:
            sheet_normal_file = Path(sheet_normal_file)

        if not complete_microstructure:
            # TODO: Could also generate fibers using ldrb
            raise RuntimeError("Missing microstructure")

        f0, s0, n0 = cls.load_microstructure(
            mesh,
            fiber_file,
            sheet_file,
            sheet_normal_file,
        )
        return cls(
            mesh=mesh,
            ffun=ffun,
            markers=markers,
            f0=f0,
            s0=s0,
            n0=n0,
        )

    @staticmethod
    def load_microstructure(mesh, fiber_file, sheet_file, sheet_normal_file):
        # Load signatures on root process
        signatures = {}
        if mesh.mpi_comm().rank == 0:
            with h5py.File(fiber_file, "r") as f:
                signatures["fiber"] = eval(f["fiber"].attrs["signature"].decode())
            with h5py.File(sheet_file, "r") as f:
                signatures["sheet"] = eval(f["sheet"].attrs["signature"].decode())
            with h5py.File(sheet_normal_file, "r") as f:
                signatures["sheet_normal"] = eval(
                    f["sheet_normal"].attrs["signature"].decode(),
                )
        # Broadcast to the other processes
        signatures = mesh.mpi_comm().bcast(signatures, root=0)
        signatures = signatures

        # Create functions and load them
        f0 = dolfin.Function(dolfin.FunctionSpace(mesh, signatures["fiber"]))
        with dolfin.HDF5File(
            mesh.mpi_comm(),
            fiber_file.as_posix(),
            "r",
        ) as f:
            f.read(f0, "fiber")

        s0 = dolfin.Function(dolfin.FunctionSpace(mesh, signatures["sheet"]))
        with dolfin.HDF5File(
            mesh.mpi_comm(),
            sheet_file.as_posix(),
            "r",
        ) as f:
            f.read(s0, "sheet")

        n0 = dolfin.Function(
            dolfin.FunctionSpace(mesh, signatures["sheet_normal"]),
        )
        with dolfin.HDF5File(
            mesh.mpi_comm(),
            sheet_normal_file.as_posix(),
            "r",
        ) as f:
            f.read(n0, "sheet_normal")

        return (f0, s0, n0)

    def save(self, fname: Union[str, Path]) -> None:
        """Save geometry to file

        Parameters
        ----------
        fname : str
            Name of file
        """
        save_geometry(fname=fname, geo=self)


class LVGeometry(HeartGeometry):
    r"""
    Create a truncated ellipsoidal geometry,
    defined through the coordinates:

    .. math::
        x &= r_l(t) \cos(\mu)\\
        y &= r_s(t) \sin(\mu) \cos(\theta)\\
        z &= r_s(t) \sin(\mu) \sin(\theta)\\

    for t in [0, 1], mu in [0, \mu_{mathrm{base}] and theta in [0, 2\pi).
    """

    @staticmethod
    def load_mesh(mesh_file):
        mesh = dolfin.Mesh()
        ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
        with dolfin.XDMFFile(mesh.mpi_comm(), mesh_file.as_posix()) as f:
            f.read(mesh)
            f.read(ffun_val)
        ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)

        # Hardcode markers from file
        markers = {
            "BASE": (3, 2),
            "ENDO": (1, 2),
            "EPI": (2, 2),
        }
        return mesh, ffun, markers

    @classmethod
    def from_parameters(
        cls,
        mesh_params: Optional[Dict[str, float]] = None,
        fiber_params: Optional[Dict[str, Union[float, str]]] = None,
    ) -> "LVGeometry":
        """Load geometry from parameters

        Parameters
        ----------
        mesh_params : Optional[Dict[str, float]], optional
            Parameters for the mesh, by default None
        fiber_params : Optional[Dict[str, Union[float, str]]], optional
            Parameters fro the fibers, by default None

        Returns
        -------
        LVGeometry
            The geometry
        """
        mesh_params = mesh_params or {}
        mesh_parameters = LVGeometry.default_mesh_parameters()
        mesh_parameters.update(mesh_params)

        fiber_params = fiber_params or {}
        fiber_parameters = LVGeometry.default_fiber_parameters()
        fiber_parameters.update(fiber_params)

        obj = cls()
        msh_name = Path("test.msh")
        create_benchmark_ellipsoid_mesh_gmsh(msh_name, **mesh_parameters)
        geo = gmsh2dolfin(msh_name)
        for key, value in geo._asdict().items():
            setattr(obj, key, value)

        msh_name.unlink()
        microstructure = create_microstructure(
            mesh=geo.mesh,
            ffun=geo.ffun,
            markers=geo.markers,
            mesh_params=mesh_parameters,
            fiber_params=fiber_parameters,
        )
        for key, value in microstructure._asdict().items():
            setattr(obj, key, value)
        return obj

    def __repr__(self):
        return f"{self.__class__.__name__}()"

    @staticmethod
    def default_mesh_parameters() -> Dict[str, float]:
        return dict(
            r_short_endo=0.025,
            r_short_epi=0.035,
            r_long_endo=0.09,
            r_long_epi=0.097,
            psize_ref=0.005,
            mu_apex_endo=-math.pi,
            mu_base_endo=-math.acos(5 / 17),
            mu_apex_epi=-math.pi,
            mu_base_epi=-math.acos(5 / 20),
            mesh_size_factor=1.0,
        )

    @staticmethod
    def default_fiber_parameters() -> Dict[str, Union[float, str]]:
        return dict(
            function_space="P_2",
            alpha_endo=-60.0,
            alpha_epi=+60.0,
        )


class BiVGeometry(HeartGeometry):
    @staticmethod
    def generate_fibers_ldrb(markers, mesh, ffun):
        import ldrb

        system = ldrb.dolfin_ldrb(
            mesh=mesh,
            fiber_space="P_2",
            ffun=ffun,
            markers=dict(
                base=markers["BASE"][0],
                rv=markers["ENDO_RV"][0],
                lv=markers["ENDO_LV"][0],
                epi=markers["EPI"][0],
            ),
            alpha_endo_lv=60,
            alpha_epi_lv=-60,
        )
        f0 = system.fiber
        s0 = system.sheet
        n0 = system.sheet_normal
        return (f0, s0, n0)

    @staticmethod
    def load_mesh(mesh_file):
        mesh = dolfin.Mesh()
        ffun_val = dolfin.MeshValueCollection("size_t", mesh, 2)
        with dolfin.XDMFFile(mesh.mpi_comm(), mesh_file.as_posix()) as f:
            f.read(mesh)
            f.read(ffun_val)
        ffun = dolfin.MeshFunction("size_t", mesh, ffun_val)

        # Hardcode markers from file
        markers = {
            "BASE": (10, 2),
            "ENDO_RV": (20, 2),
            "ENDO_LV": (30, 2),
            "EPI": (40, 2),
        }
        return mesh, ffun, markers


if __name__ == "__main__":
    geo = LVGeometry.from_parameters()
    dolfin.File("ffun.pvd") << geo.ffun
