"""The continuum junction on a mesh (mms_common/pvj_continuum.py, emitted as
reference/pvj_continuum.py): kernel masses, support integrals, cell integrals against analytic
values, and (cardiacFoam checkout with OpenFOAM) parity with the C++ class
pvjContinuumReference."""

import ast
import importlib.util
import os
import tempfile

import numpy as np
import pytest

import coupled_1d3d_mms as cp
from mms_common import pvj_continuum as pc
from mms_common import pvj_meshes as meshes

R = 0.3
KERNELS = pc.KERNELS
# the kernel profiles of KERNEL_FAMILY by cardiacFoam name
PHI = {key: phi for _, key, phi in cp.KERNEL_FAMILY}


def _reference():
    path = os.path.join(cp.HERE, "reference", "pvj_continuum.py")
    spec = importlib.util.spec_from_file_location("pvj_continuum_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _radial_moment(kernel, power, upper=1.0):
    """int_0^1 phi(s) s^power ds by Gauss-Legendre (the truncated Gaussian is smooth on [0, 1])."""
    x, w = np.polynomial.legendre.leggauss(64)
    s = 0.5*upper*(x + 1.0)
    return float(np.sum(0.5*upper*w*pc.kernel_profile(kernel, s)*s**power))


def _hex_cells(n, lo, hi, jitter=0.0, seed=0):
    points, hexes = meshes.hexahedra(meshes.jitter_inner(meshes.structured_points(n, lo, hi), jitter, seed)
                                     if jitter else meshes.structured_points(n, lo, hi))
    return pc.cells_from_hexahedra(points[hexes])


def _tet_cells(n, lo, hi, jitter=0.0, seed=0):
    grid = meshes.structured_points(n, lo, hi)
    points, tets = meshes.kuhn_tetrahedra(meshes.jitter_inner(grid, jitter, seed) if jitter else grid)
    return pc.cells_from_tetrahedra(points[tets])


# ---- the emitted module ---------------------------------------------------------------------

def test_emitted_files_up_to_date():
    assert cp.outputs_match()


def test_emitted_reference_is_numpy_only():
    text = open(os.path.join(cp.HERE, "reference", "pvj_continuum.py")).read()
    tree = ast.parse(text)
    imported = [alias.name for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom))
                for alias in node.names]
    assert sorted(imported) == ["collections", "numpy"]
    assert not any(isinstance(node, ast.ImportFrom) for node in ast.walk(tree))
    for banned in ("sympy", "scipy", "np."):
        assert banned not in text


def test_emitted_reference_equals_module():
    ref = _reference()
    cells = _hex_cells(3, (-0.5,)*3, (0.5,)*3, jitter=0.2, seed=4)
    t = (0.02, -0.03, 0.01)
    f = lambda p: 1 + p[:, 0] - 2*p[:, 1]**2
    a = ref.cell_integrals(cells, t, R, "gaussian", f, cut_size=0.125)
    b = pc.cell_integrals(cells, t, R, "gaussian", f, cut_size=0.125)
    for u, v in zip(a, b):
        np.testing.assert_allclose(u, v, rtol=1e-14, atol=0)
    assert ref.support_integral(t, R, "cosine", f, "interior") == pc.support_integral(t, R, "cosine", f)


# ---- kernels and the support integral --------------------------------------------------------

def test_kernel_profiles_and_weights():
    s = np.array([0.0, 0.25, 0.5, 1.0])
    np.testing.assert_allclose(pc.kernel_profile("uniform", s), 1.0)
    np.testing.assert_allclose(pc.kernel_profile("linear", s), 1 - s)
    np.testing.assert_allclose(pc.kernel_profile("gaussian", s), np.exp(-4.5*s**2))
    np.testing.assert_allclose(pc.kernel_profile("cosine", s), 0.5*(1 + np.cos(np.pi*s)))
    r = np.array([0.1, 0.5, 1.0, 1.0001])
    for kernel in KERNELS:
        w = pc.kernel_weight(kernel, r, 1.0)
        assert w[-1] == 0.0 and w[-2] == pc.kernel_profile(kernel, 1.0)
    with pytest.raises(ValueError, match="uniform linear gaussian cosine"):
        pc.kernel_profile("top-hat", 0.5)


def test_kernel_profiles_follow_the_documented_family():
    s = np.linspace(0, 1, 11)
    import sympy as sp
    for _, key, phi in cp.KERNEL_FAMILY:
        np.testing.assert_allclose(pc.kernel_profile(key, s),
                                   sp.lambdify(cp._s, phi, "numpy")(s)*np.ones_like(s), rtol=1e-15)


@pytest.mark.parametrize("kernel", KERNELS)
def test_interior_and_plane_mass_are_the_documented_mass(kernel):
    half = float(cp.half_ball_mass_factor(PHI[kernel]))*R**3
    t = (0.3, -0.2, 0.1)
    assert pc.support_integral(t, R, kernel) == pytest.approx(2*half, rel=1e-13)
    for normal in ((1, 0, 0), (0, -1, 0), (0.3, 0.5, -0.8)):
        assert pc.support_integral(t, R, kernel, None, "plane", normal) == pytest.approx(half, rel=1e-13)


@pytest.mark.parametrize("kernel", KERNELS)
@pytest.mark.parametrize("wall", [0.5, 2.0])
def test_cavity_sphere_mass(kernel, wall):
    # tissue outside the sphere: the solid angle at distance r is 2 pi (1 + r/(2 rho))
    expected = 2*np.pi*R**3*(_radial_moment(kernel, 2) + R/(2*wall)*_radial_moment(kernel, 3))
    for normal in ((1.0, 0.0, 0.0), (0.4, -0.7, 0.2)):
        got = pc.support_integral((0.1, 0.2, 0.3), R, kernel, None, "cavitySphere", normal, wall)
        assert got == pytest.approx(expected, rel=1e-13)


@pytest.mark.parametrize("kernel", KERNELS)
def test_support_integral_of_fields(kernel):
    t = np.array([0.1, 0.2, 0.3])
    n = np.array([0.4, -0.7, 0.2])
    n = n/np.linalg.norm(n)
    along = lambda p: (p - t) @ n
    # interior: an odd field integrates to zero, |x - t|^2 to 4 pi R^5 int phi s^4 ds
    assert abs(pc.support_integral(t, R, kernel, along)) < 1e-15
    assert pc.support_integral(t, R, kernel, lambda p: np.sum((p - t)**2, axis=1)) == pytest.approx(
        4*np.pi*R**5*_radial_moment(kernel, 4), rel=1e-13)
    # plane: int mu dOmega over the half sphere is pi
    assert pc.support_integral(t, R, kernel, along, "plane", n) == pytest.approx(
        np.pi*R**4*_radial_moment(kernel, 3), rel=1e-13)
    # cavity sphere: mu from -r/(2 rho): int mu dOmega = pi (1 - (r/(2 rho))^2)
    wall = 0.5
    expected = np.pi*R**4*(_radial_moment(kernel, 3) - (R/(2*wall))**2*_radial_moment(kernel, 5))
    assert pc.support_integral(t, R, kernel, along, "cavitySphere", n, wall) == pytest.approx(
        expected, rel=1e-13)


def test_support_integral_rejects_bad_input():
    t = (0, 0, 0)
    with pytest.raises(ValueError, match="kind"):
        pc.support_integral(t, R, "linear", kind="cone")
    with pytest.raises(ValueError, match="normal"):
        pc.support_integral(t, R, "linear", kind="plane")
    with pytest.raises(ValueError, match="wall_radius"):
        pc.support_integral(t, R, "linear", kind="cavitySphere", normal=(1, 0, 0))
    with pytest.raises(ValueError, match="kernel"):
        pc.support_integral(t, R, "hat")
    with pytest.raises(ValueError, match="positive"):
        pc.cell_integrals(_hex_cells(1, (0,)*3, (1,)*3), t, 0.0)


def test_support_integral_converges_with_the_rule():
    t, n = (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    coarse = pc.support_integral(t, R, "gaussian", None, "cavitySphere", n, 0.4, 24, 24, 48)
    fine = pc.support_integral(t, R, "gaussian", None, "cavitySphere", n, 0.4, 64, 64, 128)
    assert abs(coarse - fine) < 1e-9*fine


# ---- the cell geometry -------------------------------------------------------------------------

def test_cell_geometry_of_a_box():
    cell = pc.cells_from_hexahedra(np.array([[[0, 0, 0], [2, 0, 0], [2, 3, 0], [0, 3, 0],
                                              [0, 0, 4], [2, 0, 4], [2, 3, 4], [0, 3, 4]]], float))[0]
    centre, volume, face_centres = pc.cell_geometry(cell)
    np.testing.assert_allclose(centre, (1.0, 1.5, 2.0), atol=1e-15)
    assert volume == pytest.approx(24.0, rel=1e-15)
    assert len(face_centres) == 6
    # the tetrahedron: centroid of the vertices, volume |det|/6
    v = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0.5, 0.5, 3]], float)
    centre, volume, _ = pc.cell_geometry(pc.cells_from_tetrahedra(v[None])[0])
    np.testing.assert_allclose(centre, v.mean(axis=0), atol=1e-15)
    assert volume == pytest.approx(abs(np.linalg.det(v[1:] - v[0]))/6, rel=1e-14)


@pytest.mark.parametrize("builder", [_hex_cells, _tet_cells])
def test_warped_meshes_partition_the_box(builder):
    cells = builder(4, (-0.5,)*3, (0.5,)*3, jitter=0.3, seed=7)
    volumes = pc.cell_volumes(cells)
    assert np.all(volumes > 0)
    assert volumes.sum() == pytest.approx(1.0, rel=1e-13)


def test_curved_hexahedra_have_the_shell_volume():
    points, hexes = meshes.hexahedra(meshes.shell_sector_grid((6, 12, 12)))
    volumes = pc.cell_volumes(pc.cells_from_hexahedra(points[hexes]))
    assert np.all(volumes > 0)
    # the sector r1 <= r <= r2, theta0 <= theta <= pi - theta0, |phi| <= 0.7 has the exact volume
    exact = (1.0**3 - 0.5**3)/3*0.7*2*2*np.cos(1.0)
    assert volumes.sum() == pytest.approx(exact, rel=0.02)     # the facets cut the sphere


def test_polymesh_cells_match_the_vertex_helpers():
    grid = meshes.jitter_inner(meshes.structured_points(3, (-0.5,)*3, (0.5,)*3), 0.25, seed=2)
    points, hexes = meshes.hexahedra(grid)
    faces, owner, neighbour, n_internal = meshes.polymesh_arrays(points, hexes)
    assert n_internal == 3*3*2*3 and len(neighbour) == n_internal
    assert len(faces) == n_internal + 6*9
    from_poly = pc.cells_from_polymesh(points, faces, owner, neighbour)
    from_hex = pc.cells_from_hexahedra(points[hexes])
    assert len(from_poly) == len(from_hex) == 27
    np.testing.assert_allclose(pc.cell_volumes(from_poly), pc.cell_volumes(from_hex), rtol=1e-14)
    for a, b in zip(from_poly, from_hex):
        np.testing.assert_allclose(pc.cell_geometry(a)[0], pc.cell_geometry(b)[0], atol=1e-14)
    # a cell inside the ball is integrated whole: the face order changes nothing visible
    inside = pc.cell_integrals(from_poly, (0.0, 0.0, 0.0), 5.0, "cosine")
    again = pc.cell_integrals(from_hex, (0.0, 0.0, 0.0), 5.0, "cosine")
    np.testing.assert_allclose(inside[1], again[1], rtol=1e-12)


def test_face_orientation_changes_only_the_quadrature_error():
    cells = _hex_cells(3, (-0.5,)*3, (0.5,)*3, jitter=0.2, seed=1)
    reversed_cells = [[f[::-1] for f in cell] for cell in cells]
    t = (0.02, -0.03, 0.01)
    a = pc.cell_integrals(cells, t, R, "cosine", cut_size=0.125)
    b = pc.cell_integrals(reversed_cells, t, R, "cosine", cut_size=0.125)
    np.testing.assert_allclose(pc.cell_volumes(cells), pc.cell_volumes(reversed_cells), rtol=1e-14)
    assert np.array_equal(a[0], b[0])
    assert abs(a[1].sum() - b[1].sum()) < 1e-8*a[1].sum()       # 1.8e-9 at this cut size
    assert np.max(np.abs(a[1] - b[1])) < 1e-6*a[1].max()


# ---- the cell integrals against analytic values ----------------------------------------------------

CUT = {"uniform": 1/12, "gaussian": 1/12, "linear": 1/12, "cosine": 1/12}
# measured at cut size 1/12 on the jittered mesh below: 6e-4, 5e-5, 1.3e-6, 4.7e-6
TOLERANCE = {"uniform": 2e-3, "gaussian": 2e-4, "linear": 5e-5, "cosine": 5e-5}


@pytest.mark.parametrize("kernel", KERNELS)
def test_interior_cells_sum_to_the_support_integral(kernel):
    cells = _hex_cells(4, (-0.5,)*3, (0.5,)*3, jitter=0.2, seed=1)
    t = (0.02, -0.03, 0.01)
    f = lambda p: np.sum((p - np.asarray(t))**2, axis=1)
    index, omega, omega_f = pc.cell_integrals(cells, t, R, kernel, f, cut_size=CUT[kernel])
    assert len(index) == len(np.unique(index)) and np.all(np.diff(index) > 0)
    mass = pc.support_integral(t, R, kernel)
    assert omega.sum() == pytest.approx(mass, rel=TOLERANCE[kernel])
    assert omega_f.sum() == pytest.approx(pc.support_integral(t, R, kernel, f), rel=5*TOLERANCE[kernel])
    assert np.all(omega >= 0)


@pytest.mark.parametrize("kernel,low,high", [("uniform", 3.7, 4.3), ("linear", 6.5, 7.6)])
def test_cut_error_order(kernel, low, high):
    # clipping by the linear interpolant of the sphere: second order when the kernel does not
    # vanish at R (halving the cut size divides the error by 4), third order when it does
    # (linear, cosine: the weight is small where the interpolant is least accurate)
    cells = _hex_cells(4, (-0.5,)*3, (0.5,)*3, jitter=0.2, seed=1)
    t = (0.02, -0.03, 0.01)
    mass = pc.support_integral(t, R, kernel)
    err = [abs(pc.cell_integrals(cells, t, R, kernel, cut_size=c)[1].sum()/mass - 1)
           for c in (1/12, 1/24, 1/48)]
    assert low < err[0]/err[1] < high and low < err[1]/err[2] < high


@pytest.mark.parametrize("builder", [_hex_cells, _tet_cells])
@pytest.mark.parametrize("kernel", ["linear", "cosine"])
def test_plane_cells_sum_to_the_half_ball(builder, kernel):
    cells = builder(4, (0.0,)*3, (1.0,)*3, jitter=0.2, seed=2)
    t = (0.0, 0.5, 0.5)
    index, omega, _ = pc.cell_integrals(cells, t, R, kernel, cut_size=1/24)
    mass = pc.support_integral(t, R, kernel, None, "plane", (1, 0, 0))
    assert omega.sum() == pytest.approx(mass, rel=5e-6)


def test_gathered_value_of_a_linear_field_is_its_value_at_the_terminal():
    cells = _hex_cells(4, (-0.5,)*3, (0.5,)*3, jitter=0.2, seed=3)
    t = np.array([0.05, 0.0, -0.04])
    f = lambda p: 2.0 + 0.5*p[:, 0] - 1.5*p[:, 1] + 0.7*p[:, 2]
    result = pc.junction_source(cells, t, R, "linear", "interior", current=2.0, f=f, cut_size=1/12)
    assert result.omega_f.sum()/result.normalisation == pytest.approx(f(t[None])[0], rel=1e-5)


@pytest.mark.parametrize("builder", [_hex_cells, _tet_cells])
def test_cell_average_source_deposits_its_current(builder):
    cells = builder(4, (-0.5,)*3, (0.5,)*3, jitter=0.2, seed=5)
    t = (-0.04, 0.03, 0.02)
    out = pc.junction_source(cells, t, R, "linear", "interior", current=-1.7, cut_size=1/12)
    assert np.sum(out.source*out.volumes) == pytest.approx(-1.7, rel=5e-5)
    assert np.all(out.source < 0)
    assert out.volumes.shape == out.omega.shape == out.cells.shape
    # one cell holding the whole ball deposits I omega/(W |K|) = I/|K|
    one = pc.junction_source(pc.cells_from_hexahedra(
        np.array([[[0, 0, 0], [4, 0, 0], [4, 4, 0], [0, 4, 0], [0, 0, 4], [4, 0, 4], [4, 4, 4],
                   [0, 4, 4]]], float)), (2, 2, 2), 0.3, "cosine", current=3.0)
    assert one.source[0] == pytest.approx(3.0/64.0, rel=1e-4)


def test_shell_cells_approach_the_cavity_normalisation_at_second_order():
    # the polyhedral inner wall cuts the lens of the support: W of the true sphere is the limit
    r1, radius = 0.5, 0.2
    t = np.array([r1, 0.0, 0.0])
    mass = pc.support_integral(t, radius, "linear", None, "cavitySphere", (1, 0, 0), r1)
    err = []
    for n in ((3, 8, 8), (6, 16, 16)):
        points, hexes = meshes.hexahedra(meshes.shell_sector_grid(n))
        _, omega, _ = pc.cell_integrals(pc.cells_from_hexahedra(points[hexes]), t, radius, "linear",
                                        cut_size=1/24)
        err.append(omega.sum()/mass - 1)
    assert 0 < err[1] < err[0] < 0.03                   # the sphere lies inside its facets
    assert 3.0 < err[0]/err[1] < 5.0


def test_cells_far_from_the_ball_are_left_out():
    cells = _hex_cells(2, (0.0,)*3, (1.0,)*3)
    out = pc.cell_integrals(cells, (5.0, 5.0, 5.0), 0.3, "linear")
    assert all(len(a) == 0 for a in out)
    near = pc.cell_integrals(cells, (0.5, 0.5, 0.5), 0.1, "linear", cut_size=0.25)
    assert len(near[0]) == 8 and near[1].sum() == pytest.approx(
        pc.support_integral((0.5, 0.5, 0.5), 0.1, "linear"), rel=1e-4)


# ---- the support of a terminal --------------------------------------------------------------

def test_box_support():
    kind, normal = pc.box_support((0.0, 1/6, 1/3), (0,)*3, (1,)*3, 0.11)
    assert kind == "plane" and np.array_equal(normal, (1, 0, 0))
    kind, normal = pc.box_support((1.0, 1/6, 1/3), (0,)*3, (1,)*3, 0.11)
    assert kind == "plane" and np.array_equal(normal, (-1, 0, 0))
    assert pc.box_support((0.5, 0.5, 0.5), (0,)*3, (1,)*3, 0.11) == ("interior", pytest.approx(np.zeros(3)))
    with pytest.raises(ValueError, match="cut by a wall"):
        pc.box_support((0.05, 0.5, 0.5), (0,)*3, (1,)*3, 0.11)
    with pytest.raises(ValueError, match="edge or corner"):
        pc.box_support((0.0, 0.0, 0.5), (0,)*3, (1,)*3, 0.11)
    with pytest.raises(ValueError, match="opposite wall"):
        pc.box_support((0.0, 0.5, 0.5), (0,)*3, (0.1, 1, 1), 0.11)


def test_shell_support():
    kind, normal, wall = pc.shell_support((0.5, 0.0, 0.0), 0.5, 1.0, 0.2)
    assert kind == "cavitySphere" and np.allclose(normal, (1, 0, 0)) and wall == 0.5
    kind, normal, wall = pc.shell_support((0.0, 0.0, -0.5), 0.5, 1.0, 0.2)
    assert np.allclose(normal, (0, 0, -1))
    assert pc.shell_support((0.75, 0.0, 0.0), 0.5, 1.0, 0.2)[0] == "interior"
    with pytest.raises(ValueError, match="neither on the inner sphere"):
        pc.shell_support((0.55, 0.0, 0.0), 0.5, 1.0, 0.2)
    with pytest.raises(ValueError, match="outer sphere"):
        pc.shell_support((0.5, 0.0, 0.0), 0.5, 0.6, 0.2)
    with pytest.raises(ValueError, match="neither on the inner sphere"):
        pc.shell_support((1.0, 0.0, 0.0), 0.5, 1.0, 0.2)         # on the outer sphere: not covered


def test_self_checks_and_cli():
    assert cp.main(["--self-check"]) == 0


# ---- parity with the C++ class pvjContinuumReference ---------------------------------------------

def _parity():
    try:
        from parity import pvj_parity
    except ImportError:
        pytest.skip("parity harness absent (cardiacFoam-only)")
    ok, why = pvj_parity.parity_available()
    if not ok:
        pytest.skip(why)
    return pvj_parity


FIELDS = {
    "cosine": lambda p: np.cos(np.pi*p[:, 0])*np.cos(2*np.pi*p[:, 1])*np.cos(3*np.pi*p[:, 2]),
    "poly": lambda p: (1 + 0.3*p[:, 0] - 0.7*p[:, 1] + 0.2*p[:, 2] + p[:, 0]*p[:, 1]
                       - 0.5*p[:, 2]**2),
}


def _compare_with_cpp(points, cells, kernel, radius, terminal, support, field="poly",
                      numpy_support=None, **cpp_options):
    """Run both on one mesh; asserts agreement and returns the worst relative differences of
    (geometry, support integrals, cell integrals)."""
    pvj_parity = _parity()
    options = {"cut_size": cpp_options.get("referenceCutSize"),
               "inside_size": cpp_options.get("referenceInsideSize", pc.INSIDE_SIZE),
               "tip_size": cpp_options.get("referenceTipSize", pc.TIP_SIZE)}
    with tempfile.TemporaryDirectory() as case:
        faces, owner, neighbour = pvj_parity.write_case(
            case, points, cells, kernel, radius, terminal, support, field, **cpp_options)
        volumes, centres, (w_cpp, wf_cpp), index, omega, omega_f = pvj_parity.run_case(case)
    mesh = pc.cells_from_polymesh(points, faces, owner, neighbour)
    f = FIELDS[field]
    geometry = np.array([pc.cell_geometry(c)[:2] for c in mesh], dtype=object)
    d_volume = np.max(np.abs(pc.cell_volumes(mesh)/volumes - 1))
    d_centre = np.max(np.abs(np.array([g[0] for g in geometry]) - centres))
    index_np, omega_np, omega_f_np = pc.cell_integrals(mesh, terminal, radius, kernel, f, **options)
    assert np.array_equal(index, index_np), "the cells that meet the ball differ"
    kwargs = dict(kind=support["kind"], normal=support.get("normal"), wall_radius=support.get("wallRadius"))
    w_np = pc.support_integral(terminal, radius, kernel, None, **kwargs)
    wf_np = pc.support_integral(terminal, radius, kernel, f, **kwargs)
    assert d_volume < 1e-13 and d_centre < 1e-13
    assert abs(w_np/w_cpp - 1) < 1e-13 and abs(wf_np - wf_cpp) < 1e-13*max(abs(wf_cpp), w_cpp)
    scale = omega.max()
    assert np.max(np.abs(omega_np - omega)) < 1e-12*scale
    assert np.max(np.abs(omega_f_np - omega_f)) < 1e-12*max(np.abs(omega_f).max(), scale)
    return d_volume, d_centre, np.max(np.abs(omega_np - omega))/scale


HEX_GRID = lambda: meshes.jitter_inner(meshes.structured_points(3, (-0.5,)*3, (0.5,)*3), 0.2, seed=1)


@pytest.mark.parametrize("kernel,cut", [("uniform", 0.125), ("linear", None), ("gaussian", 0.125),
                                        ("cosine", 0.0625)])
def test_parity_interior_support_on_a_warped_hexahedral_mesh(kernel, cut):
    points, hexes = meshes.hexahedra(HEX_GRID())
    options = {"referenceCutSize": cut} if cut else {}
    _compare_with_cpp(points, hexes, kernel, R, (0.02, -0.03, 0.01), {"kind": "interior"},
                      **options)


@pytest.mark.parametrize("kernel", ["linear", "cosine"])
def test_parity_plane_support_on_a_warped_tetrahedral_mesh(kernel):
    grid = meshes.jitter_inner(meshes.structured_points(3), 0.2, seed=3)
    points, tets = meshes.kuhn_tetrahedra(grid)
    # Kuhn tetrahedra of a grid are symmetric: the octahedron diagonals of their sub-tetrahedra
    # tie exactly, and a tie is broken by the last bit of a distance, which differs between
    # compilers (a weakly cut cell then differs by 1e-6).  Moving every point by 1e-9 removes
    # the ties; the two implementations then agree to 1e-14 in every cell.
    points = points + np.random.default_rng(5).normal(scale=1e-9, size=points.shape)
    _compare_with_cpp(points, tets, kernel, 0.2573, (0.0, 0.5, 0.5),
                      {"kind": "plane", "normal": (1, 0, 0)}, referenceCutSize=0.0625)


@pytest.mark.parametrize("tetrahedra", [False, True])
def test_parity_cavity_sphere_support_on_a_curved_shell_mesh(tetrahedra):
    grid = meshes.shell_sector_grid((3, 6, 6), jitter=0.15, seed=1)
    points, cells = (meshes.kuhn_tetrahedra if tetrahedra else meshes.hexahedra)(grid)
    _compare_with_cpp(points, cells, "cosine", 0.2, (0.5, 0.0, 0.0),
                      {"kind": "cavitySphere", "normal": (1, 0, 0), "wallRadius": 0.5},
                      field="cosine", referenceCutSize=0.0625)


def test_parity_kernel_independent_quadrature_options():
    # the sizes of the C++ dictionary are the keyword arguments of cell_integrals
    points, hexes = meshes.hexahedra(HEX_GRID())
    _compare_with_cpp(points, hexes, "linear", R, (0.02, -0.03, 0.01), {"kind": "interior"},
                      referenceCutSize=0.1, referenceInsideSize=0.4, referenceTipSize=0.05)


def test_parity_of_the_support_integral_rule_options():
    pvj_parity = _parity()
    t = (0.1, 0.2, 0.3)
    n = (0.4, -0.7, 0.2)
    points, hexes = meshes.hexahedra(meshes.structured_points(1))
    with tempfile.TemporaryDirectory() as case:
        pvj_parity.write_case(case, points, hexes, "gaussian", 0.3, t,
                              {"kind": "cavitySphere", "normal": n, "wallRadius": 0.45},
                              referenceRadialPoints=20, referencePolarPoints=16,
                              referenceAzimuthalPoints=24)
        _, _, (w_cpp, _), _, _, _ = pvj_parity.run_case(case)
    assert pc.support_integral(t, 0.3, "gaussian", None, "cavitySphere", n, 0.45, 20, 16, 24) == \
        pytest.approx(w_cpp, rel=1e-13)
