import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import distort_mesh  # noqa: E402


def _write_synthetic_points(case_dir: Path, n: int) -> np.ndarray:
    lin = np.linspace(0.0, 1.0, n + 1)
    pts = np.array([[x, y, z] for x in lin for y in lin for z in lin])
    header = (
        "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
        "    class       vectorField;\n    location    \"constant/polyMesh\";\n"
        "    object      points;\n}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
    )
    body = "\n".join(f"({x:.17g} {y:.17g} {z:.17g})" for x, y, z in pts)
    footer = "\n\n\n// ************************************************************************* //\n"
    text = f"{header}{pts.shape[0]}\n(\n{body}\n)\n{footer}"
    polymesh = case_dir / "constant" / "polyMesh"
    polymesh.mkdir(parents=True, exist_ok=True)
    (polymesh / "points").write_text(text)
    return pts


def test_zero_amplitude_is_identity(tmp_path):
    original = _write_synthetic_points(tmp_path, n=2)
    points_path = tmp_path / "constant" / "polyMesh" / "points"

    _, _, coords = distort_mesh.read_points(points_path)
    np.testing.assert_allclose(coords, original)

    displaced = distort_mesh.distort(coords, amplitude=0.0, h=0.5, dimension=3, tol=1e-9)
    np.testing.assert_array_equal(displaced, coords)


def test_boundary_points_are_untouched(tmp_path):
    _write_synthetic_points(tmp_path, n=2)
    points_path = tmp_path / "constant" / "polyMesh" / "points"
    _, _, coords = distort_mesh.read_points(points_path)

    displaced = distort_mesh.distort(coords, amplitude=0.15, h=0.5, dimension=3, tol=1e-9)

    on_boundary = np.any((coords == 0.0) | (coords == 1.0), axis=1)
    np.testing.assert_allclose(displaced[on_boundary], coords[on_boundary])


def test_interior_point_matches_closed_form(tmp_path):
    _write_synthetic_points(tmp_path, n=2)
    points_path = tmp_path / "constant" / "polyMesh" / "points"
    _, _, coords = distort_mesh.read_points(points_path)

    displaced = distort_mesh.distort(coords, amplitude=0.15, h=0.5, dimension=3, tol=1e-9)

    center = np.isclose(coords, 0.5).all(axis=1)
    assert center.sum() == 1
    delta = displaced[center][0] - coords[center][0]
    expected_dx = 0.15 * 0.5 * np.sin(np.pi * 0.5) ** 3
    assert delta[0] == pytest.approx(expected_dx, abs=1e-12)
    assert delta[1] == pytest.approx(0.0, abs=1e-12)
    assert delta[2] == pytest.approx(0.0, abs=1e-12)


def test_binary_points_file_is_rejected(tmp_path):
    polymesh = tmp_path / "constant" / "polyMesh"
    polymesh.mkdir(parents=True)
    (polymesh / "points").write_text(
        "FoamFile\n{\n    format      binary;\n    class       vectorField;\n}\n"
    )

    with pytest.raises(SystemExit):
        distort_mesh.read_points(polymesh / "points")


def test_2d_never_moves_z(tmp_path):
    _write_synthetic_points(tmp_path, n=4)
    points_path = tmp_path / "constant" / "polyMesh" / "points"
    _, _, coords = distort_mesh.read_points(points_path)

    displaced = distort_mesh.distort(coords, amplitude=0.15, h=0.25, dimension=2, tol=1e-9)

    np.testing.assert_array_equal(displaced[:, 2], coords[:, 2])
