"""Release-level tests: the figure scripts."""

import os
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))

ANALYTIC_FIGURES = {
    "monodomain": ["fda_cube.png", "rotated_cube.png", "shell_monodomain.png",
                   "pseudo_ecg.png", "graph_1d3d_geometry.png", "graph_1d3d_fields.png"],
    "eikonal": ["cartesian_eikonal.png", "eikonal_ecg.png", "shell_eikonal.png",
                "shell_eikonal_c0.png"],
}


def _run(mode, outdir, *extra):
    script = os.path.join(HERE, "docs", mode, "figures.py")
    return subprocess.run([sys.executable, script, "--draft", "--outdir", str(outdir),
                           *extra], capture_output=True, text=True, cwd=HERE)


@pytest.mark.parametrize("mode", sorted(ANALYTIC_FIGURES))
def test_figures_build_without_mesh_cases(mode, tmp_path):
    result = _run(mode, tmp_path)
    assert result.returncode == 0, result.stderr
    for name in ANALYTIC_FIGURES[mode]:
        path = tmp_path / name
        assert path.exists() and path.stat().st_size > 5000, name
    assert "--mesh-case" in result.stdout


def test_mesh_case_without_polymesh_fails_loudly(tmp_path):
    empty = tmp_path / "case"
    empty.mkdir()
    result = _run("monodomain", tmp_path / "out", "--mesh-case", "fda_cube_hex=%s" % empty)
    assert result.returncode == 1
    assert "fda_cube_hex" in result.stderr and "polyMesh" in result.stderr
