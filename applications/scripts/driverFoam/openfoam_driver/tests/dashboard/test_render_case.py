from pathlib import Path
from unittest import mock

import pytest

from openfoam_driver.dashboard import render_case


def test_find_blender_prefers_explicit(tmp_path):
    fake = tmp_path / "blender"
    fake.write_text("")
    fake.chmod(0o755)
    assert render_case.find_blender(str(fake)) == str(fake)


def test_find_blender_missing_raises():
    with mock.patch.object(render_case, "_CANDIDATES", []):
        with pytest.raises(FileNotFoundError):
            render_case.find_blender(None)


def test_render_invokes_blender_subprocess(tutorial_tree, tmp_path):
    case = tutorial_tree / "PATHOS" / "RBBB"
    (case / "heart.vtu").write_text("x")
    out = tmp_path / "renders"
    with mock.patch.object(render_case.subprocess, "run") as run, \
         mock.patch.object(render_case, "find_blender", return_value="/fake/blender"), \
         mock.patch.object(render_case, "_to_stl", return_value=out / "surface.stl"):
        run.return_value = mock.Mock(returncode=0)
        render_case.render(case_dir=case, out_dir=out, blender=None)
    assert run.called
    args = run.call_args[0][0]
    assert args[0] == "/fake/blender"
    assert "--background" in args


def test_render_missing_surface_raises(tutorial_tree, tmp_path):
    case = tutorial_tree / "PATHOS" / "RBBB"  # no vtu/vtk
    with mock.patch.object(render_case, "find_blender", return_value="/fake/blender"):
        with pytest.raises(FileNotFoundError):
            render_case.render(case_dir=case, out_dir=tmp_path, blender=None)
