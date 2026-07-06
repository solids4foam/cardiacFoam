from pathlib import Path

import pytest

from openfoam_driver.dashboard.assets3d import is_final_3d, export_glb


def test_allowlist():
    assert is_final_3d("PATHOS/RBBB")
    assert is_final_3d("heartSim3D-1D/monodomainHeart")
    assert not is_final_3d("manufacturedSolutions/monodomainPseudoECG")


def test_export_glb_from_synthetic_surface(tmp_path: Path):
    pv = pytest.importorskip("pyvista")
    src = tmp_path / "surf.vtp"
    sphere = pv.Sphere()
    sphere.save(src)
    out = tmp_path / "out.glb"
    result = export_glb(src, out)
    assert result == out
    assert out.is_file() and out.stat().st_size > 0


def test_export_glb_is_cached(tmp_path: Path):
    pv = pytest.importorskip("pyvista")
    src = tmp_path / "surf.vtp"
    pv.Sphere().save(src)
    out = tmp_path / "out.glb"
    export_glb(src, out)
    first_mtime = out.stat().st_mtime_ns
    export_glb(src, out)  # cached: unchanged source -> no rewrite
    assert out.stat().st_mtime_ns == first_mtime
