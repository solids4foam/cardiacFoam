import pytest

import cardiac_benchmark


@pytest.fixture
def lvgeo():
    return cardiac_benchmark.geometry.LVGeometry.from_parameters()


@pytest.fixture
def lvproblem(lvgeo):
    material = cardiac_benchmark.material.HolzapfelOgden(
        f0=lvgeo.f0,
        s0=lvgeo.s0,
    )
    return cardiac_benchmark.problem.LVProblem(
        geometry=lvgeo,
        material=material,
    )


def test_datacollector_lv(tmp_path, lvproblem):
    result_file = tmp_path / "result.h5"

    collector = cardiac_benchmark.postprocess.DataCollector(
        result_file,
        problem=lvproblem,
    )
    collector.store(0.0)
    loader = cardiac_benchmark.postprocess.DataLoader(result_file)

    # Check mesh
    assert (
        loader.problem.geometry.mesh.coordinates()
        == lvproblem.geometry.mesh.coordinates()
    ).all()
    # and check problem type
    assert isinstance(loader.problem, type(lvproblem))


def test_datacollector_biv(tmp_path, lvproblem):
    result_file = tmp_path / "result.h5"

    collector = cardiac_benchmark.postprocess.DataCollector(
        result_file,
        problem=lvproblem,
    )
    collector.store(0.0)
    loader = cardiac_benchmark.postprocess.DataLoader(result_file)

    # Check mesh
    assert (
        loader.problem.geometry.mesh.coordinates()
        == lvproblem.geometry.mesh.coordinates()
    ).all()
    # and check problem type
    assert isinstance(loader.problem, type(lvproblem))
