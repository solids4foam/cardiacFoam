import shutil
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from typer.testing import CliRunner

from cardiac_benchmark.cli import app

here = Path(__file__).parent.absolute()


@pytest.fixture(scope="module")
def geo_path():
    yield (here / "test_geometry.h5").absolute().as_posix()


def _benchmark_1(
    fname: str,
    cmd: str,
    _geo_path: str,
    max_up0: float,
    min_up0: float,
) -> None:
    runner = CliRunner(mix_stderr=False)
    outdir = Path(f"test_{fname}")
    if outdir.is_dir():
        shutil.rmtree(outdir)

    # Use P_1 to speed up the tests
    result = runner.invoke(
        app,
        [
            *cmd,
            "--geometry-path",
            _geo_path,
            "--outdir",
            outdir.as_posix(),
            "--function-space=P_1",
        ],
    )

    assert result.exit_code == 0

    up0 = np.load(outdir / "componentwise_displacement_up0.npy")
    time_stamps = np.load(outdir / "time_stamps.npy")

    assert np.isclose(time_stamps[0], 0.001)
    assert np.isclose(time_stamps[-1], 0.999)

    computed_max = up0.max(0)
    computed_min = up0.min(0)

    assert np.allclose(computed_max, max_up0)
    assert np.allclose(computed_min, min_up0)

    shutil.rmtree(outdir)


def test_benchmark1_step0_caseA(geo_path):
    _benchmark_1(
        "benchmark1_step0_caseA",
        ["benchmark1-step0", "a", "--geometry-path", geo_path],
        geo_path,
        max_up0=(1.34564644e-05, 1.13097351e-03, 1.17750664e-03),
        min_up0=(-2.10341755e-02, 3.73614057e-07, -3.95987517e-04),
    )


def test_benchmark1_step0_caseB(geo_path):
    _benchmark_1(
        "benchmark1_step0_caseB",
        ["benchmark1-step0", "b", "--geometry-path", geo_path],
        geo_path,
        max_up0=(1.00143349e-02, 4.97783144e-05, 1.64292843e-04),
        min_up0=(-1.00857813e-05, -2.82413565e-04, -2.64611270e-07),
    )


def test_benchmark1_step1(geo_path):
    _benchmark_1(
        "benchmark1_step1",
        ["benchmark1-step1", "--geometry-path", geo_path],
        geo_path,
        max_up0=(0.00039218, 0.00151579, 0.00122652),
        min_up0=(-1.98295540e-02, 4.51570134e-07, -3.08690114e-04),
    )


def test_benchmark1_step2_caseA(geo_path):
    _benchmark_1(
        "benchmark1_step2_caseA",
        ["benchmark1-step2", "a", "--geometry-path", geo_path],
        geo_path,
        max_up0=(0.00031393, 0.00156036, 0.00175083),
        min_up0=(-2.10530866e-02, -2.32237185e-07, -2.32820346e-04),
    )


def test_benchmark2():
    runner = CliRunner(mix_stderr=False)
    data_folder = Path.cwd() / "coarse_data"
    res_download = runner.invoke(
        app,
        ["download-data-benchmark2", "coarse", "--outdir", data_folder.as_posix()],
    )
    assert res_download.exit_code == 0

    geometry_path = data_folder / "biv_geometry_coarse.h5"

    res_convert = runner.invoke(
        app,
        [
            "convert-data-benchmark2",
            data_folder.as_posix(),
            "--outpath",
            geometry_path.as_posix(),
        ],
    )
    assert res_convert.exit_code == 0

    outdir = Path.cwd() / "test_results_benchmark2"
    # Let us mock out the actual solving to speed it up
    with patch("cardiac_benchmark.solver.NonlinearSolver.solve") as mock_solve:
        mock_solve.return_value = (0, True)
        result = runner.invoke(
            app,
            [
                "benchmark2",
                "--geometry-path",
                geometry_path.as_posix(),
                "--outdir",
                outdir.as_posix(),
                "--t",
                0.05,
            ],
        )
    assert result.exit_code == 0

    assert (outdir / "result.h5").is_file()
    assert (outdir / "von_Mises_stress_sp1.npy").is_file()
