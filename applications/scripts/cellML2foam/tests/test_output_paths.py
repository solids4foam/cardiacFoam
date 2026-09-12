import subprocess
import sys
from pathlib import Path

import pytest


UTILITY_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = UTILITY_DIR / "src"
sys.path.insert(0, str(SRC_DIR))

import mapping_engine
import pipeline


SIM_C = """
void updateConstants(void) {
    AC_scale = 1.0;
}

void default_initial_values(void) {
    NV_Ith_S(y, 0) = -80.0;
}

static int rhs(void) {
    NV_Ith_S(ydot, 0) = NV_Ith_S(y, 0);
    return 0;
}
"""


def test_run_mapping_colocates_paired_headers(tmp_path, monkeypatch):
    workdir = tmp_path / "neutral-cwd"
    output_dir = tmp_path / "nested" / "generated"
    workdir.mkdir()
    output_dir.mkdir(parents=True)
    input_c = tmp_path / "sim.c"
    input_c.write_text(SIM_C, encoding="utf-8")
    monkeypatch.chdir(workdir)

    result = mapping_engine.run_mapping(
        input_c,
        output_dir / "Example_2024.H",
        mapping={0: "V"},
    )

    assert result["source"] == output_dir / "Example_2024.H"
    assert result["header"] == output_dir / "Example_2024Names.H"
    assert result["source"].is_file()
    assert result["header"].is_file()
    assert '#include "Example_2024Names.H"' in result["source"].read_text()
    assert not (workdir / "Example_2024Names.H").exists()


def test_cli_honours_nested_outdir_from_neutral_cwd(tmp_path):
    workdir = tmp_path / "neutral-cwd"
    output_dir = tmp_path / "nested" / "generated"
    input_dir = tmp_path / "input"
    workdir.mkdir()
    input_dir.mkdir()
    input_c = input_dir / "sim.c"
    input_c.write_text(SIM_C, encoding="utf-8")
    (workdir / "state_map.txt").write_text("0 V\n", encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            str(UTILITY_DIR / "cellML2Foam.py"),
            str(input_c),
            "--model",
            "Example_2024",
            "--outdir",
            str(output_dir),
        ],
        cwd=workdir,
        text=True,
        capture_output=True,
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    model_dir = output_dir / "Example"
    assert (model_dir / "Example_2024.H").is_file()
    assert (model_dir / "Example_2024Names.H").is_file()
    assert (model_dir / "Example.H").is_file()
    assert (model_dir / "Example.C").is_file()
    assert not (workdir / "Example_2024.H").exists()
    assert not (workdir / "Example_2024Names.H").exists()


def test_myokit_subprocess_failure_and_arguments_are_preserved(tmp_path, monkeypatch):
    cellml = tmp_path / "model.cellml"
    cellml.touch()
    failure = subprocess.CalledProcessError(9, ["myokit"])

    def fail(command, check):
        assert command == [
            "myokit",
            "import",
            "cellml",
            str(cellml),
            str(cellml.with_suffix(".mmt")),
        ]
        assert check is True
        raise failure

    monkeypatch.setattr(pipeline.subprocess, "run", fail)

    with pytest.raises(subprocess.CalledProcessError) as caught:
        pipeline.run_cellml_to_mmt(cellml)

    assert caught.value is failure
