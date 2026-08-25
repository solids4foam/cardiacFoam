"""Environment preflight gate tests.

Every test overrides SKIP_ENV_DIAGNOSTICS (set suite-wide in conftest.py)
and monkeypatches shutil.which / os.environ instead of touching the machine.
"""

import pytest

from openfoam_driver.core import strict_planning
from openfoam_driver.core.runtime.openfoam_environment import load_openfoam_environment
from openfoam_driver.core.strict_planning import (
    StrictDiagnostic,
    StrictPlanReport,
    _environment_diagnostics,
    _required_executables,
)


def _dag(*commands):
    """Build a workflow_dag dict. Each entry is a str command or (command, args)."""
    steps = []
    for entry in commands:
        if isinstance(entry, tuple):
            command, args = entry
        else:
            command, args = entry, []
        steps.append({"command": command, "args": list(args)})
    return {"steps": steps}


def _which_factory(present):
    present = set(present)

    def fake_which(name, *args, **kwargs):
        return f"/usr/bin/{name}" if name in present else None

    return fake_which


def _diags(workflow_dag):
    import os

    return _environment_diagnostics(workflow_dag, env=dict(os.environ))


def test_required_executables_collects_step_commands():
    reqs = _required_executables(_dag("blockMesh", "cardiacFoam", "setExprFields"))
    assert reqs.executables == ("blockMesh", "cardiacFoam", "setExprFields")
    assert reqs.is_parallel is False
    assert reqs.mpi_launcher_in_dag is False


def test_required_executables_skips_interpreters_and_dedupes():
    reqs = _required_executables(_dag("python3", "cardiacFoam", "cardiacFoam"))
    assert reqs.executables == ("cardiacFoam",)


def test_required_executables_unwraps_mpi_launcher():
    reqs = _required_executables(
        _dag(("mpirun", ["-np", "4", "cardiacFoam", "-parallel"]))
    )
    assert reqs.executables == ("mpirun", "cardiacFoam")
    assert reqs.is_parallel is True
    assert reqs.mpi_launcher_in_dag is True


def test_required_executables_parallel_via_flag_or_decompose():
    via_flag = _required_executables(_dag(("cardiacFoam", ["-parallel"])))
    assert via_flag.is_parallel is True
    assert via_flag.mpi_launcher_in_dag is False

    via_decompose = _required_executables(_dag("decomposePar", "cardiacFoam"))
    assert via_decompose.is_parallel is True


def test_required_executables_none_dag_is_empty():
    reqs = _required_executables(None)
    assert reqs.executables == ()
    assert reqs.is_parallel is False


def test_unwrap_mpi_program_edge_cases():
    from openfoam_driver.core.strict_planning import _unwrap_mpi_program

    assert _unwrap_mpi_program(()) is None
    assert _unwrap_mpi_program(("-np",)) is None            # value flag, no value
    assert _unwrap_mpi_program(("-np", "4")) is None        # launcher flags only, no program
    assert _unwrap_mpi_program(("--oversubscribe", "prog")) == "prog"
    assert _unwrap_mpi_program(("-np", "4", "cardiacFoam")) == "cardiacFoam"


def _make_exec(path):
    path.write_text("#!/bin/sh\n")
    path.chmod(0o755)


def test_build_staleness_warns_when_user_binary_older_than_source(tmp_path):
    import os

    from openfoam_driver.core.runtime.environment_preflight import (
        _build_staleness_diagnostics,
    )

    src_root = tmp_path / "src"
    src_root.mkdir()
    source = src_root / "solver.C"
    source.write_text("// code")
    os.utime(source, (2000, 2000))  # source newer than the binary below

    appbin = tmp_path / "appbin"
    appbin.mkdir()
    _make_exec(appbin / "myFoam")
    os.utime(appbin / "myFoam", (1000, 1000))  # stale: older than source

    env = {"PATH": str(appbin), "FOAM_USER_APPBIN": str(appbin)}
    diags = _build_staleness_diagnostics(_dag("myFoam"), env, src_root=src_root)

    assert len(diags) == 1
    assert diags[0].level == "warning"
    assert diags[0].code == "stale_build"
    assert diags[0].field == "myFoam"


def test_build_staleness_silent_when_binary_newer_than_source(tmp_path):
    import os

    from openfoam_driver.core.runtime.environment_preflight import (
        _build_staleness_diagnostics,
    )

    src_root = tmp_path / "src"
    src_root.mkdir()
    source = src_root / "solver.C"
    source.write_text("// code")
    os.utime(source, (1000, 1000))  # older than the binary

    appbin = tmp_path / "appbin"
    appbin.mkdir()
    _make_exec(appbin / "myFoam")
    os.utime(appbin / "myFoam", (2000, 2000))  # fresh

    env = {"PATH": str(appbin), "FOAM_USER_APPBIN": str(appbin)}
    diags = _build_staleness_diagnostics(_dag("myFoam"), env, src_root=src_root)
    assert diags == ()


def test_build_staleness_ignores_binaries_outside_user_appbin(tmp_path):
    # A core/system binary (not under FOAM_USER_APPBIN) is never flagged even if
    # older than source -- only user-compiled utilities are policed.
    import os

    from openfoam_driver.core.runtime.environment_preflight import (
        _build_staleness_diagnostics,
    )

    src_root = tmp_path / "src"
    src_root.mkdir()
    source = src_root / "solver.C"
    source.write_text("// code")
    os.utime(source, (2000, 2000))

    sysbin = tmp_path / "sysbin"
    sysbin.mkdir()
    _make_exec(sysbin / "blockMesh")
    os.utime(sysbin / "blockMesh", (1000, 1000))  # older, but not user-compiled

    env = {"PATH": str(sysbin), "FOAM_USER_APPBIN": str(tmp_path / "appbin")}
    diags = _build_staleness_diagnostics(_dag("blockMesh"), env, src_root=src_root)
    assert diags == ()


def test_build_staleness_no_src_root_is_silent(tmp_path):
    from openfoam_driver.core.runtime.environment_preflight import (
        _build_staleness_diagnostics,
    )

    env = {"PATH": "", "FOAM_USER_APPBIN": str(tmp_path)}
    assert _build_staleness_diagnostics(_dag("myFoam"), env, src_root=None) == ()


@pytest.fixture
def clean_env(monkeypatch):
    """A fully-sourced OpenFOAM env with the preflight gate enabled."""
    monkeypatch.delenv("SKIP_ENV_DIAGNOSTICS", raising=False)
    monkeypatch.setenv("WM_PROJECT_DIR", "/opt/openfoam")
    monkeypatch.setenv("WM_PROJECT_VERSION", "v2406")
    monkeypatch.setenv("FOAM_USER_LIBBIN", "/home/u/platforms/lib")
    return monkeypatch


def test_missing_executable_is_error(clean_env):
    clean_env.setattr(
        strict_planning.shutil, "which", _which_factory({"blockMesh", "cardiacFoam"})
    )
    diags = _diags(_dag("blockMesh", "cardiacFoam", "setExprFields"))
    missing = [d for d in diags if d.code == "missing_executable"]
    assert len(missing) == 1
    assert missing[0].level == "error"
    assert missing[0].field == "setExprFields"


def test_present_executables_have_no_error(clean_env):
    clean_env.setattr(
        strict_planning.shutil, "which", _which_factory({"blockMesh", "cardiacFoam"})
    )
    diags = _diags(_dag("blockMesh", "cardiacFoam"))
    assert "missing_executable" not in {d.code for d in diags}


def test_case_script_commands_are_not_path_checked(clean_env):
    # Allrun/Allclean/etc are case-local scripts (CASE_SCRIPT_COMMANDS in
    # workflow.py), resolved relative to caseRoot at execution time -- they
    # are never on PATH by design, so shutil.which() must never be asked
    # about them. Previously this produced a false-positive
    # "missing_executable" for every Allrun-routed plan (e.g. every
    # sweep-run case), even with OpenFOAM fully sourced.
    clean_env.setattr(strict_planning.shutil, "which", _which_factory(set()))
    diags = _diags(_dag("Allrun"))
    assert "missing_executable" not in {d.code for d in diags}


def test_mpi_wrapper_checks_launcher_and_program(clean_env):
    clean_env.setattr(strict_planning.shutil, "which", _which_factory({"mpirun"}))
    diags = _diags(
        _dag(("mpirun", ["-np", "4", "cardiacFoam", "-parallel"]))
    )
    missing = {d.field for d in diags if d.code == "missing_executable"}
    assert missing == {"cardiacFoam"}  # mpirun present, cardiacFoam missing
    assert "missing_mpi" not in {d.code for d in diags}  # launcher is in the dag


def test_parallel_without_launcher_reports_missing_mpi(clean_env):
    clean_env.setattr(strict_planning.shutil, "which", _which_factory({"cardiacFoam"}))
    diags = _diags(_dag(("cardiacFoam", ["-parallel"])))
    assert "missing_mpi" in {d.code for d in diags}


def test_serial_plan_has_no_mpi_error(clean_env):
    clean_env.setattr(strict_planning.shutil, "which", _which_factory({"cardiacFoam"}))
    diags = _diags(_dag("cardiacFoam"))
    assert "missing_mpi" not in {d.code for d in diags}


def test_missing_wm_project_dir_is_error(clean_env):
    clean_env.delenv("WM_PROJECT_DIR", raising=False)
    clean_env.setattr(strict_planning.shutil, "which", _which_factory({"cardiacFoam"}))
    diags = _diags(_dag("cardiacFoam"))
    errors = [d for d in diags if d.code == "missing_openfoam_env"]
    assert len(errors) == 1
    assert errors[0].level == "error"


def test_partial_openfoam_env_is_warning_only(clean_env):
    clean_env.delenv("WM_PROJECT_VERSION", raising=False)
    clean_env.setattr(strict_planning.shutil, "which", _which_factory({"cardiacFoam"}))
    diags = _diags(_dag("cardiacFoam"))
    partial = [d for d in diags if d.code == "partial_openfoam_env"]
    assert len(partial) == 1
    assert partial[0].level == "warning"
    assert partial[0].field == "WM_PROJECT_VERSION"


def test_skip_env_diagnostics_short_circuits(monkeypatch):
    monkeypatch.setenv("SKIP_ENV_DIAGNOSTICS", "1")
    monkeypatch.setattr(strict_planning.shutil, "which", _which_factory(set()))
    assert _diags(_dag("cardiacFoam")) == ()


def test_both_partial_env_vars_missing_yield_two_warnings(clean_env):
    clean_env.delenv("WM_PROJECT_VERSION", raising=False)
    clean_env.delenv("FOAM_USER_LIBBIN", raising=False)
    clean_env.setattr(strict_planning.shutil, "which", _which_factory({"cardiacFoam"}))
    diags = _diags(_dag("cardiacFoam"))
    partial = {d.field for d in diags if d.code == "partial_openfoam_env"}
    assert partial == {"WM_PROJECT_VERSION", "FOAM_USER_LIBBIN"}


def test_no_partial_warnings_when_openfoam_unsourced(clean_env):
    clean_env.delenv("WM_PROJECT_DIR", raising=False)
    clean_env.delenv("WM_PROJECT_VERSION", raising=False)
    clean_env.delenv("FOAM_USER_LIBBIN", raising=False)
    clean_env.setattr(strict_planning.shutil, "which", _which_factory({"cardiacFoam"}))
    diags = _diags(_dag("cardiacFoam"))
    codes = {d.code for d in diags}
    assert "missing_openfoam_env" in codes
    assert "partial_openfoam_env" not in codes


def test_load_openfoam_environment_sources_bashrc(tmp_path, monkeypatch):
    bashrc = tmp_path / "bashrc"
    foam_bin = tmp_path / "bin"
    foam_bin.mkdir()
    (foam_bin / "cardiacFoam").write_text("#!/bin/sh\n")
    (foam_bin / "cardiacFoam").chmod(0o755)
    bashrc.write_text(
        f"export WM_PROJECT_DIR={tmp_path}\n"
        "export WM_PROJECT_VERSION=v2412\n"
        f"export FOAM_USER_LIBBIN={tmp_path / 'lib'}\n"
        f"export PATH={foam_bin}:$PATH\n"
    )

    monkeypatch.delenv("SKIP_ENV_DIAGNOSTICS", raising=False)
    loaded = load_openfoam_environment(explicit_bashrc=bashrc, base_env={})

    assert loaded.error is None
    assert loaded.env["WM_PROJECT_DIR"] == str(tmp_path)
    assert loaded.env["WM_PROJECT_VERSION"] == "v2412"
    diags = _environment_diagnostics(_dag("cardiacFoam"), env=loaded.env)
    assert not [d for d in diags if d.level == "error"]


def test_report_to_json_contains_environment_diagnostics():
    report = StrictPlanReport(
        status="ok",
        entry="singleCell",
        resolved_entry={},
        environment_diagnostics=(
            StrictDiagnostic(level="error", code="missing_executable",
                             message="x", source="environment", field="cardiacFoam"),
        ),
    )
    payload = report.to_json()
    assert "environment_diagnostics" in payload
    entry = payload["environment_diagnostics"][0]
    assert entry["field"] == "cardiacFoam"
    assert entry["level"] == "error"
    assert entry["code"] == "missing_executable"


def test_report_to_json_environment_diagnostics_defaults_empty():
    report = StrictPlanReport(status="ok", entry="singleCell", resolved_entry={})
    assert report.to_json()["environment_diagnostics"] == []
