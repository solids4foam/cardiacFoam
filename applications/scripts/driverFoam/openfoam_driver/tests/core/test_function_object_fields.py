"""Tests for plan-time sampled-field validation of controlDict function objects."""

from pathlib import Path

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.specs.function_object_fields import (
    function_object_field_diagnostics,
)

_HEADER = (
    "FoamFile{ version 2.0; format ascii; class dictionary; object controlDict; }\n"
)


def _write_controldict(
    tmp_path: Path, functions_body: str, *, subdir: str = "system"
) -> Path:
    system = tmp_path / subdir
    system.mkdir(parents=True, exist_ok=True)
    (system / "controlDict").write_text(
        _HEADER + "functions\n{\n" + functions_body + "\n}\n"
    )
    return tmp_path


def test_known_field_no_diagnostic(tmp_path):
    root = _write_controldict(tmp_path, "probe1{ type probes; fields (Vm); }")
    diags = function_object_field_diagnostics(
        root, samplable={"electro": {"Vm"}, "solid": set()}
    )
    assert [d for d in diags if d.level == "error"] == []
    assert [d for d in diags if d.code == "unknown_sampled_field"] == []


def test_unknown_field_warns(tmp_path):
    root = _write_controldict(tmp_path, "probe1{ type probes; fields (bananas); }")
    diags = function_object_field_diagnostics(
        root, samplable={"electro": {"Vm"}, "solid": set()}
    )
    warns = [d for d in diags if d.code == "unknown_sampled_field"]
    assert len(warns) == 1
    assert warns[0].level == "warning"
    assert "bananas" in warns[0].message


def test_multiple_fields_mixed(tmp_path):
    root = _write_controldict(
        tmp_path, "probe1{ type probes; fields (Vm Iion bananas); }"
    )
    diags = function_object_field_diagnostics(
        root, samplable={"electro": {"Vm", "Iion"}, "solid": set()}
    )
    warns = [d for d in diags if d.code == "unknown_sampled_field"]
    assert [w.field for w in warns] == ["bananas"]


def test_solid_region_checked_against_solid_set(tmp_path):
    root = _write_controldict(
        tmp_path, "ta{ type volFieldValue; region solid; fields (Ta); }"
    )
    diags = function_object_field_diagnostics(
        root, samplable={"electro": {"Vm"}, "solid": {"Ta"}}
    )
    assert [d for d in diags if d.code == "unknown_sampled_field"] == []


def test_solid_field_in_electro_region_warns(tmp_path):
    root = _write_controldict(
        tmp_path, "ta{ type volFieldValue; region electro; fields (Ta); }"
    )
    diags = function_object_field_diagnostics(
        root, samplable={"electro": {"Vm"}, "solid": {"Ta"}}
    )
    warns = [d for d in diags if d.code == "unknown_sampled_field"]
    assert len(warns) == 1 and warns[0].field == "Ta"


def test_includefunc_not_flagged(tmp_path):
    root = _write_controldict(tmp_path, "#includeFunc probes")
    diags = function_object_field_diagnostics(
        root, samplable={"electro": {"Vm"}, "solid": set()}
    )
    assert diags == ()


def test_missing_controldict_is_silent(tmp_path):
    diags = function_object_field_diagnostics(
        tmp_path, samplable={"electro": set(), "solid": set()}
    )
    assert diags == ()


def test_strict_plan_exposes_field_family_and_stays_ok(monkeypatch):
    monkeypatch.setenv("SKIP_ENV_DIAGNOSTICS", "1")
    from openfoam_driver.strict_planning import strict_plan

    report = strict_plan("singleCell").to_json()
    assert "function_object_diagnostics" in report
    # A clean registered entry must not be pushed to failed by this family.
    assert report["status"] == "ok"


def test_electro_controldict_subdir_scanned(tmp_path):
    root = _write_controldict(
        tmp_path, "p{ type probes; fields (bananas); }", subdir="system/electro"
    )
    diags = function_object_field_diagnostics(
        root, samplable={"electro": {"Vm"}, "solid": set()}
    )
    warns = [d for d in diags if d.code == "unknown_sampled_field"]
    assert len(warns) == 1 and warns[0].field == "bananas"
