from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.plugin_interface import generic_openfoam_context
from openfoam_driver.core.strict_planning import strict_plan


def test_plain_allrun_case_plans_without_cardiac_dictionaries(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    (case_root / "Allrun").write_text("#!/bin/sh\nexit 0\n")

    report = strict_plan(
        "plainOpenFoamCase",
        overrides={"tutorials_root": str(tmp_path)},
    )

    assert report.status == "ok"
    assert report.resolved_entry["entry_kind"] == "case_folder"
    assert report.validation_diagnostics == ()
    assert report.readiness_score["blocked_stages"] == []
    assert {artifact.artifact_id for artifact in report.expected_artifacts} == {
        "core.workflow_state",
        "core.workflow_logs",
    }
    assert {
        artifact.artifact_id: artifact.path_pattern
        for artifact in report.expected_artifacts
    } == {
        "core.workflow_state": "postProcessing/workflow_state.json",
        "core.workflow_logs": "postProcessing/workflow_logs",
    }
    assert report.workflow_dag["steps"][0]["produces"] == []


def test_plain_allrun_case_works_with_the_no_domain_context(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    (case_root / "Allrun").write_text("#!/bin/sh\nexit 0\n")

    report = strict_plan(
        "plainOpenFoamCase",
        overrides={"tutorials_root": str(tmp_path)},
        driver_context=generic_openfoam_context(),
    )

    assert report.status == "ok"
    assert report.plugin["id"] == "org.driverfoam.generic-openfoam"
    assert report.run_document is not None
    assert report.run_document.plugin == report.plugin


def test_make_spec_accepts_generic_path_addressed_overrides_not_cardiac_kwargs() -> None:
    """make_spec's public signature must not require electro/physics-named
    keyword arguments -- a non-cardiac caller should be able to pass generic,
    path-addressed dictionary overrides and file relpaths."""
    import inspect
    from openfoam_driver.core.runtime.generic_case import make_spec

    params = set(inspect.signature(make_spec).parameters)
    cardiac_params = {
        "electro_properties_relpath", "physics_properties_relpath",
        "electro_property_overrides", "physics_property_overrides",
    }
    assert not (params & cardiac_params), (
        f"make_spec still has cardiac-named parameters: {params & cardiac_params}"
    )


# --------------------------------------------------------------------------
# P2.6: generic, path-addressed dictionary overrides
# --------------------------------------------------------------------------


class _MutationSpy:
    """Records what the generic factory hands its mutation callback."""

    def __init__(self) -> None:
        self.calls: list[dict] = []

    def __call__(self, case_root, case, **kwargs) -> None:
        self.calls.append({"case_root": case_root, "case": case, **kwargs})


def _spec(tmp_path: Path, **kwargs):
    from openfoam_driver.core.runtime.generic_case import make_spec

    return make_spec(tutorials_root=tmp_path, case_dir_name="aCase", **kwargs)


def test_generic_dict_file_overrides_reach_the_mutation_callback(tmp_path: Path) -> None:
    """A plugin names its own dictionary files; core just carries the mapping
    through to the callback without knowing what the names mean."""
    spy = _MutationSpy()
    spec = _spec(
        tmp_path,
        dict_file_relpaths={"turbulence": "constant/turbulenceProperties"},
        dict_file_overrides={"turbulence": {"simulationType": "laminar"}},
        _apply_case_mutation=spy,
    )
    case = spec.build_cases()[0]
    spec.apply_case(spec.case_root, case)

    assert spy.calls[0]["dict_file_relpaths"] == {
        "turbulence": Path("constant/turbulenceProperties"),
    }
    assert spy.calls[0]["dict_file_overrides"] == {
        "turbulence": {"simulationType": "laminar"},
    }


def test_deprecated_cardiac_kwargs_still_work_as_aliases(tmp_path: Path) -> None:
    """P2.6 keeps the historical names working -- they are advertised as common
    override keys and reach make_spec verbatim from --config/--set -- but they
    now arrive at the callback under the generic mapping."""
    spy = _MutationSpy()
    spec = _spec(
        tmp_path,
        electro_properties_relpath="constant/electro/electroProperties",
        electro_property_overrides={"a.b": "1"},
        physics_property_overrides={"type": "electroModel"},
        _apply_case_mutation=spy,
    )
    case = spec.build_cases()[0]
    spec.apply_case(spec.case_root, case)

    assert spy.calls[0]["dict_file_relpaths"] == {
        "electro": Path("constant/electro/electroProperties"),
        "physics": Path("constant/physicsProperties"),
    }
    assert spy.calls[0]["dict_file_overrides"] == {
        "electro": {"a.b": "1"},
        "physics": {"type": "electroModel"},
    }


def test_per_case_entries_accept_both_generic_and_deprecated_override_keys(
    tmp_path: Path,
) -> None:
    spec = _spec(
        tmp_path,
        dict_file_overrides={"electro": {"spec.level": "0"}},
        cases=[
            {"case_id": "generic", "dict_file_overrides": {"electro": {"a": "1"}}},
            {"case_id": "deprecated", "electro_property_overrides": {"a": "2"}},
            {"case_id": "inherits"},
        ],
        _apply_case_mutation=_MutationSpy(),
    )
    by_id = {case.case_id: case.params["dict_file_overrides"] for case in spec.build_cases()}

    assert by_id["generic"] == {"electro": {"a": "1"}}
    assert by_id["deprecated"] == {"electro": {"a": "2"}}
    assert by_id["inherits"] == {"electro": {"spec.level": "0"}}


def test_genuinely_unknown_keyword_still_raises_type_error(tmp_path: Path) -> None:
    """The deprecated-alias catch-all must not swallow typos."""
    import pytest

    with pytest.raises(TypeError, match="notAParameter"):
        _spec(tmp_path, notAParameter=1)


def test_generic_detection_keys_off_the_primary_declared_dict_file(
    tmp_path: Path,
) -> None:
    """A folder is generic while the first declared dictionary file is absent,
    and stops being generic once it appears. Later entries in the mapping do
    not decide this -- a case carrying only a secondary dictionary stays
    generic, which is what a plain OpenFOAM folder with an unrelated
    ``constant/`` file relies on."""
    case_root = tmp_path / "aCase"
    (case_root / "constant").mkdir(parents=True)
    relpaths = {"primary": "constant/primaryDict", "secondary": "constant/secondaryDict"}

    assert _spec(tmp_path, dict_file_relpaths=relpaths).metadata["generic_case"] is True

    (case_root / "constant" / "secondaryDict").write_text("")
    assert _spec(tmp_path, dict_file_relpaths=relpaths).metadata["generic_case"] is True

    (case_root / "constant" / "primaryDict").write_text("")
    assert _spec(tmp_path, dict_file_relpaths=relpaths).metadata["generic_case"] is False


def test_default_mapping_stays_generic_when_only_the_secondary_dict_exists(
    tmp_path: Path,
) -> None:
    """The real scenario the primary-file heuristic exists for: two shipped
    tutorials (``electroMechanicalNiedererEtAl2011``,
    ``monodomainTotalLagrangianEM``) keep ``electroProperties`` nested under
    ``constant/electro/electroProperties`` while still carrying a top-level
    ``constant/physicsProperties``. Under the DEFAULT dict-file mapping (primary
    ``electro`` -> ``constant/electroProperties``, secondary ``physics`` ->
    ``constant/physicsProperties``, from ``core.compatibility``), the presence
    of ``physicsProperties`` alone must not flip the case out of generic --
    only the primary file's presence may do that. This is not exercised by
    ``test_generic_detection_keys_off_the_primary_declared_dict_file`` above,
    which uses synthetic ``primary``/``secondary`` names rather than
    ``make_spec``'s real default mapping."""
    case_root = tmp_path / "aCase"
    (case_root / "constant").mkdir(parents=True)
    (case_root / "constant" / "physicsProperties").write_text("")

    assert _spec(tmp_path).metadata["generic_case"] is True


def test_declaring_no_dict_files_leaves_the_case_generic(tmp_path: Path) -> None:
    spec = _spec(tmp_path, dict_file_relpaths={})
    assert spec.metadata["dict_file_relpaths"] == {}
    assert spec.metadata["generic_case"] is True


def test_default_dict_files_come_from_the_named_compatibility_seam(
    tmp_path: Path,
) -> None:
    """Core declares no dictionary files of its own; the historical defaults
    arrive through core.compatibility so the vocabulary stays at that seam."""
    from openfoam_driver.core import compatibility

    with compatibility.track_fallback_calls() as calls:
        spec = _spec(tmp_path)

    assert "legacy_generic_case_dict_file_relpaths" in calls
    assert spec.metadata["dict_file_relpaths"] == {
        "electro": "constant/electroProperties",
        "physics": "constant/physicsProperties",
    }


def test_metadata_reports_dict_file_overrides_as_one_generic_flag(
    tmp_path: Path,
) -> None:
    without = _spec(tmp_path).metadata
    with_overrides = _spec(tmp_path, dict_file_overrides={"electro": {"a": "1"}}).metadata

    assert without["has_default_dict_file_overrides"] is False
    assert with_overrides["has_default_dict_file_overrides"] is True
    for stale in (
        "electro_properties_relpath",
        "physics_properties_relpath",
        "has_default_electro_property_overrides",
        "has_default_physics_property_overrides",
    ):
        assert stale not in without


def test_make_generic_case_spec_applies_no_solver_mutation(tmp_path: Path) -> None:
    """The dedicated generic entry point must never reach into a plugin's
    mutator, even though bare make_spec still defaults to the legacy seam."""
    from openfoam_driver.core import compatibility
    from openfoam_driver.core.runtime.generic_case import make_generic_case_spec

    spec = make_generic_case_spec(tutorials_root=tmp_path, case_dir_name="aCase")
    with compatibility.track_fallback_calls() as calls:
        spec.apply_case(spec.case_root, spec.build_cases()[0])

    assert calls == []


def test_bare_make_spec_still_defaults_to_the_legacy_cardiac_mutation(
    tmp_path: Path,
) -> None:
    """The other half of the same contract: a direct make_spec caller with no
    callback keeps the historical behaviour via the named seam."""
    from openfoam_driver.core import compatibility

    case_root = tmp_path / "aCase"
    (case_root / "constant").mkdir(parents=True)
    (case_root / "constant" / "physicsProperties").write_text("type electroModel;\n")

    spec = _spec(tmp_path, dict_file_overrides={"physics": {"type": "electroMechanicalModel"}})
    with compatibility.track_fallback_calls() as calls:
        spec.apply_case(spec.case_root, spec.build_cases()[0])

    assert "legacy_generic_case_mutation" in calls
    assert "electroMechanicalModel" in (case_root / "constant" / "physicsProperties").read_text()
