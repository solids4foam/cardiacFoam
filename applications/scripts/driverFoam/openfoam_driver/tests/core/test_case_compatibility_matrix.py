from __future__ import annotations

from pathlib import Path

import pytest

from openfoam_driver.core.runtime.registry import list_entries


def _touch(case_root: Path, relative: str) -> None:
    path = case_root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")


@pytest.mark.parametrize(
    ("files", "discovered", "runnable"),
    [
        ((), False, False),
        (("Allrun",), True, True),
        (("constant/electroProperties",), True, False),
        (("constant/electroProperties.variant",), True, False),
        (
            (
                "constant/electroProperties",
                "constant/physicsProperties",
                "system/controlDict",
                "system/fvSchemes",
                "system/fvSolution",
            ),
            True,
            True,
        ),
    ],
)
def test_existing_case_discovery_and_runnability_matrix(
    tmp_path: Path,
    files: tuple[str, ...],
    discovered: bool,
    runnable: bool,
) -> None:
    case_root = tmp_path / "candidate"
    case_root.mkdir()
    for relative in files:
        _touch(case_root, relative)

    matches = [
        entry for entry in list_entries(tmp_path)
        if entry["entry_name"] == "candidate"
    ]
    assert bool(matches) is discovered
    if discovered:
        assert matches[0]["is_runnable"] is runnable


def test_legacy_resolve_case_models_neutral_shape_has_no_cardiac_keys() -> None:
    from openfoam_driver.core.compatibility import legacy_resolve_case_models

    class NotCardiac:
        plugin_id = "org.example.notcardiac"

    result = legacy_resolve_case_models(NotCardiac(), case_root=None)
    assert result == {}


def test_legacy_samplable_fields_neutral_shape_has_no_cardiac_keys() -> None:
    from openfoam_driver.core.compatibility import legacy_samplable_fields

    class NotCardiac:
        plugin_id = "org.example.notcardiac"

    result = legacy_samplable_fields(NotCardiac(), resolved={})
    assert result == {}
