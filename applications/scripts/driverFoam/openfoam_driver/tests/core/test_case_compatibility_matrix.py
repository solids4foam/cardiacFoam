from __future__ import annotations

from pathlib import Path

import pytest

from openfoam_driver.core.runtime.registry import list_entries


def _touch(case_root: Path, relative: str) -> None:
    path = case_root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")


@pytest.mark.parametrize(
    ("files", "contract", "discovered", "runnable"),
    [
        ((), None, False, False),
        (("Allrun",), None, True, True),
        (("constant/electroProperties",), None, True, False),
        (("constant/electroProperties.variant",), None, True, False),
        (
            (
                "constant/electroProperties",
                "constant/physicsProperties",
                "system/controlDict",
                "system/fvSchemes",
                "system/fvSolution",
            ),
            None,
            True,
            True,
        ),
        ((), {"steps": [{"id": "run", "command": "Allrun"}]}, True, True),
        (
            ("Allrun",),
            {"status": {"runnable_without_substitution": False}},
            True,
            False,
        ),
    ],
)
def test_existing_case_discovery_and_runnability_matrix(
    tmp_path: Path,
    files: tuple[str, ...],
    contract: dict | None,
    discovered: bool,
    runnable: bool,
) -> None:
    case_root = tmp_path / "candidate"
    case_root.mkdir()
    for relative in files:
        _touch(case_root, relative)
    if contract is not None:
        import json

        (case_root / "workflow_contract.json").write_text(json.dumps(contract))

    matches = [
        entry for entry in list_entries(tmp_path)
        if entry["entry_name"] == "candidate"
    ]
    assert bool(matches) is discovered
    if discovered:
        assert matches[0]["is_runnable"] is runnable
