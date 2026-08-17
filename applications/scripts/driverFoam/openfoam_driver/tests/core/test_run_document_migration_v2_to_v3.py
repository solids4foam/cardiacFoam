"""v2 -> v3 RunDocument migration (P2.3): archived cardiac v2 documents must
migrate deterministically into the open-config v3 shape."""
from __future__ import annotations

from openfoam_driver.core.runtime.run_model import RunDocument


def _archived_v2_document() -> dict:
    return {
        "version": "2",
        "id": "plan-singleCell",
        "name": "singleCell",
        "status": "planned",
        "config": {
            "anatomy": {},
            "physics": {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP"},
            "stimulus": {},
            "solver": {},
        },
        "validation": {"status": "ok", "diagnostics": []},
    }


def test_migrate_v2_preserves_config_and_bumps_version() -> None:
    migrated = RunDocument.migrate_v2(_archived_v2_document())
    assert migrated.version == "3"
    assert migrated.config == _archived_v2_document()["config"]
    payload = migrated.to_json()
    assert payload["version"] == "3"


def test_migrate_v2_rejects_non_v2_input() -> None:
    import pytest

    with pytest.raises(ValueError, match="expects a RunDocument with version '2'"):
        RunDocument.migrate_v2({"version": "3", "id": "x", "name": "x", "status": "draft", "config": {}})
