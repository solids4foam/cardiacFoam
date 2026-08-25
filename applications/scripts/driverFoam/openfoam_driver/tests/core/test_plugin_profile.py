from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin
from openfoam_driver.core.plugin_profile import PluginProfile, load_plugin_profile
from openfoam_driver.plugins.cardiacfoam import runtime_profile
from openfoam_driver.plugins.cardiacfoam.runtime_profile import configure_runtime_environment
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin


def test_cardiac_profile_declares_case_files_and_cxx_provenance() -> None:
    profile = CardiacFoamPlugin().get_profile()

    assert profile.plugin_id == "org.cardiacfoam"
    assert {rule.path for rule in profile.case_files} >= {
        "system/controlDict",
        "constant/physicsProperties",
        "constant/electroProperties",
    }
    assert profile.cxx_mapping is not None
    assert profile.cxx_mapping.allowlist_path.is_file()
    assert all(path.is_dir() for path in profile.cxx_mapping.source_roots)
    assert profile.digest.startswith("sha256:")


def test_cardiac_catalog_partitions_entries_by_document() -> None:
    catalog = CardiacFoamPlugin().get_dictionary_catalog()

    assert {"electroProperties", "physicsProperties", "controlDict"} <= set(catalog.documents)
    assert {entry.driver_path for entry in catalog.entries_for("physicsProperties")} == {"type"}
    assert {entry.driver_path for entry in catalog.entries_for("controlDict")} >= {"deltaT", "endTime"}


def test_cardiac_runtime_requires_explicit_solids4foam_root(tmp_path: Path) -> None:
    del tmp_path
    env, error = CardiacFoamPlugin().configure_execution_environment({})

    assert env == {}
    assert error is not None
    assert "DRIVERFOAM_CARDIACFOAM_BACKEND" in error


def test_cardiac_runtime_exports_one_validated_solids4foam_root(tmp_path: Path) -> None:
    root = tmp_path / "solids4foam"
    header = root / "src/solids4FoamModels/physicsModel/physicsModel.H"
    ln_include = root / "src/solids4FoamModels/lnInclude/physicsModel.H"
    header.parent.mkdir(parents=True)
    ln_include.parent.mkdir(parents=True)
    header.write_text("// source header\n")
    ln_include.write_text("// generated include\n")
    manifest = tmp_path / "cardiacFoam.build.json"
    manifest.write_text(
        json.dumps({
            "backend": "full",
            "openfoam": {"root": str(tmp_path)},
            "solids4foam": {"root": str(root)},
            "linked_libraries": [
                "libsolids4FoamModels.dylib",
                "libelectroMechanicalModels.dylib",
            ],
            "artifacts": [],
        })
    )

    env, error = CardiacFoamPlugin().configure_execution_environment({
        "DRIVERFOAM_CARDIACFOAM_BACKEND": "full",
        "DRIVERFOAM_CARDIACFOAM_SOLIDS4FOAM_ROOT": str(root),
        "DRIVERFOAM_CARDIACFOAM_BUILD_MANIFEST": str(manifest),
        "WM_PROJECT_DIR": str(tmp_path),
    })

    assert error is None
    assert env["SOLIDS4FOAM_INST_DIR"] == str(root.resolve())
    assert env["DRIVERFOAM_CARDIACFOAM_SOLIDS4FOAM_ROOT"] == str(root.resolve())


def test_infer_backend_from_linked_libraries() -> None:
    options = {
        "lightweight": {
            "required_libraries": ["libphysicsModel"],
            "forbidden_libraries": ["libsolids4FoamModels", "libelectroMechanicalModels"],
        },
        "full": {
            "required_libraries": ["libsolids4FoamModels", "libelectroMechanicalModels"],
            "forbidden_libraries": ["libphysicsModel"],
        },
    }
    lightweight_linked = ("libphysicsModel.dylib", "libelectroModels.dylib")
    full_linked = ("libsolids4FoamModels.dylib", "libelectroMechanicalModels.dylib")

    assert runtime_profile._infer_backend(lightweight_linked, options) == "lightweight"
    assert runtime_profile._infer_backend(full_linked, options) == "full"
    assert runtime_profile._infer_backend(lightweight_linked + full_linked, options) is None
    assert runtime_profile._infer_backend((), options) is None


def _write_fake_library(directory: Path, bare_name: str, content: bytes) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"lib{bare_name}.dylib"
    path.write_bytes(content)
    return path


def _write_fake_lightweight_build(appbin: Path, libbin: Path, *, solver_content: bytes) -> Path:
    appbin.mkdir(parents=True, exist_ok=True)
    solver = appbin / "cardiacFoam"
    solver.write_bytes(solver_content)
    for bare_name in ("electroModels", "ionicModels", "genericWriter", "activeTensionModels", "physicsModel"):
        _write_fake_library(libbin, bare_name, f"fake-{bare_name}".encode())
    return solver


def test_build_manifest_self_generates_from_compiled_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    appbin = tmp_path / "appbin"
    libbin = tmp_path / "libbin"
    _write_fake_lightweight_build(appbin, libbin, solver_content=b"fake-solver")
    monkeypatch.setattr(runtime_profile, "_linked_library_names", lambda binary: ("libphysicsModel.dylib",))

    manifest = tmp_path / "cardiacFoam.build.json"
    assert not manifest.exists()

    env, error = configure_runtime_environment({
        "DRIVERFOAM_CARDIACFOAM_BACKEND": "lightweight",
        "DRIVERFOAM_CARDIACFOAM_BUILD_MANIFEST": str(manifest),
        "WM_PROJECT_DIR": str(tmp_path),
        "FOAM_USER_APPBIN": str(appbin),
        "FOAM_USER_LIBBIN": str(libbin),
    })

    assert error is None, error
    payload = json.loads(manifest.read_text())
    assert payload["backend"] == "lightweight"
    assert {artifact["name"] for artifact in payload["artifacts"]} == {
        "cardiacFoam",
        "libelectroModels",
        "libionicModels",
        "libgenericWriter",
        "libactiveTensionModels",
        "libphysicsModel",
    }
    solver_artifact = next(a for a in payload["artifacts"] if a["name"] == "cardiacFoam")
    assert solver_artifact["sha256"] == hashlib.sha256(b"fake-solver").hexdigest()


def test_build_manifest_self_heals_when_stale(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    appbin = tmp_path / "appbin"
    libbin = tmp_path / "libbin"
    solver = _write_fake_lightweight_build(appbin, libbin, solver_content=b"fake-solver-v2")

    manifest = tmp_path / "cardiacFoam.build.json"
    manifest.write_text(json.dumps({
        "backend": "lightweight",
        "openfoam": {"root": str(tmp_path)},
        "solids4foam": {"root": None},
        "linked_libraries": ["libphysicsModel.dylib"],
        "artifacts": [],
    }))
    stale_time = solver.stat().st_mtime - 10
    os.utime(manifest, (stale_time, stale_time))
    monkeypatch.setattr(runtime_profile, "_linked_library_names", lambda binary: ("libphysicsModel.dylib",))

    env, error = configure_runtime_environment({
        "DRIVERFOAM_CARDIACFOAM_BACKEND": "lightweight",
        "DRIVERFOAM_CARDIACFOAM_BUILD_MANIFEST": str(manifest),
        "WM_PROJECT_DIR": str(tmp_path),
        "FOAM_USER_APPBIN": str(appbin),
        "FOAM_USER_LIBBIN": str(libbin),
    })

    assert error is None, error
    payload = json.loads(manifest.read_text())
    solver_artifact = next(a for a in payload["artifacts"] if a["name"] == "cardiacFoam")
    assert solver_artifact["sha256"] == hashlib.sha256(b"fake-solver-v2").hexdigest()


def test_cardiac_runtime_file_selects_backend_and_bashrc(tmp_path: Path) -> None:
    root = tmp_path / "solids4foam"
    for relative in (
        "src/solids4FoamModels/physicsModel/physicsModel.H",
        "src/solids4FoamModels/lnInclude/physicsModel.H",
    ):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("// header\n")
    manifest = tmp_path / "cardiacFoam.build.json"
    manifest.write_text(json.dumps({
        "backend": "full",
        "openfoam": {"root": str(tmp_path)},
        "solids4foam": {"root": str(root)},
        "linked_libraries": ["libsolids4FoamModels.dylib", "libelectroMechanicalModels.dylib"],
        "artifacts": [],
    }))
    config = tmp_path / "driverfoam-runtime.yaml"
    config.write_text(
        "openfoam:\n  bashrc: /tmp/openfoam/etc/bashrc\n"
        "plugins:\n  org.cardiacfoam:\n"
        f"    backend: full\n    solids4foam_root: {root}\n"
        f"    build_manifest: {manifest}\n"
    )

    env, error = CardiacFoamPlugin().configure_execution_environment({
        "DRIVERFOAM_RUNTIME_CONFIG": str(config),
        "WM_PROJECT_DIR": str(tmp_path),
    })

    assert error is None
    assert env["DRIVERFOAM_CARDIACFOAM_BACKEND"] == "full"
    assert env["SOLIDS4FOAM_INST_DIR"] == str(root.resolve())
    assert CardiacFoamPlugin().get_openfoam_bashrc({
        "DRIVERFOAM_RUNTIME_CONFIG": str(config),
    }) == "/tmp/openfoam/etc/bashrc"


def test_generic_profile_declares_no_solver_specific_files() -> None:
    profile = GenericOpenFOAMPlugin().get_profile()

    assert profile.plugin_id == "org.driverfoam.generic-openfoam"
    assert profile.case_files == ()
    assert profile.cxx_mapping is None


def test_profile_rejects_case_path_escape(tmp_path: Path) -> None:
    path = tmp_path / "plugin.yaml"
    path.write_text(
        """schema_version: 1
plugin: {id: example.bad, api_version: '1'}
case_profile:
  dictionaries:
    - path: ../outside
      kind: openfoam_dictionary
      role: plugin.configuration
      required: always
"""
    )

    with pytest.raises(ValueError, match="escapes the case"):
        load_plugin_profile(path)


def test_profile_digest_is_stable_after_payload_mutation() -> None:
    payload = {
        "schema_version": 1,
        "plugin": {"id": "example.profile", "api_version": "1"},
        "case_profile": {"dictionaries": []},
    }
    profile = PluginProfile(
        path=Path("example-plugin.yaml"),
        plugin_id="example.profile",
        api_version="1",
        case_files=(),
        cxx_mapping=None,
        payload=payload,
    )

    digest = profile.digest
    payload["plugin"]["id"] = "example.mutated"

    assert profile.digest == digest
