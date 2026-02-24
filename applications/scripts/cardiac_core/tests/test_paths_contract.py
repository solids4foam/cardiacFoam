"""Contract tests for centralized repository paths."""

from pathlib import Path

from cardiac_core.engine.paths import resolve_paths


SCRIPTS_ROOT = Path(__file__).resolve().parents[2]


def test_resolve_paths_targets_exist() -> None:
    paths = resolve_paths(project_root=SCRIPTS_ROOT / "cardiac_core")

    assert paths.diffusivity_config.exists()
    assert paths.scar_config.exists()
    assert paths.purkinje_slab_config.exists()
    assert paths.purkinje_fractal_config.exists()
    assert paths.purkinje_fractal_runner.exists()
    assert paths.convert_runner.exists()


def test_default_outputs_root_under_files_organize() -> None:
    paths = resolve_paths(project_root=SCRIPTS_ROOT / "cardiac_core")

    assert paths.default_outputs_root.is_relative_to(paths.files_organize)
