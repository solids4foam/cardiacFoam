from __future__ import annotations

import shlex
import sys
from pathlib import Path
from typing import Any

from .core.runtime.registry import resolve_tutorial


VALID_DRIVER_ACTIONS = ("sim", "post", "all")


def describe_launch(
    action: str,
    tutorial: str,
    *,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
    tutorials_root: str | Path | None = None,
    dry_run: bool = False,
    continue_on_error: bool = False,
    python_executable: str | None = None,
) -> dict[str, Any]:
    if action not in VALID_DRIVER_ACTIONS:
        raise ValueError(f"Unsupported driver action: {action}")
    if action == "post" and dry_run:
        raise ValueError("dry_run is not valid for action='post'")

    effective_overrides = dict(overrides or {})
    if tutorials_root is not None:
        effective_overrides["tutorials_root"] = str(tutorials_root)

    resolution = resolve_tutorial(tutorial, overrides=effective_overrides)
    spec = resolution["factory"](**resolution["factory_overrides"])
    interpreter = python_executable or sys.executable

    command = [interpreter, "-m", "openfoam_driver", action, "--tutorial", tutorial]
    if dry_run:
        command.append("--dry-run")
    if continue_on_error:
        command.append("--continue-on-error")
    if config_path is not None:
        command.extend(["--config", str(config_path)])
    if tutorials_root is not None:
        command.extend(["--tutorials-root", str(tutorials_root)])

    manifest_root = spec.setup_root if dry_run else spec.output_dir
    plots_manifest_path = None
    if action in {"post", "all"} and spec.postprocess is not None:
        plots_manifest_path = str(spec.output_dir / "plots.json")

    return {
        "action": action,
        "tutorial": tutorial,
        "resolved_name": resolution["resolved_name"],
        "resolution": resolution["resolution"],
        "command": command,
        "command_display": shlex.join(command),
        "manifest_path": str(manifest_root / "run_manifest.json"),
        "plots_manifest_path": plots_manifest_path,
        "case_root": str(spec.case_root),
        "setup_root": str(spec.setup_root),
        "output_dir": str(spec.output_dir),
        "dry_run": dry_run,
        "continue_on_error": continue_on_error,
        "postprocess_available": spec.postprocess is not None,
    }


def describe_launch_matrix(
    tutorial: str,
    *,
    overrides: dict[str, Any] | None = None,
    config_path: str | Path | None = None,
    tutorials_root: str | Path | None = None,
    python_executable: str | None = None,
) -> dict[str, dict[str, Any]]:
    return {
        action: describe_launch(
            action,
            tutorial,
            overrides=overrides,
            config_path=config_path,
            tutorials_root=tutorials_root,
            python_executable=python_executable,
        )
        for action in VALID_DRIVER_ACTIONS
    }
