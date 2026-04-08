from __future__ import annotations

import importlib.util
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PostprocessTask:
    module_relpath: Path
    function_name: str = "run_postprocessing"
    kwargs: dict[str, Any] = field(default_factory=dict)


def _load_python_module(module_path: Path, *, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _resolve_kwargs(kwargs: dict[str, Any], *, setup_root: Path, output_dir: Path) -> dict[str, Any]:
    resolved: dict[str, Any] = {}
    for key, value in kwargs.items():
        if value == "$OUTPUT_DIR":
            resolved[key] = str(output_dir)
            continue
        if value == "$SETUP_ROOT":
            resolved[key] = str(setup_root)
            continue
        if isinstance(value, str) and value.startswith("$OUTPUT_DIR/"):
            resolved[key] = str(output_dir / value.removeprefix("$OUTPUT_DIR/"))
            continue
        if isinstance(value, str) and value.startswith("$SETUP_ROOT/"):
            resolved[key] = str(setup_root / value.removeprefix("$SETUP_ROOT/"))
            continue
        resolved[key] = value
    return resolved


def _coerce_artifact_items(value: Any) -> list[dict[str, Any]]:
    if value is None:
        return []
    if isinstance(value, (str, Path)):
        return [{"path": str(value)}]
    if isinstance(value, dict):
        if "artifacts" in value:
            artifacts = value["artifacts"]
            if not isinstance(artifacts, list):
                raise TypeError("Postprocess result key 'artifacts' must be a list")
            return _coerce_artifact_items(artifacts)
        if "path" in value:
            return [dict(value)]
        return []
    if isinstance(value, list):
        items: list[dict[str, Any]] = []
        for item in value:
            items.extend(_coerce_artifact_items(item))
        return items
    raise TypeError(
        "Postprocess function return must be None, path-like, dict, or list of those"
    )


def _normalize_artifacts(
    raw_items: list[dict[str, Any]],
    *,
    output_dir: Path,
    module_path: Path,
    function_name: str,
    task_index: int,
) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    for raw in raw_items:
        if "path" not in raw:
            raise KeyError("Artifact entry must include a 'path' key")
        raw_path = Path(str(raw["path"]))
        absolute_path = raw_path if raw_path.is_absolute() else (output_dir / raw_path)

        try:
            path_for_manifest = str(absolute_path.relative_to(output_dir))
        except ValueError:
            path_for_manifest = str(absolute_path)

        file_format = raw.get("format")
        if file_format is None:
            file_format = absolute_path.suffix.lstrip(".") or "unknown"

        normalized.append(
            {
                "path": path_for_manifest,
                "absolute_path": str(absolute_path),
                "exists": absolute_path.exists(),
                "kind": raw.get("kind", "plot"),
                "format": file_format,
                "label": raw.get("label", absolute_path.stem),
                "task_index": task_index,
                "module": module_path.name,
                "function": function_name,
            }
        )
    return normalized


def _write_plots_manifest(
    *,
    tutorial_name: str,
    output_dir: Path,
    artifacts: list[dict[str, Any]],
) -> Path:
    manifest = {
        "schema_version": "1.1",
        "tutorial": tutorial_name,
        "output_dir": str(output_dir),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "artifact_count": len(artifacts),
        "artifacts": artifacts,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / "plots.json"
    destination.write_text(json.dumps(manifest, indent=2))
    return destination


def run_postprocess_tasks(
    setup_root: Path,
    output_dir: Path,
    *,
    tutorial_name: str,
    tasks: list[PostprocessTask],
    strict_artifacts: bool = False,
) -> None:
    if not tasks:
        print(f"No post-processing tasks configured for tutorial '{tutorial_name}'.")
        return

    print(f"Post-processing tutorial '{tutorial_name}' with {len(tasks)} task(s).")
    artifacts: list[dict[str, Any]] = []

    for index, task in enumerate(tasks, start=1):
        module_path = setup_root / task.module_relpath
        module_name = f"{tutorial_name}_post_task_{index}"
        module = _load_python_module(module_path, module_name=module_name)

        if not hasattr(module, task.function_name):
            raise AttributeError(
                f"{module_path.name} does not expose {task.function_name}()"
            )

        function = getattr(module, task.function_name)
        resolved_kwargs = _resolve_kwargs(
            task.kwargs,
            setup_root=setup_root,
            output_dir=output_dir,
        )
        call_kwargs = {
            "output_dir": str(output_dir),
            "setup_root": str(setup_root),
            **resolved_kwargs,
        }

        print(
            f"[{index}/{len(tasks)}] {module_path.name}:{task.function_name} "
            f"with args={sorted(call_kwargs.keys())}"
        )
        task_result = function(**call_kwargs)
        raw_artifacts = _coerce_artifact_items(task_result)
        artifacts.extend(
            _normalize_artifacts(
                raw_artifacts,
                output_dir=output_dir,
                module_path=module_path,
                function_name=task.function_name,
                task_index=index,
            )
        )

    manifest_path = _write_plots_manifest(
        tutorial_name=tutorial_name,
        output_dir=output_dir,
        artifacts=artifacts,
    )
    print(f"Wrote post-processing artifact manifest: {manifest_path}")

    if strict_artifacts:
        missing_paths = [item["absolute_path"] for item in artifacts if not item["exists"]]
        if missing_paths:
            missing_fmt = "\n".join(f"  - {path}" for path in missing_paths)
            raise FileNotFoundError(
                "Strict artifact validation failed; missing generated artifacts:\n"
                f"{missing_fmt}"
            )
