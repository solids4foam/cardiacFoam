from __future__ import annotations

import ast
from pathlib import Path


def test_core_imports_cardiac_implementation_only_at_compatibility_boundary() -> None:
    core_root = Path(__file__).resolve().parents[2] / "core"
    offenders: list[str] = []
    for path in core_root.rglob("*.py"):
        if path.name == "compatibility.py":
            continue
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            module = ""
            if isinstance(node, ast.ImportFrom):
                module = node.module or ""
            elif isinstance(node, ast.Import):
                for alias in node.names:
                    if "plugins.cardiacfoam" in alias.name:
                        offenders.append(f"{path.relative_to(core_root)}:{node.lineno}")
                continue
            if "plugins.cardiacfoam" in module:
                offenders.append(f"{path.relative_to(core_root)}:{node.lineno}")

    assert offenders == []


def test_production_consumers_do_not_bypass_capability_bundle() -> None:
    package_root = Path(__file__).resolve().parents[2]
    offenders: list[str] = []
    for path in package_root.rglob("*.py"):
        if "tests" in path.parts or path.name in {"plugin_interface.py", "plugin_capabilities.py"}:
            continue
        text = path.read_text()
        if "driver_context.plugin." in text:
            offenders.append(str(path.relative_to(package_root)))
    assert offenders == []
