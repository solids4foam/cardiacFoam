"""Every declared distribution must actually be imported by the package.

An unused hard dependency is not cosmetic: `gmsh` pulled a ~100MB wheel to
satisfy a name that only ever appears as an external workflow command string.
"""
from __future__ import annotations

import ast
import pathlib
import re
import sys

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - the project requires >=3.11
    import tomli as tomllib

_DRIVER_ROOT = pathlib.Path(__file__).resolve().parents[3]
_PACKAGE_ROOT = _DRIVER_ROOT / "openfoam_driver"

# Distribution name -> the module name it actually installs, where they differ.
_IMPORT_NAME = {
    "pyyaml": "yaml",
    "prompt-toolkit": "prompt_toolkit",
}


def _declared() -> dict[str, set[str]]:
    data = tomllib.loads((_DRIVER_ROOT / "pyproject.toml").read_text())
    groups = {"dependencies": set(data["project"].get("dependencies", []))}
    for extra, items in data["project"].get("optional-dependencies", {}).items():
        groups[f"optional:{extra}"] = set(items)
    out: dict[str, set[str]] = {}
    for group, specs in groups.items():
        names = set()
        for spec in specs:
            dist = re.split(r"[<>=!~\[; ]", spec, maxsplit=1)[0].strip().lower()
            names.add(_IMPORT_NAME.get(dist, dist.replace("-", "_")))
        out[group] = names
    return out


def _imported_top_level_modules() -> set[str]:
    found: set[str] = set()
    for path in _PACKAGE_ROOT.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:  # pragma: no cover
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                found.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                found.add(node.module.split(".")[0])
    return found


def test_every_declared_dependency_is_imported_somewhere() -> None:
    imported = _imported_top_level_modules()
    unused = {
        group: sorted(names - imported)
        for group, names in _declared().items()
        if names - imported
    }
    assert unused == {}, (
        "Declared but never imported: "
        f"{unused}. Either import it or remove it from pyproject.toml. "
        "A name that is only an external binary (e.g. gmsh, blockMesh) belongs "
        "in CORE_NEUTRAL_COMMANDS, not in [project].dependencies."
    )
