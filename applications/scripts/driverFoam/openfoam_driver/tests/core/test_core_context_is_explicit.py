"""core/ must receive an explicit DriverContext, never resolve one implicitly.

`resolve_public_driver_context(None)` returns the cardiac context. A core module
that calls it silently becomes cardiacFoam for any plugin that failed to thread a
context through -- and it does so without raising, which is why this guard is
static rather than behavioural.

The cardiac default is legitimate at the public edge (openfoam_driver/*.py and
cli.py), where "no plugin supplied" genuinely means "the built-in one". It is not
legitimate inside core/.
"""
from __future__ import annotations

import ast
import pathlib

_CORE_ROOT = pathlib.Path(__file__).resolve().parents[2] / "core"

# compatibility.py defines the function; it is allowed to mention its own name.
_EXEMPT = {_CORE_ROOT / "compatibility.py"}


def _calls_resolve_public(path: pathlib.Path) -> list[int]:
    tree = ast.parse(path.read_text(), filename=str(path))
    lines = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = getattr(func, "id", None) or getattr(func, "attr", None)
        if name == "resolve_public_driver_context":
            lines.append(node.lineno)
    return lines


def test_core_never_resolves_an_implicit_driver_context() -> None:
    offenders: dict[str, list[int]] = {}
    for path in sorted(_CORE_ROOT.rglob("*.py")):
        if path in _EXEMPT or "__pycache__" in path.parts:
            continue
        hits = _calls_resolve_public(path)
        if hits:
            offenders[str(path.relative_to(_CORE_ROOT))] = hits
    assert offenders == {}, (
        "core/ modules resolving an implicit (cardiac) DriverContext:\n"
        + "\n".join(f"  {f}: lines {ls}" for f, ls in sorted(offenders.items()))
        + "\nMake driver_context a required parameter instead."
    )
