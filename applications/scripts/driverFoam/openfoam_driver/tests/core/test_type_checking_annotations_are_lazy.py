"""A module with TYPE_CHECKING imports must defer its annotations.

This guards a bug that shipped and was invisible locally for a long time.

``core/plugin_interface.py`` annotated ``get_dict_entries(self) -> tuple[DictEntry, ...]``
while importing ``DictEntry`` only under ``if TYPE_CHECKING:``. Without
``from __future__ import annotations`` that name is resolved when the class body
executes, so merely importing the module raised::

    NameError: name 'DictEntry' is not defined

on every Python before 3.14 -- that is, on 3.11, 3.12 and 3.13, three of the
four versions ``pyproject.toml`` declares support for. It took out 45 test
modules at collection time.

It was invisible to anyone developing on 3.14, where PEP 649 defers annotation
evaluation: the same tree showed a fully green suite on 3.14 and could not
collect a single test on 3.12. A version-matrix CI would catch it; so does this,
without needing one.

The rule is mechanical: if a module has an ``if TYPE_CHECKING:`` block, its
annotations may name types that do not exist at runtime, so it must defer them.
"""
from __future__ import annotations

import ast
import pathlib

_PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[2]


def _has_type_checking_block(tree: ast.Module) -> bool:
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        test = node.test
        if isinstance(test, ast.Name) and test.id == "TYPE_CHECKING":
            return True
        if isinstance(test, ast.Attribute) and test.attr == "TYPE_CHECKING":
            return True
    return False


def _defers_annotations(tree: ast.Module) -> bool:
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            if any(alias.name == "annotations" for alias in node.names):
                return True
    return False


def test_modules_with_type_checking_imports_defer_annotations() -> None:
    offenders: list[str] = []
    for path in sorted(_PACKAGE_ROOT.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        try:
            tree = ast.parse(path.read_text(), filename=str(path))
        except SyntaxError:  # pragma: no cover - not this test's business
            continue
        if _has_type_checking_block(tree) and not _defers_annotations(tree):
            offenders.append(str(path.relative_to(_PACKAGE_ROOT)))

    assert offenders == [], (
        "these modules import names under `if TYPE_CHECKING:` but do not defer "
        "annotation evaluation, so they raise NameError on import for any Python "
        "before 3.14:\n  " + "\n  ".join(offenders)
        + "\nAdd `from __future__ import annotations` at the top of each."
    )
