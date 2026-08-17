"""P2.5: modules outside plugins/cardiacfoam/ and core/compatibility.py must
not import cardiac plugin internals at module scope."""
from __future__ import annotations

import ast
from pathlib import Path

_PACKAGE_ROOT = Path(__file__).parents[2]
_EXEMPT = {
    _PACKAGE_ROOT / "core" / "compatibility.py",
}


def _module_level_imports_cardiac_plugin(path: Path) -> bool:
    tree = ast.parse(path.read_text())
    for node in tree.body:  # module-level only, not inside functions
        if isinstance(node, ast.ImportFrom) and node.module and "plugins.cardiacfoam" in node.module:
            return True
    return False


def test_dict_entries_has_no_module_level_cardiac_import() -> None:
    path = _PACKAGE_ROOT / "dict_entries.py"
    assert not _module_level_imports_cardiac_plugin(path), (
        f"{path} imports plugins.cardiacfoam at module scope"
    )
