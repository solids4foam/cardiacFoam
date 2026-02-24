"""Compatibility shim for upper-layer Purkinje slab CLI."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

def main() -> None:
    warnings.warn(
        "cardiac_preproc.cli.purkinje_slab is a Phase 1 compatibility shim. "
        "Use purkinje/slab/slab.py instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    repo_root = Path(__file__).resolve().parents[4]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    from purkinje.slab.slab import main as run  # noqa: WPS433

    run()


if __name__ == "__main__":
    main()
