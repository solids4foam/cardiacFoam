#!/usr/bin/env python3
"""Compatibility wrapper for Purkinje fractal runner."""

import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_ROOT = ROOT.parent

if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))


def main() -> None:
    from purkinje.fractal_3d.purkinje3D_fractal import main as run  # noqa: WPS433
    run()


if __name__ == "__main__":
    warnings.warn(
        "cardiac_preproc/scripts/purkinje_fractal.py is a Phase 1 compatibility wrapper. "
        "Use purkinje/fractal_3d/purkinje3D_fractal.py instead.",
        DeprecationWarning,
        stacklevel=1,
    )
    main()
