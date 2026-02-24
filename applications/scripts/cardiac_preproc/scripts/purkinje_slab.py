#!/usr/bin/env python3
"""Compatibility wrapper for Purkinje slab CLI."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from cardiac_preproc.cli.purkinje_slab import main  # noqa: E402


if __name__ == "__main__":
    warnings.warn(
        "cardiac_preproc/scripts/purkinje_slab.py is a Phase 1 compatibility wrapper. "
        "Use purkinje/slab/slab.py instead.",
        DeprecationWarning,
        stacklevel=1,
    )
    main()
