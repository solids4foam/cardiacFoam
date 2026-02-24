"""Compatibility shim for the upper-layer `cardiac_engine` package."""

from __future__ import annotations

import sys
from pathlib import Path

try:
    from cardiac_engine import (
        CardiacPreprocEngine,
        build_parser,
        run_pipeline_cli,
        run_pipeline_from_args,
    )
except ModuleNotFoundError:
    # Ensure old entrypoints still work without manually exporting PYTHONPATH.
    upper_src = Path(__file__).resolve().parents[3] / "engine" / "src"
    if str(upper_src) not in sys.path:
        sys.path.insert(0, str(upper_src))
    from cardiac_engine import (
        CardiacPreprocEngine,
        build_parser,
        run_pipeline_cli,
        run_pipeline_from_args,
    )

__all__ = [
    "CardiacPreprocEngine",
    "build_parser",
    "run_pipeline_cli",
    "run_pipeline_from_args",
]
