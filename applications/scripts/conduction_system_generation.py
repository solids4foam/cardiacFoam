#!/usr/bin/env python3
"""Top-level pipeline API for conduction system preprocessing."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS_ROOT = Path(__file__).resolve().parent
CARDIAC_CORE_ROOT = SCRIPTS_ROOT / "cardiac_core"
ENGINE_SRC = SCRIPTS_ROOT / "engine" / "src"

if str(ENGINE_SRC) not in sys.path:
    sys.path.insert(0, str(ENGINE_SRC))
if str(CARDIAC_CORE_ROOT) not in sys.path:
    sys.path.insert(0, str(CARDIAC_CORE_ROOT))

from cardiac_engine import run_pipeline_cli  # noqa: E402


def main() -> None:
    run_pipeline_cli(project_root=CARDIAC_CORE_ROOT)


if __name__ == "__main__":
    main()
