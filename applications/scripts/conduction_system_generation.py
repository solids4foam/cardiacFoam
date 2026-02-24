#!/usr/bin/env python3
"""Top-level pipeline API for conduction system preprocessing."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS_ROOT = Path(__file__).resolve().parent
CARDIAC_CORE_ROOT = SCRIPTS_ROOT / "cardiac_core"

if str(SCRIPTS_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_ROOT))

from cardiac_core.engine import run_pipeline_cli  # noqa: E402


def main() -> None:
    run_pipeline_cli(project_root=CARDIAC_CORE_ROOT)


if __name__ == "__main__":
    main()
