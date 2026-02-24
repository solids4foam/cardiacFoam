#!/usr/bin/env python3
"""Pipeline runner for conduction system preprocessing.

Inputs:
- Steps are provided via CLI (--steps).
- Each solver reads its own default input/output paths from its config.

Outputs:
- Writes outputs defined in each solver config (defaults under outputs/).

Main:
- Delegates to cardiac_preproc.conduction_system_generation.main().
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
sys.path.insert(0, str(SRC))

from cardiac_preproc.conduction_system_generation import main  # noqa: E402


if __name__ == "__main__":
    main()
