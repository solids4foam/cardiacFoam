#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_DIR))

from openfoam_driver.postprocessing.single_cell_mode_compare import main


if __name__ == "__main__":
    raise SystemExit(main())
