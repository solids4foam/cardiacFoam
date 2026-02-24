#!/usr/bin/env python3
"""Compatibility wrapper for Purkinje slab product CLI."""

from __future__ import annotations

import runpy
from pathlib import Path


def main() -> None:
    runpy.run_path(str(Path(__file__).resolve().with_name("slab.py")), run_name="__main__")


if __name__ == "__main__":
    main()
