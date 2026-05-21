"""Pytest session bootstrap.

The driverFOAM postprocessing chain imports `matplotlib.pyplot` at module
load time (`tutorials/.../post_processing_manufactured.py`). On macOS
without a display server — and in any headless CI/SSH/Docker environment
— matplotlib's default backend (`MacOSX` or `TkAgg`) attempts a Cocoa /
display connection and crashes the interpreter, aborting the test run.

This conftest pins the backend to `Agg` (non-interactive, no display
required) BEFORE pytest collects any tests, so any subsequent import of
matplotlib sees the headless backend. It is a test-environment-only
guard; production scripts should set their own backend if they need one.
"""
from __future__ import annotations

import os

# Set the backend before matplotlib is imported anywhere. The env var
# wins over rcParams and is honoured by every later `matplotlib.pyplot`
# import in the test process.
os.environ.setdefault("MPLBACKEND", "Agg")
