#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     conftest
#
# Description
#     Configures shared pytest fixtures and runtime dependencies.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

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
