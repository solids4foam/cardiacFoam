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
#     run_discovery
#
# Description
#     Discovers past execution runs stored within a directory tree.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Discover past runs under a directory tree.

`list_runs(root)` walks `root` recursively and yields one parsed
state document per `workflow_state.json` found. Malformed state files are
silently skipped so an unfinished or partially-written run does not
break agent recovery workflows.

Each yielded entry is the raw manifest dict augmented with a
``_state_path`` key carrying the absolute path to the source file --
agents use it to locate sibling sidecars: ``workflow_logs/`` (per-step
stdout/stderr logs) and, if any override was ever applied via
``--apply``, ``remediation_history.jsonl``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator


def list_runs(root: Path) -> Iterator[dict]:
    """Yield every parseable `workflow_state.json` under `root`.

    Walks recursively. Order of iteration follows ``Path.rglob`` —
    filesystem-defined and not deterministic across platforms. Callers
    that need a stable order should sort by ``_state_path`` or
    ``started_at_utc``.
    """
    root = Path(root)
    if not root.is_dir():
        return
    for state_path in root.rglob("workflow_state.json"):
        try:
            payload = json.loads(state_path.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        if not isinstance(payload, dict):
            continue
        payload["_state_path"] = str(state_path)
        yield payload
