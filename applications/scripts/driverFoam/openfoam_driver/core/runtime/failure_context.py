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
#     failure_context
#
# Description
#     Builds a self-contained, agent-readable bundle for a failed workflow
#     step: echoed identifying metadata plus bounded log tails. Pure: reads
#     files only, never mutates state, never raises on a bad/missing path.
#     Interpretation/remediation stays with the agent.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

from pathlib import Path
from typing import Any

DEFAULT_TAIL_LINES = 200
DEFAULT_TAIL_BYTES = 65536  # 64 KiB default cap on the tail read window


def _tail_file(path: str | None, *, max_lines: int, max_bytes: int) -> tuple[str, bool]:
    """Return (tail_text, truncated) for the last lines of a file.

    Bounded from the end so a diverged run that produced a huge log cannot
    exhaust memory or context. Degrades to ("", False) on any I/O problem.
    """
    if not path:
        return "", False
    file_path = Path(path)
    try:
        size = file_path.stat().st_size
        with file_path.open("rb") as handle:
            read_from_start = size <= max_bytes
            if not read_from_start:
                handle.seek(size - max_bytes)
            raw = handle.read()
    except OSError:
        return "", False

    text = raw.decode("utf-8", errors="replace")
    lines = text.splitlines()
    # After a mid-file seek the first "line" is almost certainly a fragment of a
    # longer log line; drop it so the agent never reads a partial leading line.
    if not read_from_start and lines:
        lines = lines[1:]
    truncated = not read_from_start or len(lines) > max_lines
    return "\n".join(lines[-max_lines:]), truncated


def build_failure_context(
    step_state,
    *,
    max_lines: int = DEFAULT_TAIL_LINES,
    max_bytes: int = DEFAULT_TAIL_BYTES,
) -> dict[str, Any]:
    """Return a self-contained failure bundle for a failed workflow step."""
    stdout_tail, stdout_truncated = _tail_file(
        step_state.stdout_log, max_lines=max_lines, max_bytes=max_bytes
    )
    stderr_tail, stderr_truncated = _tail_file(
        step_state.stderr_log, max_lines=max_lines, max_bytes=max_bytes
    )
    return {
        "step_id": step_state.step_id,
        "attempt": step_state.attempt,
        "exit_code": step_state.exit_code,
        "diagnostics": [dict(diagnostic) for diagnostic in step_state.diagnostics],
        "stdout_log": step_state.stdout_log,
        "stderr_log": step_state.stderr_log,
        "stdout_tail": stdout_tail,
        "stderr_tail": stderr_tail,
        "stdout_truncated": stdout_truncated,
        "stderr_truncated": stderr_truncated,
    }
