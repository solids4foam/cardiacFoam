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
#     output_collection
#
# Description
#     Generic snapshot/diff archiving of a case's postProcessing/ output.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Archives a sweep case's raw output without knowing what produced it.

Every C++ output path in this codebase -- the manufactured verifiers
(bidomain/bath/eikonal/ECG/electromechanics), the Purkinje writer, the
single-cell solver -- independently computes the same
`<time>.globalPath()/"postProcessing"` location before writing anything
(confirmed by reading each one directly). OpenFOAM's own functionObjects
(e.g. Niederer's Niedererpoints/Niedererlines, configured in system/
controlDict) use the same convention natively. So does any custom
functionObject or verification model a user adds later -- and this module
never needs to be taught about it, because it doesn't key off filenames or
solver type at all.

Instead: snapshot what's under postProcessing/ before a case runs, diff
against what's there after, and archive whatever is new or changed into a
destination directory the caller chooses. That is "this case's output" by
construction, regardless of which solver, ionic model, or functionObject
produced it, or whether driverFOAM has ever heard of it. Callers should pass
a destination that is already unique to this case -- typically the case's
own output_dir_name folder, the same directory workflow_state.json lives in
-- so which case produced what is never ambiguous and a cross-case name
collision (two cases writing the same non-case-qualified functionObject
name) is structurally impossible rather than merely detected.

For a case whose workflow_dag actually cleans first (tet mesh_family's
"clean" -> Allclean step), postProcessing/ genuinely doesn't exist yet at
snapshot time -- snapshot_postprocessing() returns {} -- so the diff
degenerates to "collect everything found after," no special-casing needed.
It's specifically the hex workflow_dags (no clean step: just mesh->solve)
where postProcessing/ persists and accumulates across sequential cases
sharing one case_root, which is exactly why a real diff (not just "collect
everything present") is needed at all rather than always trusting an empty
starting point.
"""

from __future__ import annotations

import filecmp
import os
import shutil
from pathlib import Path


class OutputCollisionError(Exception):
    """The same case_id produced different content at the same relative path
    across two collection calls -- e.g. a retry racing a partial prior
    attempt. Raised rather than silently overwritten, because silently
    overwriting would silently lose data."""


def snapshot_postprocessing(case_root: Path) -> dict[str, tuple[float, int]]:
    """Record (mtime, size) for every file currently under case_root/postProcessing.

    Returns {} if postProcessing/ doesn't exist yet (nothing has run there).
    """
    root = Path(case_root) / "postProcessing"
    if not root.is_dir():
        return {}
    snapshot: dict[str, tuple[float, int]] = {}
    for dirpath, _dirnames, filenames in os.walk(root):
        for filename in filenames:
            path = Path(dirpath) / filename
            relpath = str(path.relative_to(root))
            stat = path.stat()
            snapshot[relpath] = (stat.st_mtime, stat.st_size)
    return snapshot


def collect_new_outputs(
    case_root: Path,
    before: dict[str, tuple[float, int]],
    destination_dir: Path,
    *,
    label: str = "",
) -> list[Path]:
    """Copy every file under postProcessing/ that is new or changed since
    `before` into destination_dir, preserving its path relative to
    postProcessing/. Returns the list of archived destination paths.

    `destination_dir` is the caller's responsibility to make unique to this
    case (e.g. that case's own output_dir_name folder) -- this function does
    not append any case identifier itself. `label` is used only in the
    collision error message below, to help a human tell which case's
    archiving call raised it.

    Raises OutputCollisionError if a destination already exists with
    genuinely different content -- e.g. a retry racing a partial prior
    attempt for the same case. Re-running to identical content is not a
    collision.
    """
    root = Path(case_root) / "postProcessing"
    destination_root = Path(destination_dir)
    if not root.is_dir():
        return []

    archived: list[Path] = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for filename in filenames:
            source = Path(dirpath) / filename
            relpath = str(source.relative_to(root))
            stat = source.stat()
            current = (stat.st_mtime, stat.st_size)
            if before.get(relpath) == current:
                continue  # unchanged leftover from before this case ran

            destination = destination_root / relpath
            if destination.exists():
                if filecmp.cmp(source, destination, shallow=False):
                    continue  # identical rerun, nothing new to archive
                raise OutputCollisionError(
                    f"{relpath!r} already exists in {destination_root} "
                    f"({label or 'no label given'}) with different content -- "
                    "likely a retry racing a partial prior attempt for this "
                    "case. Refusing to overwrite."
                )
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            archived.append(destination)
    return archived
