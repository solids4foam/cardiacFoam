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
#     function_object_fields
#
# Description
#     Plan-time validation of the field names sampled by controlDict function
#     objects. Function objects themselves are OpenFOAM's (not driverFOAM's),
#     but the *fields* they can sample are this solver's — a name the solver
#     does not expose is silently dropped at solve time. This module turns that
#     silent drop into a non-blocking warning during strict planning.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import os
from pathlib import Path
from typing import Mapping, Sequence

from foamlib import FoamFile

from openfoam_driver.core.planning_types import StrictDiagnostic, diagnostic as _diagnostic

# controlDict locations to scan: top-level and the electro sub-region used by
# multi-region electromechanical cases.
_CONTROLDICT_RELPATHS = ("system/controlDict", "system/electro/controlDict")


def _diagnostics_for_path(
    path: Path, samplable: Mapping[str, set[str]], source: str
) -> list[StrictDiagnostic]:
    try:
        functions = FoamFile(path)["functions"]
    except KeyError:
        return []
    diagnostics: list[StrictDiagnostic] = []
    for name in functions.keys():
        subdict = functions[name]
        if not hasattr(subdict, "keys"):
            # A non-dict entry inside functions -- e.g. #includeFunc, which
            # carries no braces and samples nothing of its own.
            continue
        region = subdict.get("region", "electro")
        if region not in samplable:
            # A region this plugin's samplable-fields map doesn't name
            # (e.g. a not-yet-cataloged bath/torso domain) -- skip rather
            # than guess. Forcing it into "electro" would fabricate a
            # warning for a region electro never claimed to cover; the
            # same "never surface a spurious warning from a limitation"
            # rule this module already applies to parse/IO failures.
            continue
        allowed = samplable[region]
        for field_name in subdict.get("fields", []):
            if field_name not in allowed:
                diagnostics.append(
                    _diagnostic(
                        "warning",
                        "unknown_sampled_field",
                        (
                            f"Function object samples field {field_name!r} "
                            f"(region {region!r}) which the resolved model does "
                            "not expose; the solver will silently drop it."
                        ),
                        source=source,
                        field=field_name,
                    )
                )
    return diagnostics


def function_object_field_diagnostics(
    case_root: str | Path,
    *,
    samplable: Mapping[str, Sequence[str] | set[str]],
) -> tuple[StrictDiagnostic, ...]:
    """Warn (never error) about controlDict function objects sampling fields
    absent from ``samplable`` (an open ``{region_name: {field, ...}}`` map,
    typically from :func:`capability_manifest.build_capability_manifest` --
    the built-in cardiac plugin currently declares ``"electro"`` and
    ``"solid"``, but core imposes no fixed key set). A function object whose
    ``region`` isn't a key in ``samplable`` at all is skipped rather than
    checked against a guessed bucket -- see :func:`_diagnostics_for_text`.

    Degrades to silence on any parse or IO failure — a parser limitation must
    never surface as a spurious field warning. Honors
    ``SKIP_FUNCTION_OBJECT_DIAGNOSTICS`` to bypass the check entirely.
    """
    if os.environ.get("SKIP_FUNCTION_OBJECT_DIAGNOSTICS"):
        return ()
    normalized = {region: set(fields) for region, fields in samplable.items()}
    root = Path(case_root)
    diagnostics: list[StrictDiagnostic] = []
    for relpath in _CONTROLDICT_RELPATHS:
        path = root / relpath
        if not path.is_file():
            continue
        try:
            diagnostics.extend(_diagnostics_for_path(path, normalized, str(path)))
        except Exception:
            # A parse failure must not fabricate a warning -- this module's
            # entire contract is "degrade to silence," never surface a
            # parser limitation as a false positive or a crash. Catches
            # both OSError (unreadable file) and foamlib.FoamFileDecodeError
            # (a ValueError subclass, malformed file) uniformly.
            continue
    return tuple(diagnostics)
