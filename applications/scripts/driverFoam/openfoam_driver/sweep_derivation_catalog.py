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
#     sweep_derivation_catalog
#
# Description
#     Fixed, hardcoded registry of named sweep derivation functions.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from typing import Any, Callable

from .sweep_expansion import SweepValidationError


_CASE_ID_RE = re.compile(r"^[A-Za-z0-9_.=-]+$")


def _path_safe_case_id(value: str) -> str:
    if (
        not value
        or value in {".", ".."}
        or value.strip() != value
        or "/" in value
        or "\x00" in value
        or not _CASE_ID_RE.fullmatch(value)
    ):
        raise SweepValidationError(
            f"caseId {value!r} is not path-safe; use only letters, digits, '_', '-', '.', '='"
        )
    return value


def _case_id_template(values: dict[str, Any]) -> dict[str, Any]:
    """Join every named value into a single filesystem-safe label, in order."""
    return {"caseId": _path_safe_case_id("_".join(str(v) for v in values.values()))}


SWEEP_DERIVATION_CATALOG: dict[str, Callable[[dict[str, Any]], dict[str, Any]]] = {
    "case_id_template": _case_id_template,
}


def get_derivation(name: str) -> Callable[[dict[str, Any]], dict[str, Any]]:
    """Fixed-registry lookup. No getattr/eval/dynamic import off agent input."""
    try:
        return SWEEP_DERIVATION_CATALOG[name]
    except KeyError:
        known = ", ".join(sorted(SWEEP_DERIVATION_CATALOG)) or "(none registered)"
        raise SweepValidationError(
            f"Unknown derivation '{name}'. Known derivations: {known}"
        ) from None
