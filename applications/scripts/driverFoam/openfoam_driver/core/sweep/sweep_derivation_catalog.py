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

from openfoam_driver.core.sweep.sweep_expansion import SweepValidationError


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


def _label_for(value: Any) -> str:
    """Render one axis value as a label fragment.

    List/tuple values (e.g. niederer_2012.py's dx_values=[0.5]) are flattened
    element-wise rather than stringified as a Python literal (str([0.5]) ==
    "[0.5]", which is not path-safe).
    """
    if isinstance(value, (list, tuple)):
        return "-".join(str(v) for v in value)
    return str(value)


def _join_values_path_safe(values: dict[str, Any]) -> str:
    return _path_safe_case_id("_".join(_label_for(v) for v in values.values()))


def _case_id_template(values: dict[str, Any]) -> dict[str, Any]:
    """Join every named value into a single filesystem-safe label, in order."""
    return {"caseId": _join_values_path_safe(values)}


def _output_dir_name_template(values: dict[str, Any]) -> dict[str, Any]:
    """Join every named value into a single filesystem-safe output_dir_name.

    For entry-based sweeps (an existing registered tutorial's own make_spec),
    this plays the same role case_id_template plays for generic case_folder
    sweeps: naming this case's on-disk output directory.
    """
    return {"output_dir_name": _join_values_path_safe(values)}


SWEEP_DERIVATION_CATALOG: dict[str, Callable[[dict[str, Any]], dict[str, Any]]] = {
    "case_id_template": _case_id_template,
    "output_dir_name_template": _output_dir_name_template,
}


def get_derivation(name: str) -> Callable[[dict[str, Any]], dict[str, Any]]:
    """Fixed-registry lookup. No getattr/eval/dynamic import off agent input."""
    if not isinstance(name, str):
        raise SweepValidationError(
            f"'derive' must be a string naming a registered derivation, got {name!r}"
        )
    try:
        return SWEEP_DERIVATION_CATALOG[name]
    except KeyError:
        known = ", ".join(sorted(SWEEP_DERIVATION_CATALOG)) or "(none registered)"
        raise SweepValidationError(
            f"Unknown derivation '{name}'. Known derivations: {known}"
        ) from None
