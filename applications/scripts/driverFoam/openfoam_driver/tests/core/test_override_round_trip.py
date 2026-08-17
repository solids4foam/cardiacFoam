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
#     test_override_round_trip
#
# Description
#     Every catalogued override must actually land in the emitted dictionary,
#     at the nesting the solver reads. Guards the silent-drop class.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""An override that validates but emits nothing is the worst defect class here.

Two instances shipped and survived review:

* ``cellZone`` was catalogued with a bare ``driver_path``, so the emitter wrote
  it at the ``electroProperties`` **root** while the solver reads it from the
  resolved ``<solver>Coeffs`` block. Setting it did nothing, and the run
  silently used the whole mesh -- bath included -- instead of the requested
  zone.
* The dynamic-path emitter expanded only two hardcoded placeholders,
  ``<name>`` and ``<electrode>``. The ionic constant overrides used
  ``<AC_name>``, so every driver-written drug or channelopathy override was
  dropped entirely while validation reported success.

Both passed validation. Both produced plausible results. Neither was caught by
any test, because the suite checked that overrides were *accepted*, never that
they were *emitted*.

This module closes that gap generically: for every catalogued entry, set it and
assert it appears in the output at the right nesting. It is deliberately a
whole-catalog sweep rather than a handful of examples -- the two bugs above
were in different entries, and picking examples is how they were missed.
"""

from __future__ import annotations

import pytest

from openfoam_driver.core.plugin_interface import default_driver_context
from openfoam_driver.plugins.cardiacfoam.dict_builder import (
    build_electro_properties,
    select_applicable_entries,
)

# Selector contexts to try. An entry only has to emit under ONE of them, since
# many are gated by applicable_when on the solver.
_CONTEXTS = (
    {"myocardiumSolver": "monodomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
    {"myocardiumSolver": "bidomainSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
    {"myocardiumSolver": "singleCellSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
    {"myocardiumSolver": "eikonalSolver", "ionicModel": "TNNP", "tissue": "epicardialCells"},
)

_COEFFS_PREFIX = "$ELECTRO_MODEL_COEFFS."

# A plausible value per value_kind. The point is emission, not validity, so
# these only need to survive the builder.
_VALUES = {
    "scalar": "1.5",
    "integer": "3",
    "label": "3",
    "boolean": "yes",
    "word": "probeValue",
    "word_list": "(alpha beta)",
    "wordList": "(alpha beta)",
    "scalar_list": "(1 2 3)",
    "scalarList": "(1 2 3)",
    "label_list": "(0 1 2)",
    "vector3": "(1 0 0)",
    "vector3_list": "((1 0 0) (0 1 0))",
    "dimensioned_scalar_literal": "[0 0 0 0 0 0 0] 1.0",
    "dimensioned_tensor_literal": "[0 0 0 0 0 0 0] (1 0 0 1 0 1)",
    "openfoam_literal": "1.0",
}


def _concrete_path(driver_path: str) -> str:
    """Substitute each <placeholder> with a distinct concrete instance name."""
    import re

    counter = {"n": 0}

    def _sub(_match):
        counter["n"] += 1
        return f"probe{counter['n']}"

    return re.sub(r"<[A-Za-z_][A-Za-z0-9_]*>", _sub, driver_path)


def _value_for(entry) -> str:
    if entry.enum_values:
        return entry.enum_values[0]
    if entry.typical_value:
        return entry.typical_value
    return _VALUES.get(entry.value_kind, "1.0")


def _electro_entries():
    catalog = default_driver_context().capabilities.dictionaries.catalog()
    return [
        entry
        for entry in catalog.entries_for("electroProperties")
        if entry.driver_path.startswith(_COEFFS_PREFIX)
    ]


def _emit_with(entry, selectors):
    """Return (text, rejected). `rejected` means the builder raised -- which is
    acceptable: a loud refusal is not a silent drop."""
    path = _concrete_path(entry.driver_path)
    try:
        return (
            build_electro_properties(
                selectors=dict(selectors), overrides={path: _value_for(entry)}
            ),
            False,
        )
    except Exception:  # noqa: BLE001 - a refusal is a valid outcome here
        return (None, True)


def _leaf(driver_path: str) -> str:
    return _concrete_path(driver_path).split(".")[-1]


@pytest.mark.parametrize("entry", _electro_entries(), ids=lambda e: e.driver_path)
def test_an_accepted_override_is_never_silently_dropped(entry):
    """The invariant, stated precisely.

    An override must either reach the dictionary, or be **loudly refused**. What
    must never happen is the third outcome: accepted without complaint, and
    emitted nowhere. That is what made ``cellZone`` and the ionic constant
    overrides inert -- the driver reported success and wrote nothing.

    Refusal is fine and common here: many entries need a coherent multi-key
    configuration (a declared conduction-network block, a matching coupler) and
    the validator correctly rejects a single key set in isolation. Gating by
    ``applicable_when`` is fine too. Only silence is a defect.
    """
    leaf = _leaf(entry.driver_path)
    accepted_but_absent = []

    for selectors in _CONTEXTS:
        # Only assert where the entry is actually applicable. An entry gated by
        # applicable_when (a batched-only key, a TWorld-only key, a
        # verifier-specific key) is CORRECTLY absent from a context that does
        # not satisfy its gate -- that is selection working, not a silent drop.
        # Consult the driver's own predicate rather than re-implementing it.
        context = dict(selectors)
        context[_concrete_path(entry.driver_path)] = _value_for(entry)
        if entry not in select_applicable_entries(context, entries=[entry]):
            continue

        text, rejected = _emit_with(entry, selectors)
        if rejected:
            continue
        if leaf in text:
            return  # emitted somewhere: invariant satisfied
        accepted_but_absent.append(selectors["myocardiumSolver"])

    if accepted_but_absent:
        pytest.fail(
            f"{entry.driver_path} was ACCEPTED without complaint under "
            f"{accepted_but_absent} and emitted nothing. An agent setting it "
            f"would be told it succeeded while the dictionary stayed unchanged "
            f"-- the defect class that made cellZone and the ionic constant "
            f"overrides inert."
        )
    # Refused under every context: loud, therefore acceptable.


@pytest.mark.parametrize("entry", _electro_entries(), ids=lambda e: e.driver_path)
def test_coeffs_entries_land_inside_the_coeffs_block(entry):
    """Nesting, not just presence.

    The ``cellZone`` assertion generalised: a ``$ELECTRO_MODEL_COEFFS.*`` entry
    must appear INSIDE the ``<solver>Coeffs`` block. Emitted at the file root it
    parses fine and is silently ignored by the solver, which is exactly why the
    original bug survived -- the key was present, just where nothing reads it.
    """
    leaf = _leaf(entry.driver_path)
    for selectors in _CONTEXTS:
        text, rejected = _emit_with(entry, selectors)
        if rejected or leaf not in text:
            continue
        coeffs_block = f"{selectors['myocardiumSolver']}Coeffs"
        assert coeffs_block in text, f"no coeffs block emitted for {entry.driver_path}"
        assert text.index(leaf) > text.index(coeffs_block), (
            f"{entry.driver_path} was emitted at the electroProperties ROOT, "
            f"before {coeffs_block}. The solver reads it from inside the coeffs "
            f"block, so it would be silently ignored."
        )
        return
