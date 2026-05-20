"""P5d — Cross-fixture regression guard for validate_run.

For each of the 7 tutorial spec fixtures, build a representative RunDocument
that reflects the spec's solver type and assert that ``validate_run`` returns
zero *error*-level violations.

This catches accidentally over-restrictive structured constraints introduced
during the P5b migration.  Warnings are permitted; only ``level="error"``
must be empty for each fixture run.

Fixture-to-solver mapping (derived from each spec's defaults.ELECTRO_PROPERTIES_SCOPE):
    single_cell          → singleCellSolver
    manufactured_fda     → monodomainSolver  (ionic = monodomainFDAManufactured)
    manufactured_fda_bidomain → bidomainSolver (ionic = bidomainFDAManufactured)
    manufactured_fda_bath_bidomain → bidomainSolver (ionic = bathBidomainFDAManufactured)
    niederer_2012        → monodomainSolver  (ionic = TNNP, tissue = epicardialCells)
    restitution_curves   → singleCellSolver  (ionic = TNNP, tissue = epicardialCells)
    generic_case         → monodomainSolver  (representative; generic_case is
                           solver-agnostic, monodomainSolver is the most common)
"""

from __future__ import annotations

import pytest

from openfoam_driver.dict_entries import (
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.specs.validation import ValidationError, slot_key, validate_run

_PHASE_ORDER = ("anatomy", "physics", "stimulus", "solver")


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        yield from group


def _filled_run_for_solver(myocardium_solver: str, **extra_config) -> RunDocument:
    """Build a RunDocument with every required entry pre-populated.

    Uses the same pattern as ``_filled_run`` in test_validation.py, then
    applies solver-specific overrides so the correct solver is selected and
    any solver-specific required entries are populated.  ``extra_config``
    maps phase → {slot_key: value} for additional overrides.
    """
    config: dict[str, dict] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    # Pre-populate all required=True entries with plausible stubs.
    for e in _all_entries():
        if not e.required:
            continue
        ph = next((p for p in _PHASE_ORDER if p in e.phases), None)
        if ph is None:
            continue
        key = slot_key(e.driver_path)
        if e.value_kind == "enum" and e.enum_values:
            config[ph][key] = e.enum_values[0]
        else:
            config[ph][key] = "stub"

    # Apply solver-specific overrides that match the fixture's actual configuration.
    config["physics"]["myocardiumSolver"] = myocardium_solver

    for ph, slice_ in extra_config.items():
        config.setdefault(ph, {}).update(slice_)

    return RunDocument(id="r1", name="r", status="draft", config=config)


# ---------------------------------------------------------------------------
# Fixtures parameterised by spec name + representative run
# ---------------------------------------------------------------------------

# Each tuple is (spec_label, RunDocument).
# The RunDocument is built to match the spec's actual solver and ionic model.

_FIXTURE_RUNS = [
    (
        "single_cell",
        _filled_run_for_solver(
            "singleCellSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
    (
        "manufactured_fda",
        _filled_run_for_solver(
            "monodomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "monodomainFDAManufactured",
                # tissue not required for manufactured models (applicable_when excludes them)
            },
        ),
    ),
    (
        "manufactured_fda_bidomain",
        _filled_run_for_solver(
            "bidomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "bidomainFDAManufactured",
            },
        ),
    ),
    (
        "manufactured_fda_bath_bidomain",
        _filled_run_for_solver(
            "bidomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "bathBidomainFDAManufactured",
            },
        ),
    ),
    (
        "niederer_2012",
        _filled_run_for_solver(
            "monodomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
    (
        "restitution_curves",
        _filled_run_for_solver(
            "singleCellSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
    (
        "generic_case",
        _filled_run_for_solver(
            "monodomainSolver",
            physics={
                "type": "electroModel",
                "ionicModel": "TNNP",
                "tissue": "epicardialCells",
            },
        ),
    ),
]


@pytest.mark.parametrize("spec_label,run", _FIXTURE_RUNS, ids=[t[0] for t in _FIXTURE_RUNS])
def test_representative_run_has_no_validator_errors(spec_label: str, run: RunDocument):
    """validate_run must return zero error-level violations for each fixture.

    Warnings are permitted.  An error-level violation indicates an over-
    restrictive structured constraint introduced by P5b migration.
    """
    errors = [e for e in validate_run(run) if e.level == "error"]
    assert errors == [], (
        f"spec='{spec_label}': expected no validator errors for representative run, "
        f"got:\n" + "\n".join(f"  [{e.phase}] {e.field}: {e.message}" for e in errors)
    )
