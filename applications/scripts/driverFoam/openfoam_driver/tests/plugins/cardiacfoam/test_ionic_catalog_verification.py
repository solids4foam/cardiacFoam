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
#     test_ionic_catalog_verification
#
# Description
#     Unit tests for the ionic_catalog_verification report parser and
#     diff logic. No OpenFOAM dependency -- exercises parse_report_text /
#     diff_report against hand-written synthetic report text only. The
#     live, solver-backed check lives in
#     test_ionic_catalog_live_verification.py.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Unit tests: report parsing + catalog diffing, no OpenFOAM required.

Covers the discrepancy that motivated this whole module: constant naming is
not a static rule. `AC_`-prefixed convention for CellML-generated constants,
unprefixed for hand-added ones, and mixed within a single model (TWorld: 273
`AC_*` + one bare `gnalTissueScale`). See
`ionic_catalog_verification.py`'s module docstring.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from openfoam_driver.plugins.cardiacfoam.ionic_catalog_verification import (
    ModelVerificationResult,
    VerificationResult,
    diff_report,
    find_listCellModelsVariables_binary,
    parse_report_text,
    verify_ionic_catalog,
)
from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IonicModelEntry


def _entry(*, states=(), algebraic=(), constants=()) -> IonicModelEntry:
    """Minimal synthetic catalog entry -- deliberately independent of the
    real (large, drift-prone) IONIC_MODEL_CATALOG values so these tests stay
    stable across catalog edits."""
    return IonicModelEntry(
        states=tuple(states),
        algebraic=tuple(algebraic),
        constants=tuple(constants),
        recommended_exports=(),
        compatible_tissues=("myocyte",),
        compatible_solvers=("singleCellSolver",),
        species=("generic",),
        cardiac_region=("ventricle",),
        model_type="ionic",
        description="synthetic test fixture",
    )


def _make_report(*, constants=(), states=(), algebraic=(), ionic_model="Fixture") -> str:
    """Build report text in the exact documented listCellModelsVariables.C
    format (listCellModelsVariables.C:168-198)."""
    lines = [
        "",
        "========== listCellModelsVariables ==========",
        "",
        "Selected physicsModel: electroModel",
        "Selected myocardiumSolver: singleCellSolver",
        f"Selected ionicModel: {ionic_model}",
        "",
        f"Ionic constants ({len(constants)}) --> initial value",
    ]
    for i, name in enumerate(constants):
        lines.append(f"  constants [{i}] {name} --> 1.0")
    lines += ["", f"Ionic states ({len(states)}) --> initial value"]
    for i, name in enumerate(states):
        lines.append(f"  states [{i}] {name} --> 0.0")
    lines += ["", f"Ionic algebraic ({len(algebraic)})"]
    for i, name in enumerate(algebraic):
        lines.append(f"  algebraic [{i}] {name}")
    lines += ["", "No activeTensionModel configured in electroProperties.", ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# parse_report_text / diff_report -- pure, no I/O
# ---------------------------------------------------------------------------


def test_report_matching_catalog_is_reported_as_match():
    entry = _entry(constants=("R", "T"), states=("V",), algebraic=("Istim",))
    report = _make_report(constants=["R", "T"], states=["V"], algebraic=["Istim"])
    parsed = parse_report_text(report)

    result = diff_report("Fixture", parsed, entry)

    assert isinstance(result, ModelVerificationResult)
    assert result.status == "match"
    assert result.missing_from_catalog == {}
    assert result.extra_in_catalog == {}


def test_extra_runtime_constant_is_missing_from_catalog():
    """Runtime has a constant the catalog doesn't -- the dangerous direction:
    an agent using the catalog wouldn't know this override name exists."""
    entry = _entry(constants=("R", "T"))
    report = _make_report(constants=["R", "T", "F"])
    parsed = parse_report_text(report)

    result = diff_report("Fixture", parsed, entry)

    assert result.status == "mismatch"
    assert result.missing_from_catalog == {"constants": ("F",)}
    assert result.extra_in_catalog == {}


def test_catalog_constant_absent_at_runtime_is_extra_in_catalog():
    """Catalog claims a constant the solver no longer has -- an agent would
    write a FatalError-triggering override."""
    entry = _entry(constants=("R", "T", "StaleConstant"))
    report = _make_report(constants=["R", "T"])
    parsed = parse_report_text(report)

    result = diff_report("Fixture", parsed, entry)

    assert result.status == "mismatch"
    assert result.extra_in_catalog == {"constants": ("StaleConstant",)}
    assert result.missing_from_catalog == {}


def test_malformed_report_is_surfaced_as_mismatch_not_silently_ignored():
    entry = _entry(constants=("R", "T"), states=("V",), algebraic=("Istim",))
    garbage = "this is not a listCellModelsVariables report at all\nrandom noise\n"
    parsed = parse_report_text(garbage)

    assert parsed == {"ionic_model": None, "states": [], "algebraic": [], "constants": []}

    result = diff_report("Fixture", parsed, entry)

    assert result.status == "mismatch"
    assert result.extra_in_catalog["constants"] == ("R", "T")
    assert result.extra_in_catalog["states"] == ("V",)
    assert result.extra_in_catalog["algebraic"] == ("Istim",)
    assert result.missing_from_catalog == {}


def test_unprefixed_naming_convention_TNNP_style_matches_cleanly():
    """TNNP-style: no AC_ prefix at all (R, T, F, g_Na, ...)."""
    names = ("R", "T", "F", "g_Na", "g_K1")
    entry = _entry(constants=names)
    report = _make_report(constants=list(names))
    parsed = parse_report_text(report)

    result = diff_report("TNNP", parsed, entry)

    assert result.status == "match"


def test_mixed_naming_convention_TWorld_style_matches_cleanly():
    """TWorld-style: CellML-generated constants carry AC_, plus exactly one
    hand-added bare constant (the real case is gnalTissueScale). This is the
    finding that makes a static prefix rule impossible."""
    names = ("AC_CaMK0", "AC_K_Phos_CaMK", "AC_Whole_cell_PP1", "gnalTissueScale")
    entry = _entry(constants=names)
    report = _make_report(constants=list(names))
    parsed = parse_report_text(report)

    result = diff_report("TWorld", parsed, entry)

    assert result.status == "match"
    bare = [n for n in names if not n.startswith("AC_")]
    assert bare == ["gnalTissueScale"]


def test_diff_report_checks_all_three_fields_independently():
    entry = _entry(constants=("R",), states=("V", "K_i"), algebraic=("Istim",))
    report = _make_report(constants=["R"], states=["V"], algebraic=["Istim", "extraAlg"])
    parsed = parse_report_text(report)

    result = diff_report("Fixture", parsed, entry)

    assert result.status == "mismatch"
    assert result.extra_in_catalog == {"states": ("K_i",)}
    assert result.missing_from_catalog == {"algebraic": ("extraAlg",)}


# ---------------------------------------------------------------------------
# verify_ionic_catalog -- structured, "skipped" is never "verified"
# ---------------------------------------------------------------------------


def test_find_binary_returns_none_or_a_path():
    result = find_listCellModelsVariables_binary()
    assert result is None or isinstance(result, Path)


def test_verify_ionic_catalog_reports_skipped_not_verified_when_utility_absent(monkeypatch):
    import openfoam_driver.plugins.cardiacfoam.ionic_catalog_verification as mod

    monkeypatch.setattr(mod, "find_listCellModelsVariables_binary", lambda: None)

    result = verify_ionic_catalog("TNNP")

    assert isinstance(result, VerificationResult)
    assert result.utility_available is False
    assert result.results["TNNP"].status == "skipped"
    # A skip must never be mistaken for a pass.
    assert result.all_match is False


def test_verify_ionic_catalog_all_models_skipped_when_utility_absent(monkeypatch):
    import openfoam_driver.plugins.cardiacfoam.ionic_catalog_verification as mod
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG

    monkeypatch.setattr(mod, "find_listCellModelsVariables_binary", lambda: None)

    result = verify_ionic_catalog()

    assert set(result.results) == set(IONIC_MODEL_CATALOG)
    assert all(r.status == "skipped" for r in result.results.values())


def test_verify_ionic_catalog_rejects_unknown_model():
    with pytest.raises(KeyError):
        verify_ionic_catalog("NotARealIonicModel")


# ---------------------------------------------------------------------------
# The live path's non-solver half -- exercisable without OpenFOAM
# ---------------------------------------------------------------------------


def test_case_synthesis_works_without_the_solver(tmp_path):
    """The live check needs OpenFOAM, but its case-synthesis half does not.

    Exercising it here catches signature errors in the driver APIs the live
    path calls -- which would otherwise surface only on a user's first real
    run, long after this code was written. It already caught one:
    provision_mesh is keyword-only.
    """
    from openfoam_driver.plugins.cardiacfoam.ionic_catalog_verification import (
        _synthesize_case,
    )
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import (
        IONIC_MODEL_CATALOG,
    )

    entry = IONIC_MODEL_CATALOG["TNNP"]
    case = tmp_path / "case"
    _synthesize_case(case, "TNNP", entry)

    assert (case / "constant" / "electroProperties").is_file()
    assert (case / "constant" / "physicsProperties").is_file()
    # singleCellSolver gets the bundled 1-cell mesh copied in directly, so no
    # blockMesh step is needed before the utility runs.
    assert (case / "constant" / "polyMesh").is_dir()

    electro = (case / "constant" / "electroProperties").read_text()
    assert "TNNP" in electro
    assert "singleCellSolver" in electro


def test_synthesis_is_attempted_for_every_catalogued_model(tmp_path):
    """A model the driver cannot even configure is a catalog problem in its
    own right -- report it per model rather than letting it blind the run."""
    from openfoam_driver.plugins.cardiacfoam.ionic_catalog_verification import (
        _synthesize_case,
    )
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import (
        IONIC_MODEL_CATALOG,
    )

    failures = {}
    for name, entry in IONIC_MODEL_CATALOG.items():
        try:
            _synthesize_case(tmp_path / name, name, entry)
        except Exception as exc:  # noqa: BLE001 - collected, not swallowed
            failures[name] = str(exc)
    assert failures == {}, f"models the driver cannot configure: {failures}"
