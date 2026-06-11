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
#     test_ionic_heterogeneity
#
# Description
#     Tests ionic heterogeneity logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Phase 2 — driverFOAM tissue-heterogeneity wiring.

Covers the four surfaces wired in Phase 2:
  1. ionic_model_catalog: ``supports_heterogeneity`` flag.
  2. dict_entries: the seven ``ionicHeterogeneity.*`` DictEntries (gated).
  3. dict_builder: build + parse round-trip of a heterogeneity block
     (proves the generic nested-path machinery needs no builder change).
  4. validation: model-capability gate, endo<mEpi ordering, tissue compat.
"""

from __future__ import annotations

from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.specs.validation import validate_run


# --------------------------------------------------------------------------
# 1) Catalog flag
# --------------------------------------------------------------------------

def test_supports_heterogeneity_flag_for_capable_scalar_models():
    from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in ("BuenoOrovio", "TNNP", "TWorld", "ToRORd_dynCl"):
        assert IONIC_MODEL_CATALOG[name].supports_heterogeneity is True, name


def test_supports_heterogeneity_inherited_by_batched_variants():
    from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in (
        "BuenoOroviocompactBatched", "TNNPcompactBatched",
        "TWorldcompactBatched", "ToRORd_dynClcompactBatched",
    ):
        assert IONIC_MODEL_CATALOG[name].supports_heterogeneity is True, name


def test_single_tissue_models_do_not_support_heterogeneity():
    from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in ("AlievPanfilov", "Courtemanche", "Stewart", "PerisYague"):
        assert IONIC_MODEL_CATALOG[name].supports_heterogeneity is False, name


# --------------------------------------------------------------------------
# 2) dict_entries
# --------------------------------------------------------------------------

def _het_entries():
    from openfoam_driver.dict_entries import ELECTRO_PROPERTY_ENTRY_GROUPS
    return ELECTRO_PROPERTY_ENTRY_GROUPS["ionic_heterogeneity"]


def test_all_seven_heterogeneity_entries_exist():
    paths = {e.driver_path for e in _het_entries()}
    expected = {
        f"$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.{leaf}"
        for leaf in (
            "field", "mode", "endoMInterface", "mEpiInterface",
            "transitionWidth", "transitionMode", "smoothing",
        )
    }
    assert paths == expected


def test_heterogeneity_entries_gated_to_capable_models():
    from openfoam_driver.dict_entries import HETEROGENEITY_MODELS
    for e in _het_entries():
        assert e.applicable_when.get("ionicModel") == HETEROGENEITY_MODELS, e.driver_path


def test_heterogeneity_enum_values():
    by_leaf = {e.driver_path.rsplit(".", 1)[-1]: e for e in _het_entries()}
    assert by_leaf["mode"].enum_values == ("transmuralBands",)
    assert by_leaf["transitionMode"].enum_values == ("blend", "hard")
    assert by_leaf["smoothing"].enum_values == ("smoothstep",)


def test_heterogeneity_entries_are_optional():
    # Opt-in: no typical_value, not required, so default builds omit the block.
    for e in _het_entries():
        assert e.required is False
        assert e.typical_value == ""


# --------------------------------------------------------------------------
# 3) dict_builder round-trip (no builder code change required)
# --------------------------------------------------------------------------

_HET_OVERRIDES = {
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.field": "t",
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mode": "transmuralBands",
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.endoMInterface": "0.25",
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.mEpiInterface": "0.75",
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionWidth": "0.1",
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.transitionMode": "blend",
    "$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.smoothing": "smoothstep",
}


def test_build_emits_nested_heterogeneity_block():
    from openfoam_driver.specs.dict_builder import build_electro_properties
    text = build_electro_properties(
        selectors={
            "myocardiumSolver": "monodomainSolver",
            "ionicModel": "BuenoOrovio",
            "tissue": "epicardialCells",
        },
        overrides=_HET_OVERRIDES,
    )
    assert "ionicHeterogeneity" in text
    assert "endoMInterface 0.25;" in text
    assert "transitionMode blend;" in text


def test_build_then_parse_round_trips_heterogeneity(tmp_path):
    from openfoam_driver.specs.dict_builder import (
        build_electro_properties,
        parse_electro_properties,
    )
    text = build_electro_properties(
        selectors={
            "myocardiumSolver": "monodomainSolver",
            "ionicModel": "BuenoOrovio",
            "tissue": "epicardialCells",
        },
        overrides=_HET_OVERRIDES,
    )
    path = tmp_path / "electroProperties"
    path.write_text(text)

    parsed = parse_electro_properties(path)
    overrides = parsed["overrides"]
    for key, value in _HET_OVERRIDES.items():
        assert overrides.get(key) == value, key


def test_default_build_omits_heterogeneity_block():
    # Heterogeneity must be opt-in: a capable model with no het overrides
    # produces no ionicHeterogeneity block.
    from openfoam_driver.specs.dict_builder import build_electro_properties
    text = build_electro_properties(
        selectors={
            "myocardiumSolver": "monodomainSolver",
            "ionicModel": "BuenoOrovio",
            "tissue": "epicardialCells",
        },
    )
    assert "ionicHeterogeneity" not in text


# --------------------------------------------------------------------------
# 4) Validation
# --------------------------------------------------------------------------

def _run(physics: dict) -> RunDocument:
    config = {"anatomy": {}, "physics": physics, "stimulus": {}, "solver": {}}
    return RunDocument(id="r1", name="r", status="draft", config=config)


def test_heterogeneity_with_incapable_model_is_error():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "AlievPanfilov",
        "tissue": "myocyte",
        "ionicHeterogeneity.field": "t",
        "ionicHeterogeneity.mode": "transmuralBands",
    })
    errors = [e for e in validate_run(run)
              if e.level == "error" and "heterogeneity" in e.message.lower()]
    assert len(errors) == 1, [e.message for e in validate_run(run)]


def test_heterogeneity_with_capable_model_no_het_error():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "BuenoOrovio",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.field": "t",
        "ionicHeterogeneity.mode": "transmuralBands",
        "ionicHeterogeneity.endoMInterface": "0.3",
        "ionicHeterogeneity.mEpiInterface": "0.7",
    })
    het_errors = [e for e in validate_run(run)
                  if e.level == "error" and "heterogeneity" in e.message.lower()]
    assert het_errors == []


def test_endoM_must_be_less_than_mEpi():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "BuenoOrovio",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.field": "t",
        "ionicHeterogeneity.endoMInterface": "0.8",
        "ionicHeterogeneity.mEpiInterface": "0.3",
    })
    errors = [e for e in validate_run(run)
              if e.level == "error" and "endoMInterface" in e.message]
    assert len(errors) == 1


def test_endoM_less_than_mEpi_is_silent():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "BuenoOrovio",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.field": "t",
        "ionicHeterogeneity.endoMInterface": "0.3",
        "ionicHeterogeneity.mEpiInterface": "0.7",
    })
    errors = [e for e in validate_run(run) if "endoMInterface" in e.message]
    assert errors == []


def test_tissue_incompatible_with_model_is_error():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "AlievPanfilov",   # myocyte-only
        "tissue": "epicardialCells",
    })
    errors = [e for e in validate_run(run)
              if e.level == "error" and "compatible tissues" in e.message]
    assert len(errors) == 1, [e.message for e in validate_run(run)]


def test_tissue_compatible_with_model_is_silent():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "BuenoOrovio",
        "tissue": "epicardialCells",
    })
    issues = [e for e in validate_run(run) if "compatible tissues" in e.message]
    assert issues == []
