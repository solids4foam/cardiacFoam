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
  1. ionic_model_catalog: ``supports_heterogeneity`` and
     ``supports_apex_base_heterogeneity`` flags plus tissue semantics.
  2. dict_entries: the seven transmural ``ionicHeterogeneity.*`` DictEntries
     plus the five ``apexBaseBands.*`` entries (12 total, separately gated).
  3. dict_builder: build + parse round-trip of a heterogeneity block
     (proves the generic nested-path machinery needs no builder change).
  4. validation: model-capability gates, endo<mEpi ordering, apex-base
     numeric constraints, tissue compat.
"""

from __future__ import annotations

from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.core.specs.validation import validate_run
from openfoam_driver.core.plugin_interface import driver_context as _driver_context
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin

_CTX = _driver_context(CardiacFoamPlugin(), source="test")

_NATIVE_TISSUE_MODELS = ("BuenoOrovio", "TNNP", "TWorld", "ToRORd_dynCl")
_OVERRIDE_ONLY_TISSUE_MODELS = (
    "AlievPanfilov", "Courtemanche", "Fabbri", "Gaur",
    "Grandi", "PerisYague", "Stewart", "Trovato",
)


# --------------------------------------------------------------------------
# 1) Catalog flags
# --------------------------------------------------------------------------

def test_supports_heterogeneity_flag_for_capable_scalar_models():
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in ("BuenoOrovio", "TNNP", "TWorld", "ToRORd_dynCl"):
        assert IONIC_MODEL_CATALOG[name].supports_heterogeneity is True, name


def test_supports_heterogeneity_inherited_by_batched_variants():
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in (
        "BuenoOroviocompactBatched", "TNNPcompactBatched",
        "TWorldcompactBatched", "ToRORd_dynClcompactBatched",
    ):
        assert IONIC_MODEL_CATALOG[name].supports_heterogeneity is True, name


def test_single_tissue_models_do_not_support_heterogeneity():
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in (
        "monodomainFDAManufactured", "bidomainFDAManufactured",
        "bathBidomainFDAManufactured",
    ):
        assert IONIC_MODEL_CATALOG[name].supports_heterogeneity is False, name


def test_supports_apex_base_heterogeneity_for_capable_scalar_models():
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in ("BuenoOrovio", "TNNP", "TWorld", "ToRORd_dynCl"):
        assert IONIC_MODEL_CATALOG[name].supports_apex_base_heterogeneity is True, name


def test_supports_apex_base_heterogeneity_inherited_by_batched_variants():
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in (
        "BuenoOroviocompactBatched", "TNNPcompactBatched",
        "TWorldcompactBatched", "ToRORd_dynClcompactBatched",
    ):
        assert IONIC_MODEL_CATALOG[name].supports_apex_base_heterogeneity is True, name


def test_native_tissue_labels_mark_models_with_intrinsic_tissue_variants():
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    expected = ("epicardialCells", "mCells", "endocardialCells")
    for name in _NATIVE_TISSUE_MODELS:
        assert IONIC_MODEL_CATALOG[name].native_tissue_labels == expected, name
        assert IONIC_MODEL_CATALOG[name].approximate_tissue_labels == (), name


def test_override_only_models_advertise_approximate_tissue_labels_explicitly():
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    expected = ("epicardialCells", "mCells", "endocardialCells")
    for name in _OVERRIDE_ONLY_TISSUE_MODELS:
        assert IONIC_MODEL_CATALOG[name].native_tissue_labels == ("myocyte",), name
        assert IONIC_MODEL_CATALOG[name].approximate_tissue_labels == expected, name


def test_default_single_cell_tissue_map_uses_native_tissues_only():
    from openfoam_driver.plugins.cardiacfoam.tutorials.defaults.single_cell import IONIC_MODEL_TISSUE_MAP
    assert IONIC_MODEL_TISSUE_MAP["BuenoOrovio"] == (
        "epicardialCells", "mCells", "endocardialCells",
    )
    assert IONIC_MODEL_TISSUE_MAP["Gaur"] == ("myocyte",)


def test_default_restitution_tissue_map_uses_native_tissues_only():
    from openfoam_driver.plugins.cardiacfoam.tutorials.defaults.restitution_curves import IONIC_MODEL_TISSUE_MAP
    assert IONIC_MODEL_TISSUE_MAP["TNNP"] == (
        "epicardialCells", "mCells", "endocardialCells",
    )
    assert IONIC_MODEL_TISSUE_MAP["Courtemanche"] == ("myocyte",)


def test_transmural_only_models_do_not_support_apex_base_heterogeneity():
    # These models support neither transmural/named-region nor apex-base
    # heterogeneity at all — the manufactured verification models.
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG
    for name in (
        "monodomainFDAManufactured", "bidomainFDAManufactured",
        "bathBidomainFDAManufactured",
    ):
        assert IONIC_MODEL_CATALOG[name].supports_apex_base_heterogeneity is False, name


# --------------------------------------------------------------------------
# 2) dict_entries
# --------------------------------------------------------------------------

def _het_entries():
    from openfoam_driver.dict_entries import get_electro_property_entry_groups
    return get_electro_property_entry_groups()["ionic_heterogeneity"]


def _transmural_entries():
    return [
        e for e in _het_entries()
        if not e.driver_path.startswith("$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.apexBaseBands.")
    ]


def _apex_base_entries():
    return [
        e for e in _het_entries()
        if e.driver_path.startswith("$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.apexBaseBands.")
    ]


def test_all_twelve_heterogeneity_entries_exist():
    paths = {e.driver_path for e in _het_entries()}
    transmural_leaves = (
        "field", "mode", "endoMInterface", "mEpiInterface",
        "transitionWidth", "transitionMode", "smoothing",
    )
    ab_leaves = ("apexBaseBands.field", "apexBaseBands.beta",
                 "apexBaseBands.scalingMin", "apexBaseBands.scalingMax",
                 "apexBaseBands.variables")
    region_leaves = (
        "regions.<region_name>.baseline",
        "regions.<region_name>.cellZone",
        "regions.<region_name>.range",
    )
    expected = {
        f"$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.{leaf}"
        for leaf in (*transmural_leaves, *ab_leaves, *region_leaves)
    }
    assert paths == expected


def test_transmural_entries_gated_to_spatial_solvers():
    """Transmural ionicHeterogeneity entries must fire for all spatial EP solvers."""
    for e in _transmural_entries():
        assert e.applicable_when.get("$ionicHeterogeneity_supported") is True, e.driver_path


def test_apex_base_entries_gated_to_monodomain_and_bidomain():
    """apexBaseBands entries apply to monodomainSolver and bidomainSolver.

    Both dispatch through the same myocardiumDomainInterface::New() codepath
    (myocardiumDomainInterface.C) that parses ionicHeterogeneity.apexBaseBands;
    eikonalSolver returns early from that factory before ionic-model/heterogeneity
    setup runs, and singleCellSolver bypasses the factory entirely.
    """
    for e in _apex_base_entries():
        assert e.applicable_when.get("myocardiumSolver") == (
            "monodomainSolver", "bidomainSolver",
        ), e.driver_path


def test_heterogeneity_enum_values():
    by_leaf = {e.driver_path.rsplit(".", 1)[-1]: e for e in _het_entries()}
    assert set(by_leaf["mode"].enum_values) == {"transmuralBands", "namedRegions", "cellZoneRegions"}
    assert by_leaf["transitionMode"].enum_values == ("blend", "hard")
    assert by_leaf["smoothing"].enum_values == ("smoothstep",)


def test_apex_base_numeric_constraints():
    by_leaf = {e.driver_path.rsplit(".", 1)[-1]: e for e in _apex_base_entries()}
    assert any(">" in c for c in by_leaf["beta"].constraints)
    assert any("scalingMax" in c for c in by_leaf["scalingMin"].constraints)


def test_heterogeneity_entries_are_optional():
    for e in _het_entries():
        assert e.required is False


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
    from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
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
    from openfoam_driver.plugins.cardiacfoam.dict_builder import (
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
    from openfoam_driver.plugins.cardiacfoam.dict_builder import build_electro_properties
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
        "ionicModel": "monodomainFDAManufactured",
        "tissue": "manufactured",
        "ionicHeterogeneity.field": "t",
        "ionicHeterogeneity.mode": "transmuralBands",
    })
    errors = [e for e in validate_run(run, driver_context=_CTX)
              if e.level == "error" and "heterogeneity" in e.message.lower()]
    assert len(errors) == 1, [e.message for e in validate_run(run, driver_context=_CTX)]


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
    het_errors = [e for e in validate_run(run, driver_context=_CTX)
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
    errors = [e for e in validate_run(run, driver_context=_CTX)
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
    errors = [e for e in validate_run(run, driver_context=_CTX) if "endoMInterface" in e.message]
    assert errors == []


def test_tissue_incompatible_with_model_is_error():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "AlievPanfilovcompactBatched",   # myocyte-only (not yet wired for heterogeneity)
        "tissue": "epicardialCells",
    })
    errors = [e for e in validate_run(run, driver_context=_CTX)
              if e.level == "error" and "compatible tissues" in e.message]
    assert len(errors) == 1, [e.message for e in validate_run(run, driver_context=_CTX)]


def test_tissue_compatible_with_model_is_silent():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "BuenoOrovio",
        "tissue": "epicardialCells",
    })
    issues = [e for e in validate_run(run, driver_context=_CTX) if "compatible tissues" in e.message]
    assert issues == []


# --------------------------------------------------------------------------
# 5) Apex-to-base heterogeneity validation
# --------------------------------------------------------------------------

def test_apex_base_with_incapable_model_is_error():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "bidomainFDAManufactured",
        "tissue": "myocyte",
        "ionicHeterogeneity.apexBaseBands.field": "longitudinal",
        "ionicHeterogeneity.apexBaseBands.variables": "(g_Ks)",
    })
    errors = [e for e in validate_run(run, driver_context=_CTX)
              if e.level == "error" and "apex-to-base" in e.message]
    assert len(errors) == 1, [e.message for e in validate_run(run, driver_context=_CTX)]


def test_apex_base_with_capable_model_no_error():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.apexBaseBands.field": "longitudinal",
        "ionicHeterogeneity.apexBaseBands.beta": "3.0",
        "ionicHeterogeneity.apexBaseBands.scalingMin": "0.2",
        "ionicHeterogeneity.apexBaseBands.scalingMax": "5.0",
        "ionicHeterogeneity.apexBaseBands.variables": "(g_Ks)",
    })
    ab_errors = [e for e in validate_run(run, driver_context=_CTX)
                 if e.level == "error" and "apex" in e.message.lower()]
    assert ab_errors == []


def test_apex_base_beta_must_be_positive():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.apexBaseBands.field": "longitudinal",
        "ionicHeterogeneity.apexBaseBands.beta": "-1.0",
        "ionicHeterogeneity.apexBaseBands.variables": "(g_Ks)",
    })
    errors = [e for e in validate_run(run, driver_context=_CTX)
              if e.level == "error" and "beta" in e.message]
    assert len(errors) == 1


def test_apex_base_positive_beta_is_silent():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.apexBaseBands.field": "longitudinal",
        "ionicHeterogeneity.apexBaseBands.beta": "3.0",
        "ionicHeterogeneity.apexBaseBands.variables": "(g_Ks)",
    })
    errors = [e for e in validate_run(run, driver_context=_CTX) if "beta" in e.message]
    assert errors == []


def test_apex_base_scalingMin_must_not_exceed_scalingMax():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.apexBaseBands.field": "longitudinal",
        "ionicHeterogeneity.apexBaseBands.scalingMin": "8.0",
        "ionicHeterogeneity.apexBaseBands.scalingMax": "2.0",
        "ionicHeterogeneity.apexBaseBands.variables": "(g_Ks)",
    })
    errors = [e for e in validate_run(run, driver_context=_CTX)
              if e.level == "error" and "scalingMin" in e.message]
    assert len(errors) == 1


def test_apex_base_valid_scaling_range_is_silent():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.apexBaseBands.field": "longitudinal",
        "ionicHeterogeneity.apexBaseBands.scalingMin": "0.2",
        "ionicHeterogeneity.apexBaseBands.scalingMax": "5.0",
        "ionicHeterogeneity.apexBaseBands.variables": "(g_Ks)",
    })
    errors = [e for e in validate_run(run, driver_context=_CTX) if "scalingMin" in e.message]
    assert errors == []


def test_transmural_and_apex_base_can_coexist():
    run = _run({
        "myocardiumSolver": "monodomainSolver",
        "ionicModel": "TNNP",
        "tissue": "epicardialCells",
        "ionicHeterogeneity.field": "t",
        "ionicHeterogeneity.mode": "transmuralBands",
        "ionicHeterogeneity.endoMInterface": "0.3",
        "ionicHeterogeneity.mEpiInterface": "0.7",
        "ionicHeterogeneity.apexBaseBands.field": "longitudinal",
        "ionicHeterogeneity.apexBaseBands.beta": "3.0",
        "ionicHeterogeneity.apexBaseBands.variables": "(g_Ks)",
    })
    het_errors = [
        e for e in validate_run(run, driver_context=_CTX)
        if e.level == "error" and "heterogeneity" in e.message.lower()
    ]
    assert het_errors == []
