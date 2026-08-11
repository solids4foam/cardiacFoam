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
#     test_override_schema_capability
#
# Description
#     Override prose, examples, and document naming are plugin-authored.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""The config-schema prose and dict-entry document shape are plugin knowledge.

Core assembles and serializes; the plugin supplies the vocabulary. A generic
OpenFOAM plugin must therefore emit no cardiacFoam token at all.
"""

from __future__ import annotations

import json

from openfoam_driver.core.plugin_interface import (
    default_driver_context,
    generic_openfoam_context,
)

# Every token that would betray cardiac vocabulary leaking into a generic
# plugin's machine-readable schema.
_CARDIAC_TOKENS = (
    "ELECTRO_MODEL_COEFFS",
    "electroProperties",
    "physicsProperties",
    "ionicModel",
    "monodomainSolverCoeffs",
    "TNNP",
)

_MAKE_SPEC_INFO = {"parameters": {"ionic_models": {"default": ["TNNP"]}}}


def test_generic_config_schema_mentions_no_cardiac_tokens() -> None:
    schema = generic_openfoam_context().capabilities.override_schema.config_schema(
        "someTutorial", _MAKE_SPEC_INFO
    )
    blob = json.dumps(schema)
    leaked = [token for token in _CARDIAC_TOKENS if token in blob]
    assert leaked == [], f"generic config schema leaked cardiac tokens: {leaked}"


def test_cardiac_config_schema_keeps_its_documented_tokens() -> None:
    schema = default_driver_context().capabilities.override_schema.config_schema(
        "singleCell", _MAKE_SPEC_INFO
    )
    blob = json.dumps(schema)
    assert "$ELECTRO_MODEL_COEFFS" in blob
    assert "electroProperties" in blob
    # The worked example is tutorial-specific, so the name must be threaded
    # through rather than baked into the plugin.
    assert "singleCell" in schema["worked_example"]["json"]


def test_generic_dict_entry_catalog_names_no_cardiac_document() -> None:
    """The previous version of this test asserted only on ``.values()``, so a
    cardiac leak in the *keys* (``physicsProperties``/``electroProperties``)
    was invisible to it. Scan the whole structure."""
    catalog = generic_openfoam_context().capabilities.override_schema.dict_entry_catalog()
    blob = json.dumps(catalog)
    leaked = [token for token in _CARDIAC_TOKENS if token in blob]
    assert leaked == [], f"generic dict entry catalog leaked: {leaked}"
    assert all(not value for value in catalog.values()), catalog


def test_cardiac_dict_entry_catalog_keeps_its_document_shape() -> None:
    catalog = default_driver_context().capabilities.override_schema.dict_entry_catalog()
    assert "physicsProperties" in catalog
    assert "electroProperties" in catalog
    # physicsProperties is a flat sequence; electroProperties is grouped. That
    # asymmetry is cardiac document knowledge, not a core convention.
    assert isinstance(catalog["electroProperties"], dict)
    assert catalog["electroProperties"], "cardiac plugin must expose entry groups"
