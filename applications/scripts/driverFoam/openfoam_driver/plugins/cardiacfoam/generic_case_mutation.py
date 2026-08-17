"""cardiacFoam mutation hook for the legacy generic-case specification.

Core's generic-case factory addresses dictionary files through the neutral
``dict_file_relpaths``/``dict_file_overrides`` mappings. This module is the one
place that maps those neutral names back onto cardiacFoam's own documents.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from openfoam_driver.plugins.cardiacfoam.overrides import (
    apply_electro_property_overrides,
    apply_physics_property_overrides,
)

ELECTRO_DICT_FILE = "electro"
PHYSICS_DICT_FILE = "physics"

DEFAULT_ELECTRO_PROPERTIES_RELPATH = Path("constant/electroProperties")
DEFAULT_PHYSICS_PROPERTIES_RELPATH = Path("constant/physicsProperties")


def apply_case_mutation(
    case_root: Path,
    case,
    *,
    dict_file_relpaths: Mapping[str, str | Path] | None = None,
    dict_file_overrides: Mapping[str, Any] | None = None,
) -> None:
    # The overrides now arrive explicitly rather than being dug out of
    # ``case.params``; the case is still part of the callback contract.
    del case
    relpaths = dict(dict_file_relpaths or {})
    overrides = dict(dict_file_overrides or {})

    electro_properties_relpath = Path(
        relpaths.get(ELECTRO_DICT_FILE, DEFAULT_ELECTRO_PROPERTIES_RELPATH)
    )
    physics_properties_relpath = Path(
        relpaths.get(PHYSICS_DICT_FILE, DEFAULT_PHYSICS_PROPERTIES_RELPATH)
    )

    apply_electro_property_overrides(
        case_root / electro_properties_relpath,
        overrides.get(ELECTRO_DICT_FILE),
    )
    apply_physics_property_overrides(
        case_root / physics_properties_relpath,
        overrides.get(PHYSICS_DICT_FILE),
    )
