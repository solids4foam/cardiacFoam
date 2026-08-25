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
#     tutorials_display
#
# Description
#     Formats and displays metadata for exported tutorial catalogs.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Display metadata for exported tutorial catalogs.

The hardcoded backend tutorials live in
``openfoam_driver.plugins.cardiacfoam.tutorials.registry.REGISTERED_TUTORIALS``
(factory identifiers used to build a ``TutorialSpec``). Those identifiers are
fine for the CLI but unhelpful for an end-user reading a catalog. They need a
title, a one-line summary, a thumbnail, some tags, and a preset that walks the
new run into the right configuration.

This module is the thin display layer around the backend's factory registry.
It does NOT add new tutorials; it just decorates the existing ones for human
display.
The exporter cross-checks one-to-one against ``REGISTERED_TUTORIALS``
so a developer cannot ship a card without a backend factory or omit a
registered tutorial from the home page.

If we later move tutorials into a JSON-authored format, only this module and
the exporter are replaced.

The ``preset`` field shape mirrors the v1 ``applicable_when`` predicate
language used by ``report_catalog`` (flat ``"phase.field": value``
keys, no operator objects).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class TutorialDisplay:
    """One row of the home-page tutorial strip."""

    id: str  # MUST match an entry in REGISTERED_TUTORIALS
    title: str
    summary: str
    thumbnail: str
    tags: tuple[str, ...] = ()
    preset: dict[str, Any] = field(default_factory=dict)





def to_record(t: TutorialDisplay) -> dict:
    """Serialize one display entry to the JSON record shape."""
    return {
        "id": t.id,
        "title": t.title,
        "summary": t.summary,
        "thumbnail": t.thumbnail,
        "tags": list(t.tags),
        "preset": dict(t.preset),
    }
