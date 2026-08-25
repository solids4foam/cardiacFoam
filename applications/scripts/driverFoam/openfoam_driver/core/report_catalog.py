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
#     report_catalog
#
# Description
#     Maintains metadata definitions for executable reports.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Solver-neutral report-catalog infrastructure.

After a run completes, downstream tools can select report definitions whose
``applicable_when`` predicate matches the run's configuration. This module
owns the shared machinery -- the ``ReportDefinition`` record, the
``applicable_when`` predicate evaluator, and the JSON record shape -- but
not any concrete catalog: which reports exist is solver-specific data owned
by the plugin that authors them (the built-in cardiac plugin's catalog lives
at ``plugins/cardiacfoam/reports.py``). ``scripts/export-report-catalog.py``
reaches the active plugin's catalog through
``driver_context.capabilities.report_catalog.reports()`` and serializes it
to JSON.

Two design choices worth re-reading later:

1. **URLs are templates, not absolute.** v1 ships
   ``http://localhost:{port}/{kind}`` for entries served by the user's
   standalone Quart app, **4Dpapers**, plus a bundled fallback at
   ``/reports/stub.html``. ``{port}`` and ``{kind}`` are substituted by the
   consumer. **No ``{runId}``** in v1 — 4Dpapers does not route by run yet
   (the user confirmed: "it does not read IDs yet but that is a simple
   implementation in the future"). When 4Dpapers learns to route by run, the
   template becomes ``.../{runId}/{kind}`` and JSON consumers pick that up
   without a code change.

2. **``applicable_when`` is flat key-equality only.** ``None`` ⇒ always
   applicable; ``{"phase.field": value}`` ⇒ AND of equality checks via
   shallow get on the run's resolved config. Anything richer
   (``$in``, ``$gt``, regex) is intentionally rejected so v2 can layer
   in a real predicate language without silently mis-filtering v1
   docs. The same predicate language is shared by tutorials' ``preset``
   field (Spec § Backend contract alignment).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping


# --- v1 URL templates -------------------------------------------------------

#: The 4Dpapers backend (user's existing standalone Quart app) does not
#: route by ``runId`` in v1. ``{port}`` is supplied by the consumer and
#: ``{kind}`` is filled from the report's ``id``.
URL_TEMPLATE = "http://localhost:{port}/{kind}"

#: Bundled offline fallback path.
STUB_URL = "/reports/stub.html"


# --- definition record ------------------------------------------------------


@dataclass(frozen=True)
class ReportDefinition:
    """One row of the report catalog.

    Fields are emitted directly by the report catalog exporter. Naming stays
    Pythonic (``snake_case``).
    """

    id: str
    title: str
    kind: str  # v1: only "iframe"
    url_template: str
    applicable_when: Mapping[str, Any] | None = None
    show_by_default: bool = True
    description: str = ""


# --- predicate evaluator (v1: flat key-equality) ---------------------------


def _shallow_get(cfg: Mapping[str, Any], dotted: str) -> Any | _Missing:
    """Resolve a dotted ``"a.b.c"`` path against a nested mapping.

    Returns ``MISSING`` if any segment is absent — a missing path is
    treated as "does not match", not as "matches None".
    """
    cur: Any = cfg
    for seg in dotted.split("."):
        if not isinstance(cur, Mapping) or seg not in cur:
            return MISSING
        cur = cur[seg]
    return cur


class _Missing:
    """Sentinel for "this path is absent from the config".

    Distinct from ``None`` because a config field explicitly set to
    ``None`` should still be considered present (and equal to ``None``
    for matching purposes).
    """

    _instance: "_Missing | None" = None

    def __new__(cls) -> "_Missing":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __repr__(self) -> str:  # pragma: no cover — debug aid only
        return "<MISSING>"


MISSING = _Missing()


def matches(
    predicate: Mapping[str, Any] | None,
    config: Mapping[str, Any],
) -> bool:
    """Evaluate the v1 ``applicable_when`` predicate against a config.

    Rules:
    - ``predicate is None`` ⇒ always matches.
    - Each entry must be ``"dotted.path": scalar``. AND across entries.
    - A value that is itself a mapping is treated as an *operator
      object* and rejected with ``ValueError`` — v2 will introduce
      operators (``$in``, ``$gt``, …) but v1 must fail loudly so a
      forward-compat doc never silently mis-filters.
    """
    if predicate is None:
        return True
    for path, expected in predicate.items():
        if isinstance(expected, Mapping):
            raise ValueError(
                f"unsupported predicate operator at {path!r}: "
                f"v1 applicable_when is flat key-equality only "
                f"(operators like $in / $gt land in v2)"
            )
        actual = _shallow_get(config, path)
        if actual is MISSING or actual != expected:
            return False
    return True


# --- helper used by the exporter and tests ---------------------------------


def to_record(r: ReportDefinition) -> dict:
    """Serialize one definition to the JSON record shape.

    Kept here (not in the export script) so tests can call it without
    spawning a subprocess and so any future programmatic consumer
    (a CI lint, a docs renderer) gets the same shape.
    """
    return {
        "id": r.id,
        "title": r.title,
        "kind": r.kind,
        "url_template": r.url_template,
        # ``None`` is a meaningful marker on the wire — keep it.
        "applicable_when": (
            None if r.applicable_when is None else dict(r.applicable_when)
        ),
        "show_by_default": r.show_by_default,
        "description": r.description,
    }
