"""Display metadata for exported tutorial catalogs.

The hardcoded backend tutorials live in
``openfoam_driver.core.runtime.registry.REGISTERED_TUTORIALS`` (factory
identifiers used to build a ``TutorialSpec``). Those identifiers are
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


# Decorate every entry in REGISTERED_TUTORIALS. Adding a new factory in
# ``core/runtime/registry.py`` without adding a row here will fail the
# A11 cross-check test, by design.
TUTORIALS: tuple[TutorialDisplay, ...] = (
    TutorialDisplay(
        id="singleCell",
        title="Single-cell action potential",
        summary=(
            "Run a single-cell sweep over an ionic model to inspect AP "
            "morphology. Useful for drug-effect studies."
        ),
        thumbnail="/tutorials/single-cell.png",
        tags=("single-cell", "ionic-model", "AP"),
        preset={
            "anatomy.mesh": "single-cell",
            "physics.ionic_model": "TenTusscher",
        },
    ),
    TutorialDisplay(
        id="niederer2012",
        title="Niederer 2012 verification benchmark",
        summary=(
            "The Niederer et al. 2012 N-version benchmark for cardiac "
            "tissue electrophysiology. Validates monodomain solvers."
        ),
        thumbnail="/tutorials/niederer-2012.png",
        tags=("benchmark", "verification", "monodomain"),
        preset={
            "anatomy.mesh": "niederer-slab",
            "physics.ionic_model": "TenTusscher",
        },
    ),
    TutorialDisplay(
        id="manufacturedFDA",
        title="Manufactured solution (monodomain)",
        summary=(
            "Method of manufactured solutions on the monodomain "
            "equation. Used to verify spatial/temporal convergence."
        ),
        thumbnail="/tutorials/manufactured-fda.png",
        tags=("manufactured-solution", "verification", "monodomain"),
        preset={
            "anatomy.mesh": "fda-cuboid",
            "physics.ionic_model": "FentonKarma",
        },
    ),
    TutorialDisplay(
        id="manufacturedFDABidomain",
        title="Manufactured solution (bidomain)",
        summary=(
            "Same MMS verification at bidomain resolution. Pairs with "
            "the monodomain variant for cross-formulation comparison."
        ),
        thumbnail="/tutorials/manufactured-fda-bidomain.png",
        tags=("manufactured-solution", "verification", "bidomain"),
        preset={
            "anatomy.mesh": "fda-cuboid",
            "physics.ionic_model": "FentonKarma",
        },
    ),
    TutorialDisplay(
        id="manufacturedFDABathBidomain",
        title="Manufactured solution (bath bidomain)",
        summary=(
            "FDA bidomain-with-bath manufactured solution with a grounded "
            "bath electrode and bath ECG potential verification."
        ),
        thumbnail="/tutorials/manufactured-fda-bath-bidomain.png",
        tags=("manufactured-solution", "verification", "bidomain", "bath-ecg"),
        preset={
            "anatomy.mesh": "fda-bath-cuboid",
            "physics.ionic_model": "FentonKarma",
        },
    ),
    TutorialDisplay(
        id="manufacturedEikonalECG",
        title="Manufactured solution (eikonal ECG)",
        summary=(
            "Manufactured eikonal activation-time verification with template "
            "surrogate ECG and quadrature ECG reference."
        ),
        thumbnail="/tutorials/manufactured-eikonal-ecg.png",
        tags=("manufactured-solution", "verification", "eikonal", "ecg"),
        preset={
            "anatomy.mesh": "unit-domain",
            "physics.ionic_model": "none",
        },
    ),
    TutorialDisplay(
        id="restitutionCurves",
        title="Restitution curves (S1–S2 protocol)",
        summary=(
            "S1–S2 pacing protocol that traces APD restitution. Useful "
            "for arrhythmia substrate studies."
        ),
        thumbnail="/tutorials/restitution-curves.png",
        tags=("single-cell", "S1-S2", "restitution"),
        preset={
            "anatomy.mesh": "single-cell",
            "physics.ionic_model": "TenTusscher",
            "stimulus.protocol": "s1s2",
        },
    ),
)


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
