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
        id="manufacturedMonodomainTotalLagrangianEM",
        title="Manufactured electromechanics (MMS)",
        summary=(
            "Electromechanics verification on a fully coupled manufactured field. "
            "Vm, D, lambda, and Ta are rigorous MMS targets."
        ),
        thumbnail="/tutorials/manufactured-electromechanics-bc.png",
        tags=("manufactured-solution", "verification", "electromechanics"),
        preset={
            "anatomy.mesh": "unit-domain",
            "physics.ionic_model": "monodomainFDAManufactured",
        },
    ),
    TutorialDisplay(
        id="manufacturedPurkinjeGraph",
        title="Manufactured solution (Purkinje graph)",
        summary=(
            "Manufactured monodomain solution on a 1D Purkinje graph coupled "
            "to a 3D domain. Traces mesh refinement convergence on the "
            "Hines-ordered network."
        ),
        thumbnail="/tutorials/manufactured-purkinje-graph.png",
        tags=("manufactured-solution", "verification", "purkinje", "1D-3D"),
        preset={
            "anatomy.mesh": "purkinje-graph",
            "physics.ionic_model": "monodomainFDAManufactured",
        },
    ),
    TutorialDisplay(
        id="heartSolverComparison",
        title="Heart solver comparison (eikonal / monodomain / bidomain)",
        summary=(
            "Compares whole solver stacks -- eikonal, monodomain, bidomain, "
            "and a mixed monodomain-tissue/eikonal-Purkinje variant -- over "
            "one shared real heart anatomy (mesh + Purkinje graph)."
        ),
        thumbnail="/tutorials/heart-solver-comparison.png",
        tags=("real-anatomy", "purkinje", "eikonal", "monodomain", "bidomain", "solver-comparison"),
        preset={
            "anatomy.mesh": "heart-purkinje-graph",
            # No single physics.ionic_model: eikonalSolver has no ionic
            # model at all (pure activation-time PDE); only the monodomain,
            # monodomain-eikonal, and bidomain variants use BuenoOrovio.
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
    TutorialDisplay(
        id="monodomainAndEikonal1DCableCVConvergence",
        title="1D Cable CV Convergence (Monodomain & Eikonal)",
        summary=(
            "1D cable verification protocol to extract continuous conduction "
            "velocity profiles and perform mesh resolution convergence sweeps."
        ),
        thumbnail="/tutorials/cable-cv-convergence.png",
        tags=("cable", "cv", "convergence", "monodomain", "eikonal"),
        preset={
            "anatomy.mesh": "cable-1d",
            "physics.ionic_model": "BuenoOrovio",
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
