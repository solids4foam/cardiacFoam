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
#     reports
#
# Description
#     cardiacFoam's post-run report catalog: which report definitions the
#     built-in plugin offers, expressed with the solver-neutral
#     ``ReportDefinition`` record owned by ``openfoam_driver.report_catalog``.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""cardiacFoam's report catalog.

Moved out of ``openfoam_driver.report_catalog`` (P2.7): that module owns the
solver-neutral machinery (``ReportDefinition``, the ``applicable_when``
predicate evaluator, the JSON record shape); the actual list of reports is
solver-specific data and belongs with the plugin that authored it. Reached
through ``driver_context.capabilities.report_catalog.reports()`` --
``core/compatibility.py``'s ``legacy_report_catalog`` returns this tuple only
for ``org.cardiacfoam`` and ``()`` for every other plugin.
"""

from __future__ import annotations

from ...report_catalog import ReportDefinition, STUB_URL, URL_TEMPLATE

CARDIAC_REPORTS: tuple[ReportDefinition, ...] = (
    ReportDefinition(
        id="vm-field-3d",
        title="Vm field (3D)",
        kind="iframe",
        url_template=URL_TEMPLATE,
        applicable_when=None,  # always available post-completion
        show_by_default=True,
        description=(
            "Volumetric Vm field rendered by 4Dpapers from the run's "
            "foam/VTK output."
        ),
    ),
    ReportDefinition(
        id="activation-map",
        title="Activation map",
        kind="iframe",
        url_template=URL_TEMPLATE,
        applicable_when=None,
        show_by_default=True,
        description=(
            "Local activation time map. Useful for checking conduction "
            "patterns and reentry."
        ),
    ),
    ReportDefinition(
        id="stub",
        title="Stub (4Dpapers not running)",
        kind="iframe",
        url_template=STUB_URL,
        applicable_when=None,
        show_by_default=False,
        description=(
            "Bundled fallback that renders when the 4Dpapers backend is "
            "not reachable. Visible for diagnostics, not by default."
        ),
    ),
)
