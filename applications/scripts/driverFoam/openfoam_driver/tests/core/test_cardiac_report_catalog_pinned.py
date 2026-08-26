"""cardiacFoam's report catalog must survive the move off legacy_report_catalog.

The next task implements ``get_report_catalog()`` on the plugin and deletes the
``plugin_id == "org.cardiacfoam"`` branch in ``core/compatibility.py``. The
catalog CONTENT must not change as a result -- only where core looks it up.

This matters because the four export conformance tests are being merged, so
``test_report_catalog_export.py`` cannot serve as the pin. A committed fixture
holds the line across both the merge and the fallback removal, which is exactly
the pair of changes that could silently drop a report.
"""
from __future__ import annotations

import json
import pathlib

from openfoam_driver.core.plugin_interface import driver_context
from openfoam_driver.core.report_catalog import to_record
from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin

_BASELINE = (
    pathlib.Path(__file__).resolve().parents[1]
    / "fixtures"
    / "cardiac_reports_baseline.json"
)


def test_cardiac_report_catalog_matches_baseline() -> None:
    ctx = driver_context(CardiacFoamPlugin(), source="test")
    reports = ctx.capabilities.report_catalog.reports()
    actual = {"version": "1", "reports": [to_record(r) for r in reports]}
    expected = json.loads(_BASELINE.read_text())
    assert actual == expected
