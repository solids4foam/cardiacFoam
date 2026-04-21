from __future__ import annotations

import unittest

from openfoam_driver.gui_schema import describe_gui_schema


class TestGuiSchema(unittest.TestCase):
    def test_gui_schema_contains_expected_routes(self) -> None:
        payload = describe_gui_schema()
        route_paths = {route["path"] for route in payload["routes"]}
        self.assertEqual(
            route_paths,
            {
                "/entries",
                "/entries/:entryId",
                "/entries/:entryId/config",
                "/entries/:entryId/cases",
                "/runs",
                "/runs/:runId",
            },
        )

    def test_gui_schema_contains_core_view_models(self) -> None:
        payload = describe_gui_schema()
        view_model_ids = {view_model["id"] for view_model in payload["view_models"]}
        self.assertTrue(
            {
                "EntryCatalogViewModel",
                "WorkflowFamilyViewModel",
                "EntryOverviewViewModel",
                "SpecParameterViewModel",
                "DictEntryViewModel",
                "PlannedCaseViewModel",
                "LaunchRunRequest",
                "RunSummaryViewModel",
                "RunCaseResultViewModel",
            }.issubset(view_model_ids)
        )


if __name__ == "__main__":
    unittest.main()
