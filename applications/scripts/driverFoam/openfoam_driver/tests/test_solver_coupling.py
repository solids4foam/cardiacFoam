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
#     test_solver_coupling
#
# Description
#     Tests solver coupling logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Contract tests for solver_coupling.SOLVER_COMPATIBILITY_RULES.

Locks the shape of the table so future edits do not silently drop a field
that downstream consumers (currently introspection.py, the LLM-agent
describe-tutorial payload) depend on.
"""
from __future__ import annotations

import unittest

from openfoam_driver.solver_coupling import SOLVER_COMPATIBILITY_RULES


class TestSolverCompatibilityRules(unittest.TestCase):
    _REQUIRED_KEYS: frozenset[str] = frozenset({
        "myocardium_solver",
        "purkinje_solver",
        "required_coupler",
        "valid",
    })

    def test_every_rule_has_required_keys(self) -> None:
        for rule in SOLVER_COMPATIBILITY_RULES:
            missing = self._REQUIRED_KEYS - set(rule)
            self.assertEqual(
                missing, set(),
                f"rule {rule!r} is missing required keys: {sorted(missing)}",
            )

    def test_invalid_rules_carry_a_reason(self) -> None:
        for rule in SOLVER_COMPATIBILITY_RULES:
            if not rule["valid"]:
                self.assertIn(
                    "reason", rule,
                    f"invalid rule {rule!r} must explain why it is invalid",
                )
                self.assertTrue(rule["reason"], "reason must not be empty")

    def test_valid_rules_have_a_required_coupler(self) -> None:
        for rule in SOLVER_COMPATIBILITY_RULES:
            if rule["valid"]:
                self.assertIsNotNone(
                    rule["required_coupler"],
                    f"valid rule {rule!r} must declare its required_coupler",
                )

    def test_backward_compat_reexport_from_catalog(self) -> None:
        """ionic_model_catalog.py still re-exports the rules for any consumer
        that imported them from there before the extraction. Removing the
        re-export is a breaking change."""
        from openfoam_driver.ionic_model_catalog import (
            SOLVER_COMPATIBILITY_RULES as catalog_rules,
        )
        self.assertIs(catalog_rules, SOLVER_COMPATIBILITY_RULES)


if __name__ == "__main__":
    unittest.main()
