"""A non-cardiac plugin must never inherit cardiacFoam semantics.

``core/compatibility.py`` is the one place core is allowed to import the
cardiac plugin (enforced by
``test_plugin_dependency_boundary.py::test_core_imports_cardiac_implementation_only_at_compatibility_boundary``).
That import-level test cannot see *which plugin* the fallback is answering
for, though: a fallback may legally live in ``compatibility.py`` and still
route every caller into cardiac code.

Thirteen fallbacks gate on ``plugin_id == "org.cardiacfoam"`` and hand a
neutral value to everyone else. Six reachable through a capability adapter did
not, so the shipped ``GenericOpenFOAMPlugin`` -- which implements only one of
the six hooks -- answered cardiac questions about non-cardiac cases. The worst
of them wrote an ``Allrun`` invoking the ``cardiacFoam`` binary for a sweep
under a plugin that is not cardiacFoam.

These tests are behavioural on purpose. Asserting "the cardiac module was not
imported" is unreliable once any other test has imported it; asserting that a
deliberately cardiac-*looking* case is NOT claimed by the generic plugin can
only pass if the cardiac code path genuinely did not run.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from openfoam_driver.core.generic_plugin import GenericOpenFOAMPlugin
from openfoam_driver.core.plugin_capabilities import (
    CaseCompatibilityRequest,
    SweepMaterializationRequest,
    SweepRoutingRequest,
)
from openfoam_driver.core.plugin_interface import driver_context
from openfoam_driver.core.sweep.sweep_expansion import SweepValidationError


def _generic_capabilities():
    context = driver_context(
        GenericOpenFOAMPlugin(), source="test:generic-openfoam",
    )
    return context.capabilities


def _cardiac_looking_case() -> Path:
    """A case carrying every marker the cardiac fallbacks look for.

    If a cardiac fallback answers for the generic plugin, this case is
    claimed. If the fallback is correctly gated, it is not.
    """
    case_root = Path(tempfile.mkdtemp())
    constant = case_root / "constant"
    system = case_root / "system"
    constant.mkdir()
    system.mkdir()
    (constant / "electroProperties").write_text(
        "myocardiumSolver singleCellSolver;\n"
        "singleCellSolverCoeffs\n{\n    ionicModel TNNP;\n}\n"
    )
    (constant / "physicsProperties").write_text("physics electrophysiology;\n")
    for name in ("controlDict", "fvSchemes", "fvSolution"):
        (system / name).write_text("// placeholder\n")
    return case_root


class TestGenericPluginDoesNotInheritCardiacSemantics(unittest.TestCase):
    def test_does_not_claim_a_cardiac_looking_case(self) -> None:
        capabilities = _generic_capabilities()
        request = CaseCompatibilityRequest(case_root=_cardiac_looking_case())

        self.assertFalse(
            capabilities.case_compatibility.has_case_marker(request),
            "the generic plugin must not claim a case just because it has "
            "constant/electroProperties -- that is cardiac evidence",
        )

    def test_does_not_declare_a_cardiac_looking_case_runnable(self) -> None:
        capabilities = _generic_capabilities()
        request = CaseCompatibilityRequest(case_root=_cardiac_looking_case())

        self.assertFalse(
            capabilities.case_compatibility.is_runnable_without_workflow(request),
            "runnability for the generic plugin must not be decided by the "
            "presence of cardiac dictionaries",
        )

    def test_sweep_routing_refuses_rather_than_routing_cardiac_axes(self) -> None:
        capabilities = _generic_capabilities()
        request = SweepRoutingRequest(base={}, resolved_axis_values={"dt": 0.01})

        with self.assertRaises(SweepValidationError) as caught:
            capabilities.sweep_materializer.route(request, driver_context=None)

        message = str(caught.exception)
        self.assertIn("route_sweep_case_values", message)

    def test_sweep_materialization_refuses_rather_than_writing_cardiacfoam(self) -> None:
        """The defect with teeth: a wrong ``Allrun``, not a safe refusal."""
        capabilities = _generic_capabilities()
        case_dir = Path(tempfile.mkdtemp())
        request = SweepMaterializationRequest(case_dir=case_dir, routed={})

        with self.assertRaises(SweepValidationError) as caught:
            capabilities.sweep_materializer.materialize(request)

        self.assertIn("materialize_sweep_case", str(caught.exception))
        self.assertFalse(
            (case_dir / "Allrun").exists(),
            "a refused materialization must leave no Allrun behind",
        )


class TestCardiacPluginBehaviourIsUnchanged(unittest.TestCase):
    """Gating must be invisible to cardiacFoam, which implements every hook."""

    def _cardiac_capabilities(self):
        from openfoam_driver.plugins.cardiacfoam_plugin import CardiacFoamPlugin

        context = driver_context(CardiacFoamPlugin(), source="test:cardiacfoam")
        return context.capabilities

    def test_cardiac_plugin_still_claims_a_cardiac_case(self) -> None:
        capabilities = self._cardiac_capabilities()
        request = CaseCompatibilityRequest(case_root=_cardiac_looking_case())

        self.assertTrue(capabilities.case_compatibility.has_case_marker(request))

    def test_cardiac_plugin_still_declares_it_runnable(self) -> None:
        capabilities = self._cardiac_capabilities()
        request = CaseCompatibilityRequest(case_root=_cardiac_looking_case())

        self.assertTrue(
            capabilities.case_compatibility.is_runnable_without_workflow(request)
        )


if __name__ == "__main__":
    unittest.main()
