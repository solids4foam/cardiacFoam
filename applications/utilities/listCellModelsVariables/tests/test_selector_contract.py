"""Static contract checks for listCellModelsVariables dictionary selection."""

from pathlib import Path
import unittest


SOURCE = Path(__file__).parents[1] / "listCellModelsVariables.C"


class SelectorContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text(encoding="utf-8")
        start = cls.source.index("dictionary electroModelDict")
        end = cls.source.index("dictionary activeTensionDict", start)
        cls.selector = cls.source[start:end]

    def test_canonical_selector_precedes_legacy_selector(self):
        canonical = self.selector.index('found("myocardiumSolver")')
        legacy = self.selector.index('found("electroModel")')
        self.assertLess(canonical, legacy)
        self.assertIn('myocardiumSolverName + "Coeffs"', self.selector)

    def test_legacy_selector_keeps_matching_coefficients_path(self):
        legacy = self.selector.index('found("electroModel")')
        legacy_path = self.selector[legacy:]
        self.assertIn('electroModelName + "Coeffs"', legacy_path)
        self.assertIn("return electroDict.subDict(coeffsName);", legacy_path)

    def test_missing_selector_diagnostic_uses_current_contract(self):
        self.assertIn(
            '<< "Expected \'myocardiumSolver\' in electroProperties."',
            self.selector,
        )
        self.assertNotIn("flat dictionary", self.selector)


if __name__ == "__main__":
    unittest.main()
