"""The postprocessing package must carry no domain vocabulary.

Its role is plumbing -- schedule a case-owned analysis module, apply house plot
style, write tables. The analysis itself belongs to the case.
"""
from __future__ import annotations

import re
import unittest
from pathlib import Path


_DOMAIN_WORDS = re.compile(
    r"cardiac|electro|ionic|monodomain|bidomain|eikonal|purkinje|myocard|"
    r"transmembrane|phiE|activation|stimulus|conductiv",
    re.IGNORECASE,
)


def _executable_lines(path: Path) -> list[tuple[int, str]]:
    """Source lines outside the licence header and outside docstrings."""
    lines = path.read_text().splitlines()
    out: list[tuple[int, str]] = []
    in_doc = False
    for number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if stripped.count('"""') == 1:
            in_doc = not in_doc
            continue
        if in_doc or stripped.startswith('"""'):
            continue
        out.append((number, line))
    return out


class TestPlottingCommonIsGeneric(unittest.TestCase):
    def test_no_domain_vocabulary_in_executable_code(self) -> None:
        from openfoam_driver.postprocessing import plotting_common
        offenders = [
            f"{number}: {line.strip()}"
            for number, line in _executable_lines(Path(plotting_common.__file__))
            if _DOMAIN_WORDS.search(line)
        ]
        self.assertEqual(offenders, [], f"domain vocabulary found: {offenders}")


class TestParseTwoPartStem(unittest.TestCase):
    def test_maps_both_tokens(self) -> None:
        from openfoam_driver.postprocessing.plotting_common import parse_two_part_stem
        first, second = parse_two_part_stem(
            "Alpha_Beta_extra.csv",
            first_map={"Alpha": "A"},
            second_map={"Beta": "B"},
        )
        self.assertEqual((first, second), ("A", "B"))

    def test_unmapped_tokens_pass_through(self) -> None:
        from openfoam_driver.postprocessing.plotting_common import parse_two_part_stem
        self.assertEqual(
            parse_two_part_stem("X_Y.csv", first_map={}, second_map={}),
            ("X", "Y"),
        )

    def test_missing_second_token_falls_back(self) -> None:
        from openfoam_driver.postprocessing.plotting_common import parse_two_part_stem
        first, second = parse_two_part_stem("Only.csv", first_map={}, second_map={})
        self.assertEqual(first, "Only")
        self.assertEqual(second, "UnknownSecond")


if __name__ == "__main__":
    unittest.main()
