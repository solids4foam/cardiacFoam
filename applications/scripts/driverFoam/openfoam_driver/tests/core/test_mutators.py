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
#     test_mutators
#
# Description
#     Tests mutators logic and specification contracts.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from openfoam_driver.core.runtime.mutators import (
    ensure_foam_dict,
    read_foam_entry,
    remove_foam_dict,
    update_foam_entry,
    update_foam_entry_via_foamDictionary,
)
from openfoam_driver.tests.conftest import assert_foam_entry


def assert_entry_present(testcase: unittest.TestCase, text: str, key: str, value: str) -> None:
    """Assert `key <value>;` appears in `text`, tolerant of the column
    alignment foamDictionary applies when it re-serializes a whole file
    (e.g. `keep 1;` becomes `keep            1;`)."""
    pattern = rf"{re.escape(key)}\s+{re.escape(value)};"
    testcase.assertRegex(text, pattern)


class TestScopedMutators(unittest.TestCase):
    def test_updates_only_within_scope(self) -> None:
        text = "\n".join(
            [
                "myocardiumSolver monodomainSolver;",
                "",
                "monodomainSolverCoeffs",
                "{",
                "    ionicModel TNNP;",
                "}",
                "",
                "singleCellSolverCoeffs",
                "{",
                "    ionicModel BuenoOrovio;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            update_foam_entry(
                path,
                "ionicModel",
                "Gaur",
                scope="singleCellSolverCoeffs",
            )

            assert_foam_entry(
                path, "ionicModel", "TNNP", scope="monodomainSolverCoeffs"
            )
            assert_foam_entry(
                path, "ionicModel", "Gaur", scope="singleCellSolverCoeffs"
            )

    def test_nested_scope_path(self) -> None:
        text = "\n".join(
            [
                "outer",
                "{",
                "    inner",
                "    {",
                "        target 1;",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            update_foam_entry(path, "target", 2, scope=("outer", "inner"))
            assert_foam_entry(path, "target", "2", scope=("outer", "inner"))

    def test_quoted_regex_style_scope_name_is_matched(self) -> None:
        # OpenFOAM's fvSolution commonly names a solver block with a quoted
        # alternation, e.g. "phiE|phiEFinal|phiI|phiIFinal" { ... } -- the
        # scope-boundary regex's old trailing \b failed to match here because
        # both the character before and after the closing quote are
        # non-word characters, so there is no word boundary at all at that
        # position (verified: re.match(r'^\s*"foo"\b', '    "foo"\n') is
        # None). This is the quoted equivalent of test_nested_scope_path.
        # Forced off foamDictionary (which would mask the regex bug, since
        # it understands its own dictionary syntax natively) to test the
        # fallback parser's scope matching specifically.
        text = "\n".join(
            [
                "solvers",
                "{",
                '    "phiE|phiEFinal|phiI|phiIFinal"',
                "    {",
                "        tolerance 1e-06;",
                "    }",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "fvSolution"
            path.write_text(text)

            with mock.patch(
                "openfoam_driver.core.runtime.mutators.shutil.which",
                return_value=None,
            ):
                update_foam_entry(
                    path, "tolerance", 1e-15,
                    scope=("solvers", '"phiE|phiEFinal|phiI|phiIFinal"'),
                )
            self.assertIn("tolerance    1e-15;", path.read_text())

    def test_missing_scope_raises(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text("a { b 1; }\n")

            with self.assertRaises(KeyError):
                update_foam_entry(path, "b", 2, scope="missing")

    def test_python_parser_fails_on_c_style_comments(self) -> None:
        text = "\n".join(
            [
                "someDict",
                "{",
                "    /* This is a block comment with a brace { inside it */",
                "    value 1;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            # The fallback update_foam_entry uses brace counting, so the extra
            # { inside the block comment throws off the parser, causing it to
            # incorrectly raise a KeyError for unbalanced braces. Forced off
            # foamDictionary, which parses /* */ natively and would mask the
            # limitation being pinned here.
            with mock.patch(
                "openfoam_driver.core.runtime.mutators.shutil.which",
                return_value=None,
            ), self.assertRaises(KeyError):
                update_foam_entry(path, "value", 2, scope="someDict")

    def test_remove_foam_dict_removes_nested_dictionary(self) -> None:
        text = "\n".join(
            [
                "outer",
                "{",
                "    keep 1;",
                "    removeMe",
                "    {",
                "        nested",
                "        {",
                "            value 1;",
                "        }",
                "    }",
                "    after 2;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            remove_foam_dict(path, "removeMe", scope="outer")

            updated = path.read_text()
            assert_entry_present(self, updated, "keep", "1")
            assert_entry_present(self, updated, "after", "2")
            self.assertNotIn("removeMe", updated)
            self.assertNotIn("value 1;", updated)

    def test_remove_foam_dict_missing_ok_tolerates_absent_scope_without_foamDictionary(
        self,
    ) -> None:
        # The foamDictionary-backed path already returns cleanly on a missing
        # scope when missing_ok=True (remove_foam_dict_via_foamDictionary's
        # subprocess-failure branch checks it). The pure-Python fallback used
        # when foamDictionary isn't on PATH previously called
        # _resolve_search_region(lines, scope) unguarded, so it raised
        # KeyError("Scope '<name>' not found") before missing_ok was ever
        # consulted -- missing_ok only guarded the "dict_name not found
        # inside an existing scope" case, not "scope itself absent".
        text = "\n".join(
            [
                "outer",
                "{",
                "    keep 1;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            with mock.patch.object(shutil, "which", return_value=None):
                remove_foam_dict(
                    path,
                    "xMin",
                    scope=["outer", "neverExisted"],
                    missing_ok=True,
                )

                self.assertEqual(path.read_text(), text)

                with self.assertRaises(KeyError):
                    remove_foam_dict(
                        path,
                        "xMin",
                        scope=["outer", "neverExisted"],
                        missing_ok=False,
                    )

    def test_ensure_foam_dict_inserts_missing_dict_in_scope(self) -> None:
        text = "\n".join(
            [
                "root",
                "{",
                "    existing yes;",
                "}",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "dict"
            path.write_text(text)

            inserted = ensure_foam_dict(
                path,
                "newBlock",
                "    newBlock\n    {\n        value 1;\n    }\n",
                scope="root",
            )
            inserted_again = ensure_foam_dict(
                path,
                "newBlock",
                "    newBlock\n    {\n        value 2;\n    }\n",
                scope="root",
            )

            updated = path.read_text()
            self.assertTrue(inserted)
            self.assertFalse(inserted_again)
            self.assertEqual(updated.count("newBlock"), 1)
            assert_entry_present(self, updated, "value", "1")


class TestScopeDoesNotDescendIntoNestedDicts(unittest.TestCase):
    """A scope names one dictionary, not that dictionary and everything under
    it. The pure-Python fallback resolves a scope to a line *span* and then
    scans it, which without a depth check also matches keys belonging to
    nested sub-dictionaries -- so a read scoped to the parent returned a
    child's value, and a write scoped to the parent silently edited the
    child. foamDictionary is path-exact and does neither, so this divergence
    only appeared when OpenFOAM was sourced.

    These tests pin the Python implementation directly (foamDictionary forced
    absent) because it is the side that was wrong.
    """

    NESTED = "\n".join(
        [
            "monodomainSolverCoeffs",
            "{",
            "    ionicModel TNNP;",
            "    externalStimulus",
            "    {",
            "        stimulusIntensity 50000;",
            "    }",
            "}",
            "",
        ]
    )

    def setUp(self) -> None:
        patcher = mock.patch(
            "openfoam_driver.core.runtime.mutators.shutil.which",
            return_value=None,
        )
        patcher.start()
        self.addCleanup(patcher.stop)

    def _write(self, temp_dir: str) -> Path:
        path = Path(temp_dir) / "electroProperties"
        path.write_text(self.NESTED)
        return path

    def test_read_scoped_to_parent_ignores_nested_key(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = self._write(temp_dir)

            # stimulusIntensity is a child of externalStimulus, NOT of
            # monodomainSolverCoeffs -- so this scope has no such key.
            self.assertIsNone(
                read_foam_entry(
                    path, "stimulusIntensity", scope="monodomainSolverCoeffs"
                )
            )

    def test_write_scoped_to_parent_refuses_nested_key(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = self._write(temp_dir)

            with self.assertRaises(KeyError):
                update_foam_entry(
                    path, "stimulusIntensity", 99999, scope="monodomainSolverCoeffs"
                )

            # and the nested value must be left untouched
            self.assertIn("stimulusIntensity 50000;", path.read_text())

    def test_direct_child_of_scope_still_resolves(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = self._write(temp_dir)

            self.assertEqual(
                read_foam_entry(path, "ionicModel", scope="monodomainSolverCoeffs"),
                "TNNP",
            )
            update_foam_entry(
                path, "ionicModel", "Gaur", scope="monodomainSolverCoeffs"
            )
            self.assertIn("ionicModel    Gaur;", path.read_text())

    GRANDCHILD = "\n".join(
        [
            "monodomainSolverCoeffs",
            "{",
            "    outputVariables",
            "    {",
            "        ionic",
            "        {",
            "            export (Vm Jsi);",
            "        }",
            "    }",
            "}",
            "",
        ]
    )

    def test_scope_path_may_not_skip_an_intermediate_dict(self) -> None:
        # 'ionic' is a grandchild of monodomainSolverCoeffs, reachable only
        # through outputVariables -- a scope path that omits that level names
        # a dictionary which does not exist.
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(self.GRANDCHILD)
            skipping = ["monodomainSolverCoeffs", "ionic"]

            self.assertIsNone(read_foam_entry(path, "export", scope=skipping))
            with self.assertRaises(KeyError):
                update_foam_entry(path, "export", "(Vm)", scope=skipping)

            # the fully-qualified path still works
            full = ["monodomainSolverCoeffs", "outputVariables", "ionic"]
            self.assertEqual(read_foam_entry(path, "export", scope=full), "(Vm Jsi)")

    def test_nested_key_still_reachable_via_full_scope_path(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = self._write(temp_dir)
            scope = ["monodomainSolverCoeffs", "externalStimulus"]

            self.assertEqual(
                read_foam_entry(path, "stimulusIntensity", scope=scope), "50000"
            )
            update_foam_entry(path, "stimulusIntensity", 75000, scope=scope)
            self.assertIn("stimulusIntensity    75000;", path.read_text())


class TestReadFoamEntryIsEnvironmentIndependent(unittest.TestCase):
    """Reading a dict must not depend on whether OpenFOAM is sourced.

    foamDictionary parses each value into a double and re-serialises it, so
    reading through it respells the source text (``0.0`` -> ``0``,
    ``5.5e-3`` -> ``0.0055``, ``(a b c)`` -> ``( a b c )``). Those respelt
    values flow into build_electro_properties, which made generated dicts --
    and therefore run documents and provenance digests -- differ between a
    sourced and an unsourced shell.

    Reading through foamDictionary also *evaluates* the dictionary: a
    ``#calc`` / ``#codeStream`` entry is compiled, linked and executed to
    produce the value. Reading a case must never run code, least of all
    because override values are written verbatim into these dicts.

    These tests hold in either environment; before the fix the second one
    also left a ``dynamicCode/`` build directory behind.
    """

    def test_returns_the_literal_spelling_from_the_file(self) -> None:
        text = "\n".join(
            [
                "coeffs",
                "{",
                "    activationThreshold 0.0;",
                "    stimulusLocationMin (0 0 5.5e-3);",
                "    initialODEStep 1e-6;",
                "}",
                "",
            ]
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "electroProperties"
            path.write_text(text)

            self.assertEqual(
                read_foam_entry(path, "activationThreshold", scope="coeffs"), "0.0"
            )
            self.assertEqual(
                read_foam_entry(path, "stimulusLocationMin", scope="coeffs"),
                "(0 0 5.5e-3)",
            )
            self.assertEqual(
                read_foam_entry(path, "initialODEStep", scope="coeffs"), "1e-6"
            )

    def test_resolves_a_scope_written_as_an_inline_block(self) -> None:
        # `solvers { V { tolerance 1e-5; } }` is legal OpenFOAM and appears in
        # 10 tracked tutorial dicts. The scope machinery works on whole lines,
        # so an inline block used to resolve to a degenerate range and read as
        # None -- previously masked because foamDictionary parsed these.
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "fvSolution"
            path.write_text(
                "solvers { V { tolerance 1e-5; } p { tolerance 1e-7; } }\n"
            )

            self.assertEqual(
                read_foam_entry(path, "tolerance", scope=["solvers", "V"]), "1e-5"
            )
            self.assertEqual(
                read_foam_entry(path, "tolerance", scope=["solvers", "p"]), "1e-7"
            )

    def test_writes_into_a_scope_written_as_an_inline_block(self) -> None:
        # The read path handles inline blocks; the write path indexes real
        # lines, so it must splice the new entry back into the original line
        # rather than reformat the file. Forced off foamDictionary, which
        # re-serialises the whole file and so cannot preserve the layout --
        # minimal-diff writing is a property of the Python writer alone.
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "fvSolution"
            path.write_text(
                "solvers { V { tolerance 1e-5; } p { tolerance 1e-7; } }\n"
            )

            with mock.patch(
                "openfoam_driver.core.runtime.mutators.shutil.which",
                return_value=None,
            ):
                update_foam_entry(
                    path, "tolerance", "1e-9", scope=["solvers", "V"]
                )

            self.assertEqual(
                read_foam_entry(path, "tolerance", scope=["solvers", "V"]), "1e-9"
            )
            # the sibling block must be untouched, and the file must still be
            # one line -- no wholesale reformat
            self.assertEqual(
                read_foam_entry(path, "tolerance", scope=["solvers", "p"]), "1e-7"
            )
            self.assertEqual(len(path.read_text().splitlines()), 1)

    def test_does_not_evaluate_dictionary_directives(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "testDict"
            path.write_text('a #calc "3.0 * 7.0";\n')

            # The literal entry, never the evaluated 21.
            self.assertEqual(read_foam_entry(path, "a"), '#calc "3.0 * 7.0"')


class TestUpdateFoamEntryPrefersFoamDictionary(unittest.TestCase):
    """update_foam_entry is the one sibling of read_foam_entry/remove_foam_dict/
    ensure_foam_dict that skipped the shutil.which("foamDictionary")
    preference -- these tests pin down that it now matches its siblings."""

    def test_prefers_foamdictionary_when_available(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "controlDict"
            path.write_text("deltaT 1e-06;\n")

            with mock.patch(
                "openfoam_driver.core.runtime.mutators.shutil.which",
                return_value="/usr/bin/foamDictionary",
            ), mock.patch(
                "openfoam_driver.core.runtime.mutators.update_foam_entry_via_foamDictionary"
            ) as mock_via_foamdictionary:
                update_foam_entry(path, "deltaT", 0.0001)

            mock_via_foamdictionary.assert_called_once_with(path, "deltaT", 0.0001, scope=None)

    def test_falls_back_to_regex_when_foamdictionary_raises(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "controlDict"
            path.write_text("deltaT 1e-06;\n")

            with mock.patch(
                "openfoam_driver.core.runtime.mutators.shutil.which",
                return_value="/usr/bin/foamDictionary",
            ), mock.patch(
                "openfoam_driver.core.runtime.mutators.update_foam_entry_via_foamDictionary",
                side_effect=RuntimeError("boom"),
            ):
                update_foam_entry(path, "deltaT", 0.0001)

            self.assertIn("deltaT    0.0001;", path.read_text())

    def test_uses_regex_directly_when_foamdictionary_unavailable(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "controlDict"
            path.write_text("deltaT 1e-06;\n")

            with mock.patch(
                "openfoam_driver.core.runtime.mutators.shutil.which",
                return_value=None,
            ):
                update_foam_entry(path, "deltaT", 0.0001)

            self.assertIn("deltaT    0.0001;", path.read_text())


class TestFoamDictionarySilentTruncationGuard(unittest.TestCase):
    """A malformed OpenFOAM header comment (missing the closing ``\\*---*/``
    line) makes ``foamDictionary`` treat the whole file as an empty dict; it
    then exits 0 after silently rewriting the file with only the newly-set
    key. update_foam_entry_via_foamDictionary must detect that and refuse to
    leave the file gutted."""

    @unittest.skipUnless(shutil.which("foamDictionary"), "foamDictionary not available")
    def test_detects_and_reverts_silent_truncation(self) -> None:
        text = "\n".join(
            [
                "/*--------------------------------*- C++ -*----------------------------------*\\",
                "FoamFile",
                "{",
                "    version     2.0;",
                "    format      ascii;",
                "    class       dictionary;",
                "    object      controlDict;",
                "}",
                "",
                "application     cardiacFoam;",
                "startFrom       startTime;",
                "startTime       0;",
                "stopAt          endTime;",
                "endTime         1.0;",
                "deltaT          1e-06;",
                "writeControl    adjustableRunTime;",
                "writeInterval   0.01;",
                "purgeWrite      0;",
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "controlDict"
            path.write_text(text)

            with self.assertRaises(RuntimeError):
                update_foam_entry_via_foamDictionary(path, "deltaT", 0.0001)

            reverted = path.read_text()
            self.assertIn("application", reverted)
            self.assertIn("writeInterval", reverted)


if __name__ == "__main__":
    unittest.main()
