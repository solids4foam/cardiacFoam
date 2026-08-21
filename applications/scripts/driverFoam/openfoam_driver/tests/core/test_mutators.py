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
import tempfile
import unittest
from pathlib import Path

import pytest

from openfoam_driver.core.runtime.mutators import (
    ensure_foam_dict,
    read_foam_entry,
    remove_foam_dict,
    update_foam_entry,
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

    def test_falls_back_for_c_style_comments(self) -> None:
        """Tier 1's brace counting is comment-unaware (it only strips ``//``),
        so a ``/* ... { ... */`` block comment's literal brace throws off its
        depth tracking and it reports the enclosing scope as having
        unbalanced braces. That KeyError is exactly the class the foamlib
        tier exists to catch: a real parser understands ``/* */`` natively,
        so the update now succeeds instead of raising.
        """
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

            update_foam_entry(path, "value", 2, scope="someDict")

            updated = path.read_text()
            self.assertIn("value    2;", updated)
            self.assertIn("a brace { inside it", updated)

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

    def test_remove_foam_dict_missing_ok_tolerates_absent_scope(
        self,
    ) -> None:
        # _resolve_search_region(lines, scope) previously raised
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

    These tests pin the Python implementation directly, which is the side
    that was wrong.
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
        # rather than reformat the file -- minimal-diff writing is a property
        # of the Python writer alone.
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "fvSolution"
            path.write_text(
                "solvers { V { tolerance 1e-5; } p { tolerance 1e-7; } }\n"
            )

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


class TestUpdateFoamEntryUsesLineTier(unittest.TestCase):
    """update_foam_entry runs the line-based tier 1 unconditionally now --
    there is no foamDictionary preference or availability check left to pin."""

    def test_updates_an_existing_scalar_entry(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "controlDict"
            path.write_text("deltaT 1e-06;\n")

            update_foam_entry(path, "deltaT", 0.0001)

            self.assertIn("deltaT    0.0001;", path.read_text())


def test_update_foam_entry_falls_back_for_brace_in_quoted_string(tmp_path):
    """Tier 1 cannot parse this; the foamlib tier must pick it up.

    A brace inside a quoted value defeats brace counting, so the line scanner
    reports the enclosing scope as missing. Before the foamlib backend this
    raised KeyError whenever OpenFOAM was not sourced.
    """
    path = tmp_path / "d"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object d; }\n"
        'note  "a value with { an unbalanced brace";\n'
        "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n"
    )
    update_foam_entry(path, "tolerance", 1e-12, scope=["solvers", "Vm"])
    text = path.read_text()
    assert "1e-12" in text
    assert 'note  "a value with { an unbalanced brace";' in text


def test_remove_foam_dict_falls_back_for_brace_in_quoted_string(tmp_path):
    """Same defeats-the-scanner fixture as the update_foam_entry test above,

    but for remove_foam_dict's own early _resolve_search_region call. Before
    this fix, remove_foam_dict raised the same KeyError update_foam_entry
    used to raise here, with no fallback -- a real regression versus the old
    foamDictionary-first behaviour in a sourced environment.
    """
    path = tmp_path / "d"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object d; }\n"
        'note  "a value with { an unbalanced brace";\n'
        "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n"
    )
    remove_foam_dict(path, "Vm", scope=["solvers"])
    text = path.read_text()
    assert "Vm" not in text
    assert 'note  "a value with { an unbalanced brace";' in text


def test_remove_foam_dict_falls_back_when_name_matches_a_scalar_entry(tmp_path):
    """remove_foam_dict must handle a name that's a scalar, not a block.

    Reproduced directly: before this fix, remove_foam_dict("xMin", ...)
    raised KeyError("... has no opening brace") unconditionally -- not even
    honoring missing_ok=True -- whenever the matched name turned out to be a
    plain `xMin -1.0;` entry rather than a `{ ... }` block. This is exactly
    the case plugins/cardiacfoam/tutorials/manufactured_bath_bidomain.py
    worked around by calling the separate remove_foam_entry function
    instead. foamlib's `del` has no such shape restriction.
    """
    path = tmp_path / "d"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object d; }\n"
        "groundPatches\n{\n    xMin -1.0;\n    yMin 0.0;\n}\n"
    )
    remove_foam_dict(path, "xMin", scope=["groundPatches"])
    text = path.read_text()
    assert "xMin" not in text
    assert "yMin 0.0;" in text


def test_remove_foam_dict_missing_ok_still_uses_fallback(tmp_path):
    """missing_ok=True must not disable the foamlib fallback entirely.

    Before this fix, both of remove_foam_dict's delegation points checked
    `if missing_ok: return` before ever calling foam_backend -- so a caller
    that set missing_ok=True got a silent no-op even when the target
    genuinely existed in a file the line scanner couldn't parse. This
    fixture is deletable via foam_backend (the same brace-in-quoted-string
    case the other fallback tests use); it must actually be removed, not
    silently skipped, when missing_ok=True.
    """
    path = tmp_path / "d"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object d; }\n"
        'note  "a value with { an unbalanced brace";\n'
        "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n"
    )
    remove_foam_dict(path, "Vm", scope=["solvers"], missing_ok=True)
    text = path.read_text()
    assert "Vm" not in text
    assert 'note  "a value with { an unbalanced brace";' in text


def test_scope_resolution_does_not_treat_a_scalar_entry_as_a_block(tmp_path):
    """A scalar entry sharing a scope name must not be silently treated as

    a block by walking forward into an unrelated sibling's braces.
    Reproduced directly: scope=["outer", "Vm"] against a scalar "Vm 5;"
    line followed by an unrelated "unrelatedBlock { tolerance 1e-9; }"
    previously returned unrelatedBlock's tolerance as if it belonged to
    Vm's scope -- a silent wrong answer, not even a KeyError.
    """
    path = tmp_path / "d"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object d; }\n"
        "outer\n{\n"
        "    Vm 5;\n"
        "    unrelatedBlock\n    {\n        tolerance 1e-9;\n    }\n"
        "}\n"
    )
    assert read_foam_entry(path, "tolerance", scope=["outer", "Vm"]) is None


REGEX_KEYED = (
    "FoamFile { version 2.0; class dictionary; object fvSolution; }\n"
    "solvers\n"
    "{\n"
    '    "Vm|VmFinal|u|uFinal"\n'
    "    {\n"
    "        solver          PCG;\n"
    "        tolerance       1e-11;\n"
    "    }\n"
    '    "psi|psiFinal"\n'
    "    {\n"
    "        solver          PCG;\n"
    "        tolerance       1e-11;\n"
    "    }\n"
    "}\n"
)


def test_regex_key_resolves_member_name(tmp_path):
    path = tmp_path / "fvSolution"
    path.write_text(REGEX_KEYED)
    assert read_foam_entry(path, "tolerance", scope=["solvers", "Vm"]) == "1e-11"


def test_regex_key_write_targets_the_matching_block(tmp_path):
    path = tmp_path / "fvSolution"
    path.write_text(REGEX_KEYED)
    update_foam_entry(path, "tolerance", 1e-12, scope=["solvers", "psi"])
    text = path.read_text()
    assert text.count("1e-12") == 1
    assert text.index("1e-12") > text.index('"psi|psiFinal"')


def test_literal_key_wins_over_pattern(tmp_path):
    path = tmp_path / "fvSolution"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object fvSolution; }\n"
        "solvers\n{\n"
        '    "Vm|psi" { tolerance 1e-11; }\n'
        "    Vm       { tolerance 1e-09; }\n"
        "}\n"
    )
    assert read_foam_entry(path, "tolerance", scope=["solvers", "Vm"]) == "1e-09"


def test_last_declared_pattern_wins(tmp_path):
    path = tmp_path / "fvSolution"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object fvSolution; }\n"
        "solvers\n{\n"
        '    "V.*"  { tolerance 1e-11; }\n'
        '    "Vm.*" { tolerance 1e-09; }\n'
        "}\n"
    )
    assert read_foam_entry(path, "tolerance", scope=["solvers", "Vm"]) == "1e-09"


def test_name_matching_no_pattern_still_fails_closed(tmp_path):
    path = tmp_path / "fvSolution"
    path.write_text(REGEX_KEYED)
    assert read_foam_entry(path, "tolerance", scope=["solvers", "nope"]) is None


_SINGLE_PATTERN_BLOCK = (
    "FoamFile { version 2.0; class dictionary; object fvSolution; }\n"
    "solvers\n{\n"
    '    "Vm|VmFinal"\n    {\n        tolerance 1e-11;\n    }\n'
    "}\n"
)


def test_remove_foam_dict_resolves_a_pattern_keyed_member(tmp_path):
    """remove_foam_dict's own target-name matching must resolve patterns too.

    _find_dict_block_bounds (used for scope-PATH segments) already resolves
    quoted-regex headers. remove_foam_dict has its own separate scan for the
    dict being removed, which did not reuse that resolution -- so removing
    "Vm" from a block keyed only by "Vm|VmFinal" fell through to the
    foamlib fallback and raised KeyError, even though "Vm" plainly resolves
    against the pattern by OpenFOAM's own rules.
    """
    path = tmp_path / "fvSolution"
    path.write_text(_SINGLE_PATTERN_BLOCK)
    remove_foam_dict(path, "Vm", scope=["solvers"])
    text = path.read_text()
    assert "Vm|VmFinal" not in text
    assert "tolerance" not in text


def test_ensure_foam_dict_does_not_duplicate_a_pattern_covered_member(tmp_path):
    """ensure_foam_dict must recognize a name already covered by a pattern.

    Without pattern resolution, ensure_foam_dict("Vm", ...) on a block keyed
    only by "Vm|VmFinal" reports "Vm" as missing and inserts a duplicate
    literal Vm block -- which OpenFOAM's literal-beats-pattern precedence
    then silently prefers over the existing, intentionally-shared one.
    """
    path = tmp_path / "fvSolution"
    path.write_text(_SINGLE_PATTERN_BLOCK)
    inserted = ensure_foam_dict(
        path, "Vm", "Vm\n{\n    tolerance 1e-9;\n}\n", scope=["solvers"]
    )
    assert inserted is False
    text = path.read_text()
    assert text.count("tolerance") == 1
    assert "1e-11" in text


@pytest.mark.parametrize(
    "payload",
    [
        "1e-6;  rogue  1",
        '#calc "2*3"',
        "#codeStream { code #{ os << 1; #}; }",
        "0.3;\n    endTime 99",
    ],
)
def test_update_foam_entry_rejects_injected_value(tmp_path, payload):
    path = tmp_path / "controlDict"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object controlDict; }\n"
        "deltaT 1e-06;\n"
    )
    original = path.read_text()
    with pytest.raises(ValueError):
        update_foam_entry(path, "deltaT", payload)
    assert path.read_text() == original


def test_update_foam_entry_still_accepts_ordinary_values(tmp_path):
    path = tmp_path / "controlDict"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object controlDict; }\n"
        "deltaT 1e-06;\n"
    )
    update_foam_entry(path, "deltaT", "5e-06")
    assert "5e-06" in path.read_text()


if __name__ == "__main__":
    unittest.main()
