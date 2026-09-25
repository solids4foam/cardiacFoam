"""
Tests for applications/scripts/cardiacCompareCases.

Run with:  python3 -m unittest discover applications/scripts/tests
"""

import contextlib
import importlib.machinery
import importlib.util
import io
import os
import struct
import tempfile
import unittest

SCRIPT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "cardiacCompareCases"
)
loader = importlib.machinery.SourceFileLoader("cardiacCompareCases", SCRIPT)
spec = importlib.util.spec_from_loader(loader.name, loader)
ccc = importlib.util.module_from_spec(spec)
loader.exec_module(ccc)


HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Version:  {version}
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      {fmt};
    {arch}class       {cls};
    location    "0.1";
    object      Vm;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""


def asciiField(values, version="v2512"):
    body = "\n".join(repr(v) for v in values)
    return (
        HEADER.format(version=version, fmt="ascii", arch="", cls="volScalarField")
        + f"dimensions [0 2 -3 0 0 -1 0];\n\ninternalField nonuniform "
        f"List<scalar>\n{len(values)}\n(\n{body}\n)\n;\n\n"
        "boundaryField\n{\n    walls\n    {\n        type zeroGradient;\n"
        "    }\n}\n"
    ).encode()


def binaryField(values, vectors=()):
    head = HEADER.format(
        version="v2512", fmt="binary",
        arch='arch        "LSB;label=32;scalar=64";\n    ',
        cls="volScalarField",
    ).encode()
    data = head + b"dimensions [0 2 -3 0 0 -1 0];\n\ninternalField nonuniform "
    data += f"List<scalar> {len(values)}(".encode()
    data += struct.pack(f"<{len(values)}d", *values) + b")\n;\n\n"
    data += b"boundaryField\n{\n    inlet\n    {\n        type fixedValue;\n"
    flat = [c for v in vectors for c in v]
    data += f"        value nonuniform List<vector> {len(vectors)}(".encode()
    data += struct.pack(f"<{len(flat)}d", *flat) + b");\n    }\n}\n"
    return data


class CaseDirs:
    """Two temporary case directories populated from {relPath: bytes}."""

    def __init__(self, filesA, filesB):
        self.tmp = tempfile.TemporaryDirectory()
        self.a = os.path.join(self.tmp.name, "A")
        self.b = os.path.join(self.tmp.name, "B")
        for root, files in ((self.a, filesA), (self.b, filesB)):
            os.makedirs(root)
            for rel, data in files.items():
                path = os.path.join(root, rel)
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with open(path, "wb") as f:
                    f.write(data)

    def run(self, *options):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = ccc.main([*options, self.a, self.b])
        self.tmp.cleanup()
        return rc, out.getvalue()


class TestCompareCases(unittest.TestCase):

    def test_identical_cases_match(self):
        files = {"0.1/Vm": asciiField([1.0, 2.0]), "postProcessing/a.dat": b"# t v\n0 1\n"}
        rc, out = CaseDirs(files, files).run()
        self.assertEqual(rc, 0, out)
        self.assertIn("2 identical", out)

    def test_header_comment_differences_are_ignored(self):
        rc, out = CaseDirs(
            {"0.1/Vm": asciiField([1.0], version="v2412")},
            {"0.1/Vm": asciiField([1.0], version="v2512")},
        ).run()
        self.assertEqual(rc, 0, out)

    def test_number_formatting_is_ignored(self):
        rc, out = CaseDirs(
            {"a.dat": b"0 1e-05 2\n"}, {"a.dat": b"0 1.0e-5 2.000\n"}
        ).run()
        self.assertEqual(rc, 0, out)

    def test_value_difference_fails_by_default(self):
        rc, out = CaseDirs(
            {"0.1/Vm": asciiField([1.0, 2.0])},
            {"0.1/Vm": asciiField([1.0, 2.0 + 1e-13])},
        ).run()
        self.assertEqual(rc, 1, out)
        self.assertIn("DIFFER  0.1/Vm: 1/", out)

    def test_value_difference_within_tolerance_passes(self):
        rc, out = CaseDirs(
            {"0.1/Vm": asciiField([1.0, 2.0])},
            {"0.1/Vm": asciiField([1.0, 2.0 + 1e-13])},
        ).run("--rtol", "1e-12")
        self.assertEqual(rc, 0, out)
        self.assertIn("1 within tolerance", out)

    def test_atol_allows_differences_near_zero(self):
        files = lambda v: {"a.dat": f"0 {v}\n".encode()}
        self.assertEqual(CaseDirs(files(0.0), files(1e-20)).run("--rtol", "1e-12")[0], 1)
        self.assertEqual(CaseDirs(files(0.0), files(1e-20)).run("--atol", "1e-15")[0], 0)

    def test_structural_difference_fails(self):
        rc, out = CaseDirs(
            {"a.dat": b"# time Vm\n0 1\n"}, {"a.dat": b"# time phiE\n0 1\n"}
        ).run("--rtol", "1")
        self.assertEqual(rc, 1, out)
        self.assertIn("'Vm' vs 'phiE'", out)

    def test_different_lengths_fail(self):
        rc, out = CaseDirs({"a.dat": b"0 1\n1 2\n"}, {"a.dat": b"0 1\n"}).run()
        self.assertEqual(rc, 1, out)
        self.assertIn("different lengths", out)

    def test_missing_files_fail(self):
        rc, out = CaseDirs(
            {"0.1/Vm": b"1", "0.1/phiE": b"2"}, {"0.1/Vm": b"1", "0.1/Iion": b"3"}
        ).run()
        self.assertEqual(rc, 1, out)
        self.assertIn("ONLY A  0.1/phiE", out)
        self.assertIn("ONLY B  0.1/Iion", out)

    def test_logs_and_processors_excluded_by_default(self):
        rc, out = CaseDirs(
            {"log.cardiacFoam": b"ExecutionTime = 1 s", "processor0/0.1/Vm": b"1"},
            {"log.cardiacFoam": b"ExecutionTime = 2 s", "processor0/0.1/Vm": b"2"},
        ).run()
        self.assertEqual(rc, 0, out)

    def test_processors_option_includes_them(self):
        rc, out = CaseDirs(
            {"processor0/0.1/Vm": b"1"}, {"processor0/0.1/Vm": b"2"}
        ).run("--processors")
        self.assertEqual(rc, 1, out)

    def test_exclude_option(self):
        rc, out = CaseDirs(
            {"constant/electroProperties": b"a 1;", "a.dat": b"1"},
            {"constant/electroProperties": b"a 2;", "a.dat": b"1"},
        ).run("--exclude", "constant/*")
        self.assertEqual(rc, 0, out)

    def test_ignore_lines_option(self):
        rc, out = CaseDirs(
            {"a.dat": b"# Wall time 3.2\n0 1\n"}, {"a.dat": b"# Wall time 4.1\n0 1\n"}
        ).run("--ignore-lines", "Wall time")
        self.assertEqual(rc, 0, out)

    def test_nan_equals_nan(self):
        rc, out = CaseDirs({"a.dat": b"0 nan\n"}, {"a.dat": b"0 NaN\n"}).run()
        self.assertEqual(rc, 0, out)
        rc, out = CaseDirs({"a.dat": b"0 nan\n"}, {"a.dat": b"0 1\n"}).run("--rtol", "1")
        self.assertEqual(rc, 1, out)

    def test_binary_field_identical(self):
        field = binaryField([1.0, 2.0, 3.0], [(0.0, 1.0, 2.0)])
        rc, out = CaseDirs({"0.1/Vm": field}, {"0.1/Vm": field}).run()
        self.assertEqual(rc, 0, out)

    def test_binary_field_value_difference(self):
        a = binaryField([1.0, 2.0, 3.0], [(0.0, 1.0, 2.0)])
        b = binaryField([1.0, 2.0, 3.0], [(0.0, 1.0, 2.0 + 1e-10)])
        rc, out = CaseDirs({"0.1/Vm": a}, {"0.1/Vm": b}).run()
        self.assertEqual(rc, 1, out)
        self.assertIn("List<vector>", out)
        self.assertIn("item 2", out)
        rc, out = CaseDirs({"0.1/Vm": a}, {"0.1/Vm": b}).run("--rtol", "1e-9")
        self.assertEqual(rc, 0, out)

    def test_binary_field_size_difference(self):
        rc, out = CaseDirs(
            {"0.1/Vm": binaryField([1.0, 2.0])}, {"0.1/Vm": binaryField([1.0, 2.0, 3.0])}
        ).run("--rtol", "1")
        self.assertEqual(rc, 1, out)

    def test_other_binary_files_compared_bytewise(self):
        rc, out = CaseDirs(
            {"0.1/TNNPState": b"\0\1\2"}, {"0.1/TNNPState": b"\0\1\3"}
        ).run("--rtol", "1")
        self.assertEqual(rc, 1, out)
        self.assertIn("binary contents differ", out)


if __name__ == "__main__":
    unittest.main()
