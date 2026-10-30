"""Release-level tests: the LaTeX notation, the figure scripts and the
self-contained release folder."""

import filecmp
import glob
import os
import re
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))

GENERATORS = ["fda_monodomain_mms.py", "rotated_monodomain_mms.py", "shell_monodomain_mms.py",
              "pseudo_ecg_mms.py", "graph_monodomain_mms.py", "coupled_1d3d_mms.py",
              "cartesian_eikonal_mms.py", "eikonal_ecg_mms.py", "shell_eikonal_mms.py",
              "template_mms.py"]
DOCUMENTS = ["docs/monodomain/monodomain_mms.qmd", "docs/eikonal/eikonal_mms.qmd"]

# Code spellings that must not reach the printed formulas (the paper writes
# V_m, psi, sigma, C_m and D, r for the template duration and rise).
FORBIDDEN = [r"\bVm\b", r"\btau\b", r"\bcm\b", r"\brise\b", r"\bduration\b"]


def _without_boldsymbol(text):
    return re.sub(r"\\boldsymbol\s*\\?\{?[A-Za-z]+\}?", "", text)


def test_no_code_symbols_in_latex():
    fragments = sorted(glob.glob(os.path.join(HERE, "latex", "*.tex")))
    assert fragments
    bad = []
    for path in fragments:
        with open(path) as fh:
            text = fh.read()
        for pattern in FORBIDDEN:
            if re.search(pattern, text):
                bad.append((os.path.basename(path), pattern))
        if re.search(r"(?<![\\A-Za-z])G\b", _without_boldsymbol(text)):
            bad.append((os.path.basename(path), "G"))
    assert not bad, bad


INCLUDE = re.compile(r"\{\{< include (\S+) >\}\}")


def _included_files(doc):
    """Every file reached from a document through include shortcodes. Quarto
    resolves a nested include against the directory of the entry file."""
    root = os.path.dirname(os.path.join(HERE, doc))
    found, stack = [], [os.path.join(HERE, doc)]
    while stack:
        path = stack.pop()
        if not path.endswith(".qmd") or not os.path.exists(path):
            continue
        with open(path) as fh:
            for rel in INCLUDE.findall(fh.read()):
                target = os.path.normpath(os.path.join(root, rel))
                if target not in found:
                    found.append(target)
                    stack.append(target)
    return found


def test_documents_include_existing_fragments():
    for doc in DOCUMENTS:
        missing = [p for p in _included_files(doc) if not os.path.exists(p)]
        assert not missing, (doc, missing)


def test_every_fragment_is_included():
    included = set()
    for doc in DOCUMENTS:
        included.update(os.path.basename(p) for p in _included_files(doc)
                        if os.path.basename(os.path.dirname(p)) == "latex")
    emitted = {os.path.basename(p) for p in glob.glob(os.path.join(HERE, "latex", "*.tex"))}
    assert not emitted - included, sorted(emitted - included)


def _copy_release(dest):
    shutil.copytree(HERE, dest, ignore=shutil.ignore_patterns(
        "parity", ".venv", "__pycache__", ".pytest_cache", "*_files", ".quarto"))


@pytest.mark.slow
def test_release_is_self_contained(tmp_path):
    root = tmp_path / "generators"
    _copy_release(root)
    venv = tmp_path / "venv"
    subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
    python = str(venv / "bin" / "python")
    subprocess.run([python, "-m", "pip", "install", "-q", "-r", "requirements.txt"],
                   cwd=root, check=True)
    for gen in GENERATORS:
        for flag in ("--emit", "--self-check"):
            r = subprocess.run([python, gen, flag], cwd=root, capture_output=True, text=True)
            assert r.returncode == 0, (gen, flag, r.stdout, r.stderr)
    emitted = glob.glob(os.path.join(HERE, "latex", "*.tex")) + \
        glob.glob(os.path.join(HERE, "reference", "*.py"))
    for path in emitted:
        rel = os.path.relpath(path, HERE)
        assert filecmp.cmp(path, root / rel, shallow=False), "%s differs after --emit" % rel
    for doc in DOCUMENTS:
        pdf = root / doc.replace(".qmd", ".pdf")
        pdf.unlink()
        r = subprocess.run(["quarto", "render", doc, "--to", "pdf"], cwd=root,
                           capture_output=True, text=True)
        assert r.returncode == 0, (doc, r.stderr[-2000:])
        assert pdf.exists() and pdf.stat().st_size > 100000, doc

ANALYTIC_FIGURES = {
    "monodomain": ["fda_cube.png", "rotated_cube.png", "shell_monodomain.png",
                   "pseudo_ecg.png", "graph_1d3d_geometry.png", "graph_1d3d_fields.png"],
    "eikonal": ["cartesian_eikonal.png", "eikonal_ecg.png", "shell_eikonal.png",
                "shell_eikonal_c0.png"],
}


def _run(mode, outdir, *extra):
    script = os.path.join(HERE, "docs", mode, "figures.py")
    return subprocess.run([sys.executable, script, "--draft", "--outdir", str(outdir),
                           *extra], capture_output=True, text=True, cwd=HERE)


@pytest.mark.parametrize("mode", sorted(ANALYTIC_FIGURES))
def test_figures_build_without_mesh_cases(mode, tmp_path):
    result = _run(mode, tmp_path)
    assert result.returncode == 0, result.stderr
    for name in ANALYTIC_FIGURES[mode]:
        path = tmp_path / name
        assert path.exists() and path.stat().st_size > 5000, name
    assert "--mesh-case" in result.stdout


def test_mesh_case_without_polymesh_fails_loudly(tmp_path):
    empty = tmp_path / "case"
    empty.mkdir()
    result = _run("monodomain", tmp_path / "out", "--mesh-case", "fda_cube_hex=%s" % empty)
    assert result.returncode == 1
    assert "fda_cube_hex" in result.stderr and "polyMesh" in result.stderr
