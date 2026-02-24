"""Smoke tests: each CLI script must respond to --help without error."""
import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"


def run_help(script_name):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS_DIR / script_name), "--help"],
        capture_output=True,
        text=True,
    )
    return result.returncode


def test_diffusivity_vtk_help():
    assert run_help("diffusivity_vtk.py") == 0


def test_scar_vtk_help():
    assert run_help("scar_vtk.py") == 0


def test_purkinje_slab_help():
    assert run_help("purkinje_slab.py") == 0
