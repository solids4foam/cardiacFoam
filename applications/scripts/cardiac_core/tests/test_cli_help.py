"""Smoke tests for canonical CLI entrypoints."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run_help(script_relpath):
    result = subprocess.run(
        [sys.executable, str(ROOT / script_relpath), "--help"],
        capture_output=True,
        text=True,
    )
    return result.returncode


def test_diffusivity_help():
    assert run_help("diffusivity/diffusivity.py") == 0


def test_scar_help():
    assert run_help("scar/scar.py") == 0


def test_purkinje_slab_help():
    assert run_help("purkinje/slab/slab.py") == 0


def test_pipeline_help():
    assert run_help("conduction_system_generation.py") == 0
