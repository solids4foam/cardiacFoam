from __future__ import annotations

from pathlib import Path

import pytest


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


@pytest.fixture
def tutorial_tree(tmp_path: Path) -> Path:
    """A miniature tutorials/ tree with one ran ECG case and one un-run case."""
    root = tmp_path / "tutorials"

    # Case A: a "ran" pseudo-ECG case under manufacturedSolutions.
    a = root / "manufacturedSolutions" / "monodomainPseudoECG"
    _write(a / "system" / "controlDict", "application myocardiumSolver;\n")
    _write(a / "constant" / "electroProperties",
           "myocardiumSolver monodomainSolver;\nionicModel TNNP;\n")
    _write(a / "constant" / "physicsProperties", "physicsModel electrophysiology;\n")
    _write(a / "README.md",
           "# Manufactured Pseudo ECG\n\nVerifies the monodomain pseudo-ECG on a "
           "manufactured domain.\n\n## Stack\n")
    _write(a / "postProcessing" / "pseudoECG.dat",
           "# time  V1  V2\n"
           "0.0  1.0  2.0\n"
           "1.0  3.0  -4.0\n")
    _write(a / "0" / "Vm", "0\n")
    _write(a / "1" / "Vm", "1\n")  # a numeric time dir > 0 => "ran"

    # Case B: an un-run 3D case under PATHOS (allowlisted for 3D).
    b = root / "PATHOS" / "RBBB"
    _write(b / "system" / "controlDict", "application myocardiumSolver;\n")
    _write(b / "constant" / "electroProperties",
           "myocardiumSolver monodomainSolver;\nionicModel BuenoOrovio;\n")
    _write(b / "README.md", "# RBBB\n\nRight bundle branch block tissue model.\n")

    return root


@pytest.fixture
def store_db(tmp_path: Path) -> Path:
    return tmp_path / "store.db"
