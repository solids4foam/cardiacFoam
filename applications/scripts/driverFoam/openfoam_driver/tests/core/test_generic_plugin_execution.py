from __future__ import annotations

import os
import stat
from pathlib import Path

from openfoam_driver.cli import main


def test_generic_plugin_executes_plain_allrun_case(tmp_path: Path) -> None:
    case_root = tmp_path / "plainOpenFoamCase"
    case_root.mkdir()
    allrun = case_root / "Allrun"
    allrun.write_text("#!/bin/sh\nprintf complete > generic-proof.txt\n")
    allrun.chmod(allrun.stat().st_mode | stat.S_IXUSR)

    exit_code = main([
        "run",
        "--strict",
        "--plugin", "none",
        "--entry", "plainOpenFoamCase",
        "--tutorials-root", str(tmp_path),
    ])

    assert exit_code == 0
    assert (case_root / "generic-proof.txt").read_text() == "complete"
    assert (case_root / "postProcessing" / "workflow_state.json").is_file()
