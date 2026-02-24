#!/usr/bin/env python3
"""Compatibility runner that delegates CLI orchestration to the engine layer."""

from __future__ import annotations

from pathlib import Path

from cardiac_preproc.engine import run_pipeline_cli

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    run_pipeline_cli(project_root=PROJECT_ROOT)


if __name__ == "__main__":
    main()
