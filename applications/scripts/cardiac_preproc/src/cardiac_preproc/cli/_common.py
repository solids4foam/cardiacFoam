"""Shared helpers for core CLI modules."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

PROJECT_ROOT = Path(__file__).resolve().parents[3]
CONFIGS_DIR = PROJECT_ROOT / "configs"


def load_config_module(filename: str, module_name: str | None = None) -> ModuleType:
    path = CONFIGS_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    name = module_name or path.stem
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load config from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
