"""Discover past runs under a directory tree.

`list_runs(root)` walks `root` recursively and yields one parsed
manifest per `run_manifest.json` found. Malformed manifests are
silently skipped so an unfinished or partially-written run does not
break agent recovery workflows.

Each yielded entry is the raw manifest dict augmented with a
``_manifest_path`` key carrying the absolute path to the source file —
agents use it to locate sibling sidecars (artifacts_manifest.json,
artifacts_realized.json, action_events.jsonl, run_report.md).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterator


def list_runs(root: Path) -> Iterator[dict]:
    """Yield every parseable `run_manifest.json` under `root`.

    Walks recursively. Order of iteration follows ``Path.rglob`` —
    filesystem-defined and not deterministic across platforms. Callers
    that need a stable order should sort by ``_manifest_path`` or
    ``started_at_utc``.
    """
    root = Path(root)
    if not root.is_dir():
        return
    for manifest_path in root.rglob("run_manifest.json"):
        try:
            payload = json.loads(manifest_path.read_text())
        except (json.JSONDecodeError, OSError):
            continue
        if not isinstance(payload, dict):
            continue
        payload["_manifest_path"] = str(manifest_path)
        yield payload
