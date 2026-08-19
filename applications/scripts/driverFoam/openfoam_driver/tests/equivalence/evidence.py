"""L1 evidence: what each condition started from, and where it ran.

E1 applies no dictionary overrides, so all three conditions must start from a
byte-identical staged tree. These helpers make that checkable, and record the
environment identity that could otherwise explain a numerical difference.
"""
from __future__ import annotations

import hashlib
import os
import shutil
from collections.abc import Mapping
from pathlib import Path


# Paths excluded from the input-tree hash: outputs and run state, not inputs.
# Open Question 1 in the spec -- extend this empirically as cases are added.
VOLATILE_PATHS: tuple[str, ...] = (
    "postProcessing",
    "workflow_logs",
    "workflow_state.json",
    "workflow_state.json",
    "logs",
)

# Environment variables that can plausibly change a numerical result.
_TRACKED_ENV: tuple[str, ...] = (
    "WM_PROJECT_DIR",
    "WM_PROJECT_VERSION",
    "WM_OPTIONS",
    "WM_COMPILER",
    "WM_PRECISION_OPTION",
    "WM_LABEL_SIZE",
    "FOAM_USER_LIBBIN",
)


def _is_volatile(relpath: str) -> bool:
    head = relpath.split("/", 1)[0]
    return head in VOLATILE_PATHS


def _sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hash_input_tree(case_root: Path) -> dict[str, str]:
    """Map every non-volatile file under `case_root` to its sha256.

    Keys are POSIX-style paths relative to `case_root`, so the result is
    comparable across machines and staging directories.
    """
    hashes: dict[str, str] = {}
    for path in sorted(case_root.rglob("*")):
        if not path.is_file():
            continue
        relpath = path.relative_to(case_root).as_posix()
        if _is_volatile(relpath):
            continue
        hashes[relpath] = _sha256_of(path)
    return hashes


def tree_digest(file_hashes: Mapping[str, str]) -> str:
    """Collapse a per-file hash map into one digest over sorted entries."""
    digest = hashlib.sha256()
    for relpath in sorted(file_hashes):
        digest.update(relpath.encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_hashes[relpath].encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def environment_manifest(env: Mapping[str, str] | None = None) -> dict[str, str | None]:
    """Record the environment identity a numerical difference could hide in.

    Absent variables are recorded as None rather than omitted, so two manifests
    always have the same key set and compare directly.
    """
    source = os.environ if env is None else env
    manifest: dict[str, str | None] = {
        name: source.get(name) for name in _TRACKED_ENV
    }
    solver = shutil.which("cardiacFoam", path=source.get("PATH"))
    manifest["solver_path"] = solver
    manifest["solver_sha256"] = _sha256_of(Path(solver)) if solver else None
    return manifest
