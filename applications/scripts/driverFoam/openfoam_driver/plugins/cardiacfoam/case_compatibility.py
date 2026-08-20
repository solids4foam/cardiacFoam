"""Legacy cardiacFoam case-folder evidence behind a plugin-owned seam.

Why this exists
---------------
Older cardiacFoam tutorial folders may not provide ``Allrun``. driverFOAM
historically discovered them from ``electroProperties*`` and considered them
runnable when the matching physics and standard OpenFOAM system dictionaries
were present.

Activation
----------
The rule applies only while the cardiacFoam plugin is selected and core cannot
establish runnability from an executable ``Allrun`` first.

Compatibility
-------------
Filesystem-ingest, registry, and RunDocument execution tests preserve this
exact rule during Plan 1.  Plan 2 may replace it with broader user-facing
pipeline resolution; this module is the named boundary for that decision.
"""

from __future__ import annotations

from pathlib import Path


_REQUIRED_FILES = (
    "constant/physicsProperties",
    "system/controlDict",
    "system/fvSchemes",
    "system/fvSolution",
)


def has_case_marker(case_root: Path) -> bool:
    constant_root = case_root / "constant"
    if not constant_root.is_dir():
        return False
    return any(
        candidate.is_file() and candidate.name.startswith("electroProperties")
        for candidate in constant_root.rglob("electroProperties*")
    )


def is_runnable_without_workflow(case_root: Path) -> bool:
    return has_case_marker(case_root) and all(
        (case_root / relpath).exists() for relpath in _REQUIRED_FILES
    )
