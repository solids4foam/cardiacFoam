"""Single source of truth for run-document.json (P2.8).

applications/scripts/driverFoam/schemas/run-document.json is hand-authored.
openfoam_driver/schemas/run-document.json (the packaged resource) is
generated from it byte-for-byte. Run this after any hand-edit to the source
copy; CI (or test_packaged_schema_resource_matches_fixture_schema) fails if
they drift.
"""
from __future__ import annotations

import shutil
from pathlib import Path

_SOURCE = Path(__file__).parent / "run-document.json"
_PACKAGED = (
    Path(__file__).parents[1] / "openfoam_driver" / "schemas" / "run-document.json"
)


def main() -> None:
    shutil.copyfile(_SOURCE, _PACKAGED)
    print(f"generated {_PACKAGED} from {_SOURCE}")


if __name__ == "__main__":
    main()
