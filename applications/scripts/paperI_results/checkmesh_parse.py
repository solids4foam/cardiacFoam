#----------------------------------------------------------------------------#
# Module
#     checkmesh_parse
#
# Description
#     Extracts mesh-validity and non-orthogonality metrics from a
#     checkMesh log, for the non-orthogonal MMS distortion sweep.
#----------------------------------------------------------------------------#

from __future__ import annotations

import re
from dataclasses import dataclass

_NON_ORTHO_RE = re.compile(
    r"non-orthogonality\s+Max:\s*([\d.eE+-]+)\s+average:\s*([\d.eE+-]+)",
    re.IGNORECASE,
)
_NEGATIVE_VOLUME_RE = re.compile(
    r"cells with negative volumes|has negative volume", re.IGNORECASE
)
_MESH_OK_RE = re.compile(r"^Mesh OK\.\s*$", re.MULTILINE)


@dataclass
class CheckMeshResult:
    mesh_ok: bool
    max_non_orthogonality: float | None
    average_non_orthogonality: float | None
    has_negative_volume: bool


def parse_checkmesh_log(text: str) -> CheckMeshResult:
    non_ortho = _NON_ORTHO_RE.search(text)
    return CheckMeshResult(
        mesh_ok=bool(_MESH_OK_RE.search(text)),
        max_non_orthogonality=float(non_ortho.group(1)) if non_ortho else None,
        average_non_orthogonality=float(non_ortho.group(2)) if non_ortho else None,
        has_negative_volume=bool(_NEGATIVE_VOLUME_RE.search(text)),
    )
