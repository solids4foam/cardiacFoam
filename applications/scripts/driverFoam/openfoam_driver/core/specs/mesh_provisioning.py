#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     mesh_provisioning
#
# Description
#     Renders the generic default blockMeshDict for cases built from scratch,
#     and the dx-to-cell-count arithmetic it shares with tutorial meshes.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Generic default geometry for a case built from scratch.

A case materialized purely from selectors/overrides has no geometry concept
at all, yet OpenFOAM needs a real mesh before any solver can run. This module
renders a sane generic default `system/blockMeshDict` (a small slab, "walls"
patch) that `blockMesh` turns into a mesh at run time. It is not a
scientifically tuned geometry for any specific tutorial -- callers that need
particular dimensions should still author their own blockMeshDict.

Which solver actually wants this default (versus a static polyMesh, or a real
anatomical mesh) is a solver-vocabulary decision and belongs to the active
plugin; cardiacFoam's lives in
`plugins/cardiacfoam/mesh_provisioning.py::provision_mesh`.

The default `blockMeshDict` is generated fresh from our own template each
time (like the plugin's own `system_templates.py`) -- there is no
pre-existing author file to parse or risk corrupting, so this does not need
`mutators.py`'s dictionary mutation machinery (that's for patching
values into an *already-written* file). `dx` (metres, isotropic cell size)
derives the cell count via `cell_counts_from_dx`, a small pure function
factored out so `niederer_2012.py`'s own `_replace_blockmesh_resolution`
(which *does* patch an existing author-provided file, a genuinely different
problem) can share the exact same divide-or-error math instead of duplicating
it.
"""

from __future__ import annotations

from collections.abc import Sequence

_DEFAULT_BLOCK_MESH_DICT = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v1912                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Generic default geometry for a from-scratch case_folder with no
// author-supplied blockMeshDict. Not tuned to any specific tutorial's
// science -- author your own blockMeshDict if geometry matters.

scale   0.001;

vertices
(
    (0 0 0)
    (2 0 0)
    (2 2 0)
    (0 2 0)
    (0 0 2)
    (2 0 2)
    (2 2 2)
    (0 2 2)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (__CELLS__ __CELLS__ __CELLS__) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    walls
    {
        type patch;
        faces
        (
            (3 7 6 2)
            (0 4 7 3)
            (2 6 5 1)
            (1 5 4 0)
            (0 3 2 1)
            (4 5 6 7)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
"""


# Fixed physical size of the generic default slab, in metres (matches the
# `scale 0.001` + 0..2 vertex extent above: 2 * 0.001 = 0.002 m). `dx` sweeps
# this domain's resolution by deriving a cell count from it -- there is no
# independent geometry concept here, just this one slab's coarseness.
_DEFAULT_SLAB_SIZE_M: tuple[float, float, float] = (0.002, 0.002, 0.002)
_DEFAULT_CELLS = 4


def cell_counts_from_dx(dx: float, slab_size: Sequence[float]) -> tuple[int, ...]:
    """Exact integer cell count per axis for isotropic cell size `dx` over
    each of `slab_size`'s axis lengths (same unit for both -- this function
    does no conversion).

    Raises `ValueError` if `dx` does not evenly divide every axis length:
    deliberately no rounding. A requested `dx` that doesn't fit the domain is
    a caller error to surface, not something to approximate quietly (mirrors
    the rigor `niederer_2012.py`'s own `_replace_blockmesh_resolution`
    already established for the same problem, on a different domain).
    """
    if dx <= 0:
        raise ValueError(f"dx must be positive; got {dx}")
    counts: list[int] = []
    for axis_length in slab_size:
        raw_cells = float(axis_length) / dx
        rounded_cells = round(raw_cells)
        if abs(raw_cells - rounded_cells) > 1e-9:
            raise ValueError(
                f"dx={dx} does not evenly divide slab axis length {axis_length}"
            )
        if rounded_cells <= 0:
            raise ValueError(
                f"Computed non-positive cell count for axis length {axis_length} with dx={dx}"
            )
        counts.append(int(rounded_cells))
    return tuple(counts)


def default_block_mesh_dict_text(*, dx_m: float | None = None) -> str:
    """Generic default `system/blockMeshDict` text (a small slab, "walls" patch).

    `dx_m` (metres, isotropic cell size) derives the cell count for the
    fixed `_DEFAULT_SLAB_SIZE_M` domain via `cell_counts_from_dx` -- since
    that domain is a cube, all three axes always get the same count. Omit
    `dx_m` for the fixed default cell count.
    """
    if dx_m is None:
        cells = _DEFAULT_CELLS
    else:
        counts = cell_counts_from_dx(dx_m, _DEFAULT_SLAB_SIZE_M)
        cells = counts[0]
    return _DEFAULT_BLOCK_MESH_DICT.replace("__CELLS__", str(cells))
