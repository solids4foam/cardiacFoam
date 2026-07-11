#!/usr/bin/env python3
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
#     distort_mesh
#
# Description
#     Applies a smooth, boundary-vanishing sinusoidal displacement to an
#     already-generated blockMesh points file on the unit [0,1]^d box, to
#     produce a non-orthogonal mesh for MMS convergence testing while
#     preserving the exact domain shape and boundary flux compatibility.
#
#     The mesh must be ASCII-formatted first:
#         foamFormatConvert -case <caseDir>
#         python3 distort_mesh.py <caseDir> -N 40 -A 0.15
#         checkMesh -case <caseDir>
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np

POINTS_BLOCK_RE = re.compile(
    r"\n(?P<count>\d+)\n\(\n(?P<body>(?:\([^()\n]*\)\n)+)\)\n"
)
POINT_RE = re.compile(r"\(([^()]*)\)")


def read_points(points_path: Path) -> tuple[str, str, np.ndarray]:
    text = points_path.read_text()

    if re.search(r"\bformat\s+binary\s*;", text):
        raise SystemExit(
            f"{points_path} is binary-formatted. Run "
            f"'foamFormatConvert -case <caseDir>' first."
        )

    match = POINTS_BLOCK_RE.search(text)
    if match is None:
        raise SystemExit(f"could not locate points block in {points_path}")

    count = int(match.group("count"))
    rows = POINT_RE.findall(match.group("body"))
    if len(rows) != count:
        raise SystemExit(
            f"points header declares {count} points but parsed {len(rows)}"
        )

    coords = np.array(
        [[float(v) for v in row.split()] for row in rows], dtype=np.float64
    )
    return text[: match.start()], text[match.end():], coords


def write_points(
    points_path: Path, preamble: str, trailer: str, coords: np.ndarray
) -> None:
    count = coords.shape[0]
    body = "\n".join(f"({x:.17g} {y:.17g} {z:.17g})" for x, y, z in coords)
    points_path.write_text(f"{preamble}\n{count}\n(\n{body}\n)\n{trailer}")


def distort(
    coords: np.ndarray,
    amplitude: float,
    h: float,
    dimension: int,
    tol: float,
) -> np.ndarray:
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    pi = np.pi

    if dimension == 3:
        dx = amplitude * h * np.sin(pi * x) * np.sin(pi * y) * np.sin(pi * z)
        dy = amplitude * h * np.sin(2 * pi * x) * np.sin(pi * y) * np.sin(pi * z)
        dz = amplitude * h * np.sin(pi * x) * np.sin(2 * pi * y) * np.sin(pi * z)
    elif dimension == 2:
        dx = amplitude * h * np.sin(pi * x) * np.sin(pi * y)
        dy = amplitude * h * np.sin(2 * pi * x) * np.sin(pi * y)
        dz = np.zeros_like(z)
    else:
        raise SystemExit("distortion is only meaningful for dimension in {2, 3}")

    on_boundary = np.zeros(coords.shape[0], dtype=bool)
    for comp in (x, y, z):
        on_boundary |= np.isclose(comp, 0.0, atol=tol) | np.isclose(comp, 1.0, atol=tol)

    dx[on_boundary] = 0.0
    dy[on_boundary] = 0.0
    dz[on_boundary] = 0.0

    displaced = coords.copy()
    displaced[:, 0] += dx
    displaced[:, 1] += dy
    displaced[:, 2] += dz
    return displaced


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument(
        "-N", "--resolution", type=int, required=True,
        help="cells per direction for this mesh (10/20/40/80/...); "
             "sets h=1/N so distortion amplitude scales with local cell size",
    )
    parser.add_argument(
        "-A", "--amplitude", type=float, default=0.15,
        help="dimensionless distortion severity; 0 reproduces the "
             "untouched orthogonal mesh exactly (default: 0.15)",
    )
    parser.add_argument("--dimension", type=int, choices=(2, 3), default=3)
    parser.add_argument(
        "--tol", type=float, default=1e-9,
        help="boundary-detection tolerance for masking displacement to zero "
             "on the six domain faces",
    )
    args = parser.parse_args()

    points_path = args.case_dir / "constant" / "polyMesh" / "points"
    preamble, trailer, coords = read_points(points_path)

    h = 1.0 / args.resolution
    displaced = distort(coords, args.amplitude, h, args.dimension, args.tol)

    delta = displaced - coords
    max_disp = float(np.max(np.linalg.norm(delta, axis=1)))
    n_moved = int(np.count_nonzero(np.any(delta != 0.0, axis=1)))

    print(
        f"{points_path}: {n_moved}/{coords.shape[0]} points displaced, "
        f"max |delta| = {max_disp:.6e} (h = {h:.6e}, A = {args.amplitude}, "
        f"dimension = {args.dimension}D)"
    )

    write_points(points_path, preamble, trailer, displaced)


if __name__ == "__main__":
    main()
