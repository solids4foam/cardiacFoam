"""Small meshes for testing the continuum junction (mms_common/pvj_continuum.py).

Point grids of a box or of a sector of a spherical shell, their hexahedra and the six-tetrahedra
(Kuhn) split of each hexahedron, and the polyhedral-mesh arrays (faces, owner, neighbour) of a
mesh of hexahedra or tetrahedra, in the form cells_from_polymesh takes and OpenFOAM stores.
"""

import numpy as np

HEX_FACES = ((0, 4, 7, 3), (1, 2, 6, 5), (0, 1, 5, 4), (3, 7, 6, 2), (0, 3, 2, 1), (4, 5, 6, 7))
TET_FACES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))


def structured_points(n, lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0)):
    """Grid (n[0]+1, n[1]+1, n[2]+1, 3) of points of the box [lo, hi]; n an int or 3 ints."""
    n = (n,)*3 if np.isscalar(n) else tuple(n)
    axes = [np.linspace(lo[i], hi[i], n[i] + 1) for i in range(3)]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)


def jitter_inner(grid, amount, seed=0):
    """Move the interior points of a grid by up to `amount` of the local grid spacing."""
    rng = np.random.default_rng(seed)
    n = np.array(grid.shape[:3]) - 1
    spacing = (grid[-1, -1, -1] - grid[0, 0, 0])/n
    inner = np.ones(grid.shape[:3], dtype=bool)
    inner[[0, -1], :, :] = inner[:, [0, -1], :] = inner[:, :, [0, -1]] = False
    return grid + amount*spacing*rng.uniform(-1, 1, grid.shape)*inner[..., None]


def hexahedra(grid):
    """(points (m, 3), cells (n, 8)) of the hexahedra of a point grid, in OpenFOAM order."""
    nx, ny, nz = (s - 1 for s in grid.shape[:3])
    index = np.arange(grid[..., 0].size).reshape(grid.shape[:3])
    i, j, k = (a.reshape(-1) for a in np.meshgrid(range(nx), range(ny), range(nz), indexing="ij"))
    cells = np.stack([index[i, j, k], index[i + 1, j, k], index[i + 1, j + 1, k],
                      index[i, j + 1, k], index[i, j, k + 1], index[i + 1, j, k + 1],
                      index[i + 1, j + 1, k + 1], index[i, j + 1, k + 1]], axis=1)
    return grid.reshape(-1, 3), cells


def kuhn_tetrahedra(grid):
    """(points, cells (6 n, 4)): each cube of a point grid split into the six tetrahedra of
    the Kuhn (Freudenthal) triangulation, which conform across cubes."""
    points, hexes = hexahedra(grid)
    # hexahedron vertex of binary offset (a, b, c): index into the OpenFOAM order
    vertex = {(0, 0, 0): 0, (1, 0, 0): 1, (1, 1, 0): 2, (0, 1, 0): 3,
              (0, 0, 1): 4, (1, 0, 1): 5, (1, 1, 1): 6, (0, 1, 1): 7}
    cells = []
    for axes in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
        corner = np.zeros(3, dtype=int)
        path = [tuple(corner)]
        for ax in axes:
            corner[ax] = 1
            path.append(tuple(corner))
        cells.append(hexes[:, [vertex[c] for c in path]])
    return points, np.concatenate(cells)


def shell_sector_grid(n, r1=0.5, r2=1.0, theta=(1.0, np.pi - 1.0), phi=(-0.7, 0.7), jitter=0.0,
                      seed=0):
    """Point grid of a sector of the spherical shell r1 <= |x| <= r2: the box (r, theta, phi)
    mapped to Cartesian coordinates, so the hexahedra are curved; the polar axis is z, so the
    default sector contains the points (r, 0, 0) of the x axis, and a terminal on the inner
    sphere at (r1, 0, 0) lies at the middle of its inner face."""
    grid = structured_points(n, (r1, theta[0], phi[0]), (r2, theta[1], phi[1]))
    if jitter:
        grid = jitter_inner(grid, jitter, seed)
    r, th, ph = grid[..., 0], grid[..., 1], grid[..., 2]
    return np.stack([r*np.sin(th)*np.cos(ph), r*np.sin(th)*np.sin(ph), r*np.cos(th)], axis=-1)


def _outward(points, cell, faces):
    """Faces of a cell as vertex-index tuples whose area vectors point out of the cell."""
    centre = points[list(cell)].mean(axis=0)
    out = []
    for face in faces:
        idx = [cell[i] for i in face]
        p = points[idx]
        area = 0.5*sum(np.cross(p[i], p[(i + 1) % len(p)]) for i in range(len(p)))
        out.append(tuple(idx) if area @ (p.mean(axis=0) - centre) >= 0 else tuple(idx[::-1]))
    return out


def polymesh_arrays(points, cells):
    """(faces, owner, neighbour, n_internal) of a mesh of hexahedra (cells (n, 8)) or
    tetrahedra (cells (n, 4)) in polyMesh order: internal faces first, by owner then
    neighbour, each stored once with its vertex order pointing out of the owner (lower
    numbered cell); then the boundary faces."""
    points = np.asarray(points, dtype=float)
    faces_of = HEX_FACES if cells.shape[1] == 8 else TET_FACES
    by_key = {}
    for c, cell in enumerate(cells.tolist()):
        for face in _outward(points, cell, faces_of):
            by_key.setdefault(tuple(sorted(face)), []).append((c, face))
    internal, boundary = [], []
    for entries in by_key.values():
        if len(entries) == 2:
            (c0, f0), (c1, _) = sorted(entries)
            internal.append((c0, c1, f0))
        elif len(entries) == 1:
            boundary.append((entries[0][0], 0, entries[0][1]))
        else:
            raise ValueError("a face is shared by %d cells" % len(entries))
    internal.sort(key=lambda e: (e[0], e[1]))
    boundary.sort(key=lambda e: (e[0], e[2]))
    ordered = internal + boundary
    return ([f for _, _, f in ordered], [c0 for c0, _, _ in ordered],
            [c1 for _, c1, _ in internal], len(internal))
