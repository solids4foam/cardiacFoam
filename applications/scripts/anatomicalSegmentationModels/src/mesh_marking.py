"""Mesh topology helpers and interface component extraction."""

from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np


def iter_cell_point_ids(mesh) -> Iterable[np.ndarray]:
    """Yield per-cell point-id arrays for common PyVista cell layouts."""
    if hasattr(mesh, "offset") and hasattr(mesh, "cell_connectivity"):
        offsets = np.asarray(mesh.offset, dtype=np.int64).reshape(-1)
        conn = np.asarray(mesh.cell_connectivity, dtype=np.int64).reshape(-1)
        if offsets.size >= 2:
            for i in range(offsets.size - 1):
                a = int(offsets[i])
                b = int(offsets[i + 1])
                if b > a:
                    yield conn[a:b]
            return

    if hasattr(mesh, "cells"):
        cells = np.asarray(mesh.cells, dtype=np.int64).reshape(-1)
        i = 0
        n = cells.size
        while i < n:
            npts = int(cells[i])
            i += 1
            if npts <= 0 or (i + npts) > n:
                break
            yield cells[i : i + npts]
            i += npts
        return

    raise ValueError("Unsupported mesh cell layout: cannot iterate cell point ids.")


class UnionFind:
    def __init__(self, items: Iterable[int]) -> None:
        self.parent: Dict[int, int] = {i: i for i in items}
        self.rank: Dict[int, int] = {i: 0 for i in items}

    def find(self, x: int) -> int:
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def extract_intraventricular_interface_components(
    mesh,
    intraventricular: np.ndarray,
    longitudinal_z: np.ndarray | None = None,
    z_min_exclusive: float | None = None,
    z_cap_max: float | None = None,
    lv_tag: float = -1.0,
    rv_tag: float = 1.0,
    tag_atol: float = 1.0e-8,
    min_component_points: int = 10,
) -> List[np.ndarray]:
    """Return connected LV/RV interface components (largest first)."""
    intrav = np.asarray(intraventricular).reshape(-1)
    z_mask = None
    if longitudinal_z is not None and z_cap_max is not None:
        z_vals = np.asarray(longitudinal_z).reshape(-1)
        if z_vals.shape[0] != intrav.shape[0]:
            raise ValueError("longitudinal_z and intraventricular must have same length.")
        z_mask = z_vals <= float(z_cap_max)
        if z_min_exclusive is not None:
            z_mask = z_mask & (z_vals > float(z_min_exclusive))

    interface_cells: List[np.ndarray] = []
    interface_point_ids: set[int] = set()
    for cell_ids in iter_cell_point_ids(mesh):
        ids = np.asarray(cell_ids, dtype=np.int64).reshape(-1)
        if z_mask is not None:
            ids = ids[z_mask[ids]]
        if ids.size == 0:
            continue
        vals = intrav[ids]
        has_lv = bool(np.any(np.isclose(vals, lv_tag, atol=tag_atol)))
        has_rv = bool(np.any(np.isclose(vals, rv_tag, atol=tag_atol)))
        if has_lv and has_rv:
            interface_cells.append(ids)
            interface_point_ids.update(int(i) for i in ids.tolist())

    if not interface_point_ids:
        return []

    uf = UnionFind(interface_point_ids)
    for ids in interface_cells:
        if ids.size <= 1:
            continue
        root_id = int(ids[0])
        for pid in ids[1:]:
            uf.union(root_id, int(pid))

    comps: Dict[int, List[int]] = {}
    for pid in interface_point_ids:
        r = uf.find(pid)
        comps.setdefault(r, []).append(pid)

    comp_lists = [np.asarray(c, dtype=np.int64) for c in comps.values()]
    if min_component_points > 1:
        filtered = [c for c in comp_lists if c.size >= int(min_component_points)]
        if len(filtered) >= 2:
            comp_lists = filtered
    comp_lists = sorted(comp_lists, key=lambda a: int(a.size), reverse=True)
    return comp_lists
