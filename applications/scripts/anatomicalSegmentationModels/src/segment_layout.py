"""Ring boundary construction and configurable theta->segment mapping."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from .calc_utils import TAU, circular_mean, cw_is_between, cw_unwrap, wrap_angle
from .groove_calculator import (
    groove_angle_from_point_mask,
    infer_groove_angles_debug_from_intraventricular_interface,
    infer_groove_angles_from_intraventricular_interface,
)
from .mesh_marking import extract_intraventricular_interface_components


@dataclass(frozen=True)
class GrooveAngles:
    """Anterior and posterior interventricular groove angles (radians)."""

    anterior: float
    posterior: float

    def normalized(self) -> "GrooveAngles":
        return GrooveAngles(float(wrap_angle(self.anterior)), float(wrap_angle(self.posterior)))


def groove_angles_from_point_mask(
    points_xyz: np.ndarray,
    rotational_theta: np.ndarray,
    mask: np.ndarray,
) -> float:
    """Compatibility wrapper (points_xyz kept for API compatibility)."""
    _ = points_xyz
    return groove_angle_from_point_mask(rotational_theta, mask)


def groove_angles_from_intraventricular_interface(
    mesh,
    rotational_theta: np.ndarray,
    intraventricular: np.ndarray,
    longitudinal_z: np.ndarray | None = None,
    z_min_exclusive: float | None = None,
    z_cap_max: float | None = None,
    rot_start: float = -np.pi,
    lv_tag: float = -1.0,
    rv_tag: float = 1.0,
    tag_atol: float = 1.0e-8,
    min_component_points: int = 10,
) -> GrooveAngles:
    ant, post = infer_groove_angles_from_intraventricular_interface(
        mesh=mesh,
        rotational_theta=rotational_theta,
        intraventricular=intraventricular,
        longitudinal_z=longitudinal_z,
        z_min_exclusive=z_min_exclusive,
        z_cap_max=z_cap_max,
        rot_start=rot_start,
        lv_tag=lv_tag,
        rv_tag=rv_tag,
        tag_atol=tag_atol,
        min_component_points=min_component_points,
    )
    return GrooveAngles(ant, post)


def groove_angles_debug_from_intraventricular_interface(
    mesh,
    rotational_theta: np.ndarray,
    intraventricular: np.ndarray,
    longitudinal_z: np.ndarray | None = None,
    z_min_exclusive: float | None = None,
    z_cap_max: float | None = None,
    rot_start: float = -np.pi,
    lv_tag: float = -1.0,
    rv_tag: float = 1.0,
    tag_atol: float = 1.0e-8,
    min_component_points: int = 10,
):
    ant, post, debug = infer_groove_angles_debug_from_intraventricular_interface(
        mesh=mesh,
        rotational_theta=rotational_theta,
        intraventricular=intraventricular,
        longitudinal_z=longitudinal_z,
        z_min_exclusive=z_min_exclusive,
        z_cap_max=z_cap_max,
        rot_start=rot_start,
        lv_tag=lv_tag,
        rv_tag=rv_tag,
        tag_atol=tag_atol,
        min_component_points=min_component_points,
    )
    return GrooveAngles(ant, post), debug


def interface_components_from_intraventricular(
    mesh,
    intraventricular: np.ndarray,
    longitudinal_z: np.ndarray | None = None,
    z_min_exclusive: float | None = None,
    z_cap_max: float | None = None,
    lv_tag: float = -1.0,
    rv_tag: float = 1.0,
    tag_atol: float = 1.0e-8,
    min_component_points: int = 10,
):
    return extract_intraventricular_interface_components(
        mesh=mesh,
        intraventricular=intraventricular,
        longitudinal_z=longitudinal_z,
        z_min_exclusive=z_min_exclusive,
        z_cap_max=z_cap_max,
        lv_tag=lv_tag,
        rv_tag=rv_tag,
        tag_atol=tag_atol,
        min_component_points=min_component_points,
    )


def build_anatomical_boundaries_6(
    grooves: GrooveAngles,
) -> Dict[str, Tuple[float, float]]:
    g = grooves.normalized()
    a, p = g.anterior, g.posterior
    sept_len = cw_unwrap(a, p)
    sep_mid = wrap_angle(a + 0.5 * sept_len)
    free_len = TAU - sept_len
    f1 = wrap_angle(p + 1.0 / 4.0 * free_len)
    f2 = wrap_angle(p + 2.0 / 4.0 * free_len)
    f3 = wrap_angle(p + 3.0 / 4.0 * free_len)
    return {
        "anteroseptal": (a, sep_mid),
        "inferoseptal": (sep_mid, p),
        "inferior": (p, f1),
        "inferolateral": (f1, f2),
        "anterolateral": (f2, f3),
        "anterior": (f3, a),
    }


def build_anatomical_boundaries_4(
    grooves: GrooveAngles,
) -> Dict[str, Tuple[float, float]]:
    g = grooves.normalized()
    a, p = g.anterior, g.posterior
    sept_len = cw_unwrap(a, p)
    free_len = TAU - sept_len
    f1 = wrap_angle(p + 1.0 / 3.0 * free_len)
    f2 = wrap_angle(p + 2.0 / 3.0 * free_len)
    return {
        "septal": (a, p),
        "inferior": (p, f1),
        "lateral": (f1, f2),
        "anterior": (f2, a),
    }


def segment_ids_from_theta(
    theta: np.ndarray,
    ring: str,
    grooves: GrooveAngles,
) -> np.ndarray:
    t = wrap_angle(theta)
    ring = ring.lower()
    seg = np.full(t.shape, -1, dtype=np.int32)
    if ring in ("basal", "mid"):
        b = build_anatomical_boundaries_6(grooves)
        if ring == "basal":
            name_to_id = {
                "anterior": 1,
                "anteroseptal": 2,
                "inferoseptal": 3,
                "inferior": 4,
                "inferolateral": 5,
                "anterolateral": 6,
            }
        else:
            name_to_id = {
                "anterior": 7,
                "anteroseptal": 8,
                "inferoseptal": 9,
                "inferior": 10,
                "inferolateral": 11,
                "anterolateral": 12,
            }
        for nm, (s, e) in b.items():
            seg[cw_is_between(t, s, e)] = name_to_id[nm]
        return seg
    if ring == "apical":
        b = build_anatomical_boundaries_4(grooves)
        name_to_id = {"anterior": 13, "septal": 14, "inferior": 15, "lateral": 16}
        for nm, (s, e) in b.items():
            seg[cw_is_between(t, s, e)] = name_to_id[nm]
        return seg
    raise ValueError(f"Unknown ring '{ring}'. Expected one of: basal, mid, apical.")


def build_ring_boundaries(
    grooves: GrooveAngles,
    septal_parts: int,
    freewall_parts: int,
) -> List[Tuple[float, float]]:
    if septal_parts <= 0 or freewall_parts <= 0:
        raise ValueError("septal_parts and freewall_parts must be >= 1")
    g = grooves.normalized()
    a, p = g.anterior, g.posterior
    sept_len = cw_unwrap(a, p)
    free_len = TAU - sept_len
    bounds: List[Tuple[float, float]] = []
    for i in range(septal_parts):
        s = wrap_angle(a + (i / septal_parts) * sept_len)
        e = wrap_angle(a + ((i + 1) / septal_parts) * sept_len)
        bounds.append((float(s), float(e)))
    for i in range(freewall_parts):
        s = wrap_angle(p + (i / freewall_parts) * free_len)
        e = wrap_angle(p + ((i + 1) / freewall_parts) * free_len)
        bounds.append((float(s), float(e)))
    return bounds


def segment_ids_from_theta_layout(
    theta: np.ndarray,
    ring: str,
    grooves: GrooveAngles,
    ring_layout: Dict[str, Tuple[int, int]],
    ring_id_map: Dict[str, List[int]],
) -> np.ndarray:
    t = wrap_angle(theta)
    r = ring.lower()
    if r not in ring_layout:
        raise ValueError(f"Ring '{r}' missing in ring_layout.")
    if r not in ring_id_map:
        raise ValueError(f"Ring '{r}' missing in ring_id_map.")
    septal_parts, freewall_parts = ring_layout[r]
    ids = ring_id_map[r]
    n_expected = int(septal_parts) + int(freewall_parts)
    if len(ids) != n_expected:
        raise ValueError(
            f"Ring '{r}' expects {n_expected} ids from layout but got {len(ids)}."
        )
    seg = np.full(t.shape, -1, dtype=np.int32)
    bounds = build_ring_boundaries(grooves, septal_parts, freewall_parts)
    for seg_id, (s, e) in zip(ids, bounds):
        seg[cw_is_between(t, s, e)] = int(seg_id)
    return seg


def infer_ring_from_longitudinal(
    z: np.ndarray,
    z_apical_mid: float,
    z_mid_basal: float,
    z_apex_cap: float,
) -> np.ndarray:
    out = np.empty(z.shape, dtype=object)
    out[:] = "basal"
    out[z < z_mid_basal] = "mid"
    out[z < z_apical_mid] = "apical"
    out[z <= z_apex_cap] = "apex"
    return out
