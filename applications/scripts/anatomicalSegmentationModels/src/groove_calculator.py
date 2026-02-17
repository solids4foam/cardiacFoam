"""Groove-angle estimation from masks and LV/RV interface geometry."""

from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

from .calc_utils import TAU, circular_mean, wrap_angle
from .mesh_marking import extract_intraventricular_interface_components


def _cw_dist(start: float, end: float) -> float:
    d = float(wrap_angle(end)) - float(wrap_angle(start))
    if d < 0.0:
        d += TAU
    return d


def _cw_dist_vec(start: float, theta: np.ndarray) -> np.ndarray:
    return np.mod(np.asarray(theta) - float(wrap_angle(start)), TAU)


def groove_angle_from_point_mask(
    rotational_theta: np.ndarray,
    mask: np.ndarray,
) -> float:
    """Estimate one groove angle as circular mean over selected points."""
    if mask.dtype != bool:
        mask = mask.astype(bool)
    sel = rotational_theta[mask]
    if sel.size == 0:
        raise ValueError("Groove mask selected 0 points; cannot estimate groove angle.")
    return circular_mean(sel)


def infer_groove_angles_debug_from_intraventricular_interface(
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
) -> Tuple[float, float, Dict[str, Any]]:
    """Infer (anterior, posterior) groove angles and debug point ids."""
    theta = np.asarray(rotational_theta).reshape(-1)
    intrav = np.asarray(intraventricular).reshape(-1)
    if theta.shape[0] != intrav.shape[0]:
        raise ValueError("rotational_theta and intraventricular must have same length.")

    comp_lists = extract_intraventricular_interface_components(
        mesh=mesh,
        intraventricular=intrav,
        longitudinal_z=longitudinal_z,
        z_min_exclusive=z_min_exclusive,
        z_cap_max=z_cap_max,
        lv_tag=lv_tag,
        rv_tag=rv_tag,
        tag_atol=tag_atol,
        min_component_points=min_component_points,
    )
    if len(comp_lists) == 0:
        raise ValueError(
            "Could not detect LV/RV interface elements containing both LV and RV tags."
        )

    mean_angles: List[float] = []
    debug: Dict[str, Any] = {
        "source": "intraventricular_interface",
        "z_min_exclusive": z_min_exclusive,
        "z_cap_max": z_cap_max,
    }
    ant_pts = np.array([], dtype=np.int64)
    post_pts = np.array([], dtype=np.int64)
    if len(comp_lists) >= 2:
        comp_lists = comp_lists[:2]
        debug["method"] = "two_components"
        for comp in comp_lists:
            mean_angles.append(circular_mean(theta[comp]))
        c0, c1 = comp_lists[0], comp_lists[1]
    else:
        comp_theta = theta[comp_lists[0]]
        d = _cw_dist_vec(rot_start, comp_theta)
        ant_mask = d < np.pi
        post_mask = ~ant_mask
        if np.any(ant_mask) and np.any(post_mask):
            debug["method"] = "single_component_split"
            mean_angles.append(circular_mean(comp_theta[ant_mask]))
            mean_angles.append(circular_mean(comp_theta[post_mask]))
            c0 = comp_lists[0][ant_mask]
            c1 = comp_lists[0][post_mask]
        else:
            debug["method"] = "single_component_extrema"
            i_min = int(np.argmin(d))
            i_max = int(np.argmax(d))
            mean_angles.append(float(wrap_angle(comp_theta[i_min])))
            mean_angles.append(float(wrap_angle(comp_theta[i_max])))
            c0 = np.asarray([comp_lists[0][i_min]], dtype=np.int64)
            c1 = np.asarray([comp_lists[0][i_max]], dtype=np.int64)

    a0, a1 = float(mean_angles[0]), float(mean_angles[1])
    if _cw_dist(rot_start, a0) <= _cw_dist(rot_start, a1):
        anterior, posterior = a0, a1
        ant_pts, post_pts = c0, c1
    else:
        anterior, posterior = a1, a0
        ant_pts, post_pts = c1, c0

    debug["anterior_point_ids"] = np.asarray(ant_pts, dtype=np.int64)
    debug["posterior_point_ids"] = np.asarray(post_pts, dtype=np.int64)
    debug["component_count"] = len(comp_lists)
    return anterior, posterior, debug


def infer_groove_angles_from_intraventricular_interface(
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
) -> Tuple[float, float]:
    """Infer (anterior, posterior) groove angles from LV/RV interface."""
    ant, post, _ = infer_groove_angles_debug_from_intraventricular_interface(
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
    return ant, post
