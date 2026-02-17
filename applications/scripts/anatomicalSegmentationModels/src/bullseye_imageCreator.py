from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Wedge


def _auto_ring_id_map(ring_layout):
    out = {}
    cur = 1
    for ring in ("basal", "mid", "apical"):
        septal_parts, freewall_parts = ring_layout[ring]
        n = int(septal_parts) + int(freewall_parts)
        out[ring] = list(range(cur, cur + n))
        cur += n
    return out


def draw_bullseye_reference(
    reference_name,
    reference_cfg,
    out_path=None,
    ring_widths=(0.28, 0.24, 0.18),  # (basal, mid, apical)
    start_angle_deg=90,              # top then counterclockwise
    color_mode="segment",            # "segment" or "artery"
    show_labels=True,
    dpi=200,
):
    """Draw one bullseye based on a division reference config."""
    ring_layout = reference_cfg["ring_layout"]
    ring_id_map = reference_cfg.get("ring_id_map") or _auto_ring_id_map(ring_layout)
    use_apex_cap = bool(reference_cfg.get("use_apex_cap", False))
    apex_id = reference_cfg.get("apex_id", None)
    colors = list(reference_cfg.get("colors", ["#cccccc"]))
    artery_colors = dict(reference_cfg.get("artery_colors", {}))
    segment_artery_map = dict(reference_cfg.get("segment_artery_map", {}))

    R = 1.0
    w_basal, w_mid, w_apical = ring_widths
    r3_outer = R
    r3_inner = R - w_basal
    r2_outer = r3_inner
    r2_inner = r2_outer - w_mid
    r1_outer = r2_inner
    r1_inner = r1_outer - w_apical

    # Default inner radius for apical wedges (used only when apex cap exists)
    r0 = r1_inner

    # If there is no apex cap (or no apex_id), apical wedges should reach the center (no hole)
    if (not use_apex_cap) or (apex_id is None):
        r0 = 0.0

    rings = [
        ("basal", (r3_inner, r3_outer)),
        ("mid", (r2_inner, r2_outer)),
        ("apical", (r0, r1_outer)),  # NOTE: apical now goes to r0 (0 if no apex cap)
    ]

    fig, ax = plt.subplots(figsize=(6, 6), dpi=dpi)
    ax.set_aspect("equal")
    ax.axis("off")

    for ring_name, (r_in, r_out) in rings:
        ids = ring_id_map[ring_name]
        n = len(ids)
        dtheta = 360.0 / n
        for k, seg_id in enumerate(ids):
            # Keep segment ID order consistent with 3D tagging order:
            # from start angle, progressing by increasing angle.
            theta1 = start_angle_deg + k * dtheta
            theta2 = start_angle_deg + (k + 1) * dtheta
            if (
                color_mode == "artery"
                and int(seg_id) in segment_artery_map
                and segment_artery_map[int(seg_id)] in artery_colors
            ):
                face = artery_colors[segment_artery_map[int(seg_id)]]
            else:
                face = colors[(int(seg_id) - 1) % len(colors)]
            wedge = Wedge(
                center=(0, 0),
                r=r_out,
                theta1=theta1,
                theta2=theta2,
                width=(r_out - r_in),
                facecolor=face,
                edgecolor="black",
                linewidth=1.2,
            )
            ax.add_patch(wedge)
            if show_labels:
                theta_mid = np.deg2rad((theta1 + theta2) / 2.0)
                r_mid = (r_in + r_out) / 2.0
                x = r_mid * np.cos(theta_mid)
                y = r_mid * np.sin(theta_mid)
                ax.text(
                    x, y, str(seg_id), ha="center", va="center",
                    fontsize=12, fontweight="bold", color="black"
                )

    # Draw apex cap only when enabled AND an apex_id is provided
    if use_apex_cap and (apex_id is not None):
        if (
            color_mode == "artery"
            and int(apex_id) in segment_artery_map
            and segment_artery_map[int(apex_id)] in artery_colors
        ):
            face = artery_colors[segment_artery_map[int(apex_id)]]
        else:
            face = colors[(int(apex_id) - 1) % len(colors)]
        apex = Circle((0, 0), radius=r1_inner, facecolor=face, edgecolor="black", linewidth=1.2)
        ax.add_patch(apex)
        if show_labels:
            ax.text(0, 0, str(apex_id), ha="center", va="center", fontsize=12, fontweight="bold")

    pad = 0.0
    ax.set_xlim(-(R + pad), (R + pad))
    ax.set_ylim(-(R + pad), (R + pad))
    if out_path:
        fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
        ax.set_position([0.0, 0.0, 1.0, 1.0])
        fig.savefig(
            out_path,
            dpi=dpi,
            bbox_inches=None,
            pad_inches=0.0,
            facecolor="white",
        )
    return fig, ax


def draw_all_bullseyes(division_references, output_dir=None, show=False):
    """Draw all available bullseye references from config."""
    out_dir = None if output_dir is None else Path(output_dir)
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
    figs = []
    for name, cfg in division_references.items():
        out_path = None if out_dir is None else out_dir / f"{name}_bullseye.png"
        fig, _ = draw_bullseye_reference(
            name, cfg, out_path=str(out_path) if out_path else None, color_mode="segment"
        )
        figs.append(fig)
        if cfg.get("artery_colors") and cfg.get("segment_artery_map") and out_dir is not None:
            out_path_art = out_dir / f"{name}_bullseye_artery.png"
            fig2, _ = draw_bullseye_reference(
                name, cfg, out_path=str(out_path_art), color_mode="artery"
            )
            figs.append(fig2)
    if show:
        plt.show()
    return figs


def ensure_bullseye_png(
    division_reference,
    division_references,
    output_dir,
    overwrite=True,
):
    """Generate and return the bullseye PNG path for one active reference."""
    if division_reference not in division_references:
        raise KeyError(f"Unknown division reference '{division_reference}'.")
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{division_reference}_bullseye.png"
    if overwrite or (not out_path.exists()):
        draw_bullseye_reference(
            division_reference,
            division_references[division_reference],
            out_path=str(out_path),
            color_mode="segment",
            show_labels=True,
        )
        cfg = division_references[division_reference]
        if cfg.get("artery_colors") and cfg.get("segment_artery_map"):
            out_path_art = out_dir / f"{division_reference}_bullseye_artery.png"
            draw_bullseye_reference(
                division_reference,
                cfg,
                out_path=str(out_path_art),
                color_mode="artery",
                show_labels=True,
            )
        plt.close("all")
    return out_path


if __name__ == "__main__":
    from configs.lv_division_config import division_references

    draw_all_bullseyes(division_references, output_dir=".", show=True)
