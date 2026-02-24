from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping


DEFAULT_AHA17_REFERENCE_CFG: dict[str, object] = {
    "ring_layout": {
        "basal": [2, 4],
        "mid": [2, 4],
        "apical": [1, 3],
    },
    "ring_id_map": {
        "basal": [1, 2, 3, 4, 5, 6],
        "mid": [7, 8, 9, 10, 11, 12],
        "apical": [13, 14, 15, 16],
    },
    "apex_id": 17,
    "use_apex_cap": True,
}


@dataclass(frozen=True)
class BullseyeSegmentGeometry:
    segment_id: int
    ring_name: str
    r_inner: float
    r_outer: float
    theta1_deg: float
    theta2_deg: float
    is_apex_cap: bool = False


def _auto_ring_id_map(ring_layout: Mapping[str, Iterable[int]]) -> dict[str, list[int]]:
    out: dict[str, list[int]] = {}
    cur = 1
    for ring in ("basal", "mid", "apical"):
        parts = list(ring_layout[ring])
        n = int(parts[0]) + int(parts[1])
        out[ring] = list(range(cur, cur + n))
        cur += n
    return out


def segment_ids_from_reference(reference_cfg: Mapping[str, object]) -> list[int]:
    ring_layout = dict(reference_cfg["ring_layout"])
    ring_id_map = reference_cfg.get("ring_id_map")
    if ring_id_map is None:
        ring_id_map = _auto_ring_id_map(ring_layout)
    ring_id_map = dict(ring_id_map)

    ordered: list[int] = []
    for ring in ("basal", "mid", "apical"):
        ordered.extend(int(x) for x in ring_id_map[ring])

    use_apex_cap = bool(reference_cfg.get("use_apex_cap", False))
    apex_id = reference_cfg.get("apex_id")
    if use_apex_cap and (apex_id is not None):
        ordered.append(int(apex_id))

    return ordered


def build_bullseye_segment_geometry(
    reference_cfg: Mapping[str, object],
    ring_widths: tuple[float, float, float] = (0.28, 0.24, 0.18),
    start_angle_deg: float = 120.0,
    apical_start_angle_deg: float | None = 45.0,
) -> list[BullseyeSegmentGeometry]:
    ring_layout = dict(reference_cfg["ring_layout"])
    ring_id_map = reference_cfg.get("ring_id_map")
    if ring_id_map is None:
        ring_id_map = _auto_ring_id_map(ring_layout)
    ring_id_map = dict(ring_id_map)

    use_apex_cap = bool(reference_cfg.get("use_apex_cap", False))
    apex_id = reference_cfg.get("apex_id")

    r_outer = 1.0
    w_basal, w_mid, w_apical = ring_widths
    r_basal_in = r_outer - w_basal
    r_mid_out = r_basal_in
    r_mid_in = r_mid_out - w_mid
    r_apical_out = r_mid_in
    r_apical_in = r_apical_out - w_apical
    if (not use_apex_cap) or (apex_id is None):
        r_apical_in = 0.0

    rings = [
        ("basal", r_basal_in, r_outer),
        ("mid", r_mid_in, r_mid_out),
        ("apical", r_apical_in, r_apical_out),
    ]

    out: list[BullseyeSegmentGeometry] = []
    for ring_name, r_in, r_out in rings:
        ids = [int(x) for x in ring_id_map[ring_name]]
        n = len(ids)
        dtheta = 360.0 / float(n)
        ring_start = float(start_angle_deg)
        if ring_name == "apical" and apical_start_angle_deg is not None:
            ring_start = float(apical_start_angle_deg)
        for i, seg_id in enumerate(ids):
            theta1 = ring_start + i * dtheta
            theta2 = ring_start + (i + 1) * dtheta
            out.append(
                BullseyeSegmentGeometry(
                    segment_id=int(seg_id),
                    ring_name=ring_name,
                    r_inner=float(r_in),
                    r_outer=float(r_out),
                    theta1_deg=float(theta1),
                    theta2_deg=float(theta2),
                    is_apex_cap=False,
                )
            )

    if use_apex_cap and (apex_id is not None):
        out.append(
            BullseyeSegmentGeometry(
                segment_id=int(apex_id),
                ring_name="apex",
                r_inner=0.0,
                r_outer=float(r_apical_in),
                theta1_deg=0.0,
                theta2_deg=360.0,
                is_apex_cap=True,
            )
        )
    return out


def validate_segment_percentages(
    segment_percentages: Mapping[int, float],
    reference_cfg: Mapping[str, object],
    allow_partial: bool = True,
    min_value: float = 0.0,
    max_value: float = 100.0,
) -> None:
    if max_value <= min_value:
        raise ValueError("max_value must be greater than min_value.")

    expected_ids = set(segment_ids_from_reference(reference_cfg))
    provided_ids = set(int(k) for k in segment_percentages.keys())

    unknown = sorted(provided_ids - expected_ids)
    if unknown:
        raise ValueError(f"Unknown segment ids: {unknown}")

    if not allow_partial:
        missing = sorted(expected_ids - provided_ids)
        if missing:
            raise ValueError(f"Missing percentage values for segments: {missing}")

    for seg_id, value in segment_percentages.items():
        v = float(value)
        if v < min_value or v > max_value:
            raise ValueError(
                f"Segment {int(seg_id)} has invalid percentage {v}. "
                f"Expected [{min_value}, {max_value}]."
            )


def collect_segment_percentages_interactive(
    reference_cfg: Mapping[str, object],
    initial_percentages: Mapping[int, float] | None = None,
    title: str = "Bullseye Segment Input",
    prompt_missing_after_close: bool = True,
    input_mode: str = "widget",
    start_angle_deg: float = 120.0,
    apical_start_angle_deg: float | None = 45.0,
    ring_widths: tuple[float, float, float] = (0.28, 0.24, 0.18),
    max_percentage_value: float = 10.0,
) -> dict[int, float]:
    """
    Open a clickable bullseye plot.

    User flow:
    - Click a segment.
    - Enter percentage in the plot widget panel (default) or terminal fallback.
    - Segment label updates to show "<id>\\n<value>%".
    - Close window when done.
    - Optionally fill missing segments (0.0 auto-fill in widget mode).
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        from matplotlib.patches import Circle, Wedge
    except Exception as exc:
        raise RuntimeError(
            "Matplotlib is required for interactive bullseye input."
        ) from exc

    input_mode = str(input_mode).strip().lower()
    if input_mode not in {"widget", "terminal"}:
        raise ValueError("input_mode must be 'widget' or 'terminal'.")
    max_percentage_value = float(max_percentage_value)
    if max_percentage_value <= 0.0:
        raise ValueError("max_percentage_value must be > 0.")
    range_hint = f"[0, {max_percentage_value:g}]"

    geometry = build_bullseye_segment_geometry(
        reference_cfg=reference_cfg,
        ring_widths=ring_widths,
        start_angle_deg=start_angle_deg,
        apical_start_angle_deg=apical_start_angle_deg,
    )
    ordered_ids = [g.segment_id for g in geometry]

    percentages: dict[int, float] = {}
    if initial_percentages is not None:
        percentages = {int(k): float(v) for k, v in initial_percentages.items()}
    validate_segment_percentages(
        percentages,
        reference_cfg,
        allow_partial=True,
        min_value=0.0,
        max_value=max_percentage_value,
    )

    fig, ax = plt.subplots(figsize=(8.8, 7), dpi=140)
    center_x = -0.24
    center_y = 0.0
    direction_radius_side = 1.30
    direction_radius_vertical = 1.18
    if input_mode == "widget":
        fig.subplots_adjust(right=0.74)
    fig.suptitle(title, fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.67, 1.13)
    ax.set_ylim(-1.25, 1.25)
    ax.text(
        center_x - direction_radius_side,
        center_y,
        "SEPTAL",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        center_x + direction_radius_side,
        center_y,
        "LATERAL",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        center_x,
        center_y + direction_radius_vertical,
        "POSTERIOR",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        center_x,
        center_y - direction_radius_vertical,
        "ANTERIOR",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

    normalizer = Normalize(vmin=0.0, vmax=max_percentage_value)
    cmap = plt.get_cmap("YlOrRd")
    mappable = ScalarMappable(norm=normalizer, cmap=cmap)

    patches: dict[int, object] = {}
    labels: dict[int, object] = {}
    default_face = "#d9d9d9"
    selected = {"seg_id": None}

    def _segment_text(seg_id: int) -> str:
        if seg_id in percentages:
            return f"{seg_id}\n{percentages[seg_id]:.1f}%"
        return str(seg_id)

    for g in geometry:
        if g.is_apex_cap:
            patch = Circle(
                (center_x, center_y),
                radius=g.r_outer,
                facecolor=default_face,
                edgecolor="black",
                linewidth=1.2,
            )
            x, y = center_x, center_y
        else:
            patch = Wedge(
                center=(center_x, center_y),
                r=g.r_outer,
                theta1=g.theta1_deg,
                theta2=g.theta2_deg,
                width=(g.r_outer - g.r_inner),
                facecolor=default_face,
                edgecolor="black",
                linewidth=1.2,
            )
            theta_mid = (g.theta1_deg + g.theta2_deg) * 0.5
            r_mid = (g.r_inner + g.r_outer) * 0.5
            import math

            x = center_x + r_mid * math.cos(math.radians(theta_mid))
            y = center_y + r_mid * math.sin(math.radians(theta_mid))

        patch.set_picker(True)
        setattr(patch, "_seg_id", int(g.segment_id))
        ax.add_patch(patch)
        txt = ax.text(
            x,
            y,
            _segment_text(g.segment_id),
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color="black",
        )
        patches[g.segment_id] = patch
        labels[g.segment_id] = txt

    def _apply_colors_and_labels() -> None:
        for seg_id in ordered_ids:
            patch = patches[seg_id]
            label = labels[seg_id]
            if seg_id in percentages:
                patch.set_facecolor(cmap(normalizer(percentages[seg_id])))
            else:
                patch.set_facecolor(default_face)
            label.set_text(_segment_text(seg_id))
        fig.canvas.draw_idle()

    if input_mode == "widget":
        from matplotlib.widgets import Button, TextBox

        status_ax = fig.add_axes([0.76, 0.77, 0.22, 0.14])
        status_ax.axis("off")
        status_text = status_ax.text(
            0.0,
            0.9,
            "Selected: none",
            fontsize=10,
            va="top",
        )

        help_ax = fig.add_axes([0.76, 0.56, 0.22, 0.16])
        help_ax.axis("off")
        help_ax.text(
            0.0,
            1.0,
            (
                f"Per-segment value (%), not global sum\n"
                f"Allowed range: {range_hint}\n"
                "Click segment -> type -> Enter/Apply"
            ),
            fontsize=9,
            va="top",
        )

        # Input controls placed lower-right for unobstructed reading.
        textbox_ax = fig.add_axes([0.78, 0.30, 0.18, 0.06])
        value_box = TextBox(textbox_ax, "Value % ", initial="")

        apply_ax = fig.add_axes([0.78, 0.23, 0.08, 0.05])
        clear_ax = fig.add_axes([0.88, 0.23, 0.08, 0.05])
        apply_btn = Button(apply_ax, "Apply")
        clear_btn = Button(clear_ax, "Clear")

        def _set_status(msg: str) -> None:
            status_text.set_text(msg)
            fig.canvas.draw_idle()

        def _select_segment(seg_id: int) -> None:
            selected["seg_id"] = int(seg_id)
            current = percentages.get(int(seg_id))
            value_box.set_val("" if current is None else f"{current:.1f}")
            _set_status(f"Selected: {seg_id}")

        def _apply_value(raw_value: str | None = None) -> None:
            seg_id = selected["seg_id"]
            if seg_id is None:
                _set_status("Select a segment first.")
                return
            text = (value_box.text if raw_value is None else raw_value).strip()
            if text == "":
                _set_status("Empty value. Use Clear to remove.")
                return
            try:
                value = float(text)
            except ValueError:
                _set_status(f"Invalid value '{text}'.")
                return
            if value < 0.0 or value > max_percentage_value:
                _set_status(f"Out of range. Use {range_hint}.")
                return
            percentages[int(seg_id)] = float(value)
            _apply_colors_and_labels()
            _set_status(f"Segment {seg_id} = {value:.1f}%")

        def _clear_value(_evt=None) -> None:
            seg_id = selected["seg_id"]
            if seg_id is None:
                _set_status("Select a segment first.")
                return
            if int(seg_id) in percentages:
                del percentages[int(seg_id)]
            value_box.set_val("")
            _apply_colors_and_labels()
            _set_status(f"Segment {seg_id} cleared.")

        value_box.on_submit(_apply_value)
        apply_btn.on_clicked(lambda _evt: _apply_value())
        clear_btn.on_clicked(_clear_value)

    def _on_pick(event) -> None:
        patch = getattr(event, "artist", None)
        seg_id = getattr(patch, "_seg_id", None)
        if seg_id is None:
            return

        if input_mode == "widget":
            selected["seg_id"] = int(seg_id)
            current = percentages.get(int(seg_id))
            # Text is set in widget mode; no terminal prompt needed.
            try:
                value_box.set_val("" if current is None else f"{current:.1f}")  # type: ignore[name-defined]
                _set_status(f"Selected: {seg_id}")  # type: ignore[name-defined]
            except Exception:
                pass
            fig.canvas.draw_idle()
            return

        raw = input(f"Segment {seg_id} percentage {range_hint}, blank=keep: ").strip()
        if raw == "":
            return
        try:
            value = float(raw)
        except ValueError:
            print(f"Invalid input for segment {seg_id}. Please enter a numeric value.")
            return
        if value < 0.0 or value > max_percentage_value:
            print(f"Value out of range for segment {seg_id}. Use {range_hint}.")
            return

        percentages[int(seg_id)] = float(value)
        _apply_colors_and_labels()

    fig.canvas.mpl_connect("pick_event", _on_pick)
    def _on_close(_event) -> None:
        missing = [seg_id for seg_id in ordered_ids if seg_id not in percentages]
        if missing:
            if prompt_missing_after_close:
                print(
                    "Window closed with missing segments: "
                    f"{missing}. They will be auto-filled with 0.0."
                )
            else:
                print(
                    "Window closed with missing segments: "
                    f"{missing}. Returning partial values."
                )

    fig.canvas.mpl_connect("close_event", _on_close)
    fig.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04, label="Per-segment value (%)")
    _apply_colors_and_labels()
    plt.show()

    if prompt_missing_after_close:
        missing = [seg_id for seg_id in ordered_ids if seg_id not in percentages]
        if input_mode == "widget":
            for seg_id in missing:
                percentages[seg_id] = 0.0
        else:
            for seg_id in missing:
                while True:
                    raw = input(
                        f"Missing segment {seg_id}. Enter percentage {range_hint} (or blank=0): "
                    ).strip()
                    if raw == "":
                        percentages[seg_id] = 0.0
                        break
                    try:
                        value = float(raw)
                    except ValueError:
                        print("Invalid number. Try again.")
                        continue
                    if value < 0.0 or value > max_percentage_value:
                        print(f"Out of range. Use {range_hint}.")
                        continue
                    percentages[seg_id] = value
                    break

    validate_segment_percentages(
        percentages,
        reference_cfg,
        allow_partial=not prompt_missing_after_close,
        min_value=0.0,
        max_value=max_percentage_value,
    )
    return {seg_id: percentages[seg_id] for seg_id in ordered_ids if seg_id in percentages}
