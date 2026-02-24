from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .interactive_bullseye import build_bullseye_segment_geometry, segment_ids_from_reference
from .mesh_context import build_segment_mesh_context, save_mesh_ascii
from .vtk_bullseye_workbench import (
    load_division_reference,
    map_percentages_to_segments,
)


def launch_unified_workbench(
    vtk_path: Path,
    output_vtk_path: Path | None = None,
    reference_name: str | None = None,
    config_source: str = "anatomical",
    config_file: Path | None = None,
    segment_field: str | None = None,
    segment_association: str = "auto",
    purkinje_field_name: str = "purkinje_density",
    pmj_field_name: str = "pmj_per_mm2",
    thickness_bundles_field_name: str = "thickness_bundles",
    endocardial_area_field_name: str = "endocardial_surface_area_by_tag",
    endocardial_selector_field: str | None = None,
    endocardial_selector_value: Any = 1,
    purkinje_max_value: float = 1.0,
    pmj_max_value: float = 1.0,
    thickness_bundles_max_value: float = 200.0,
    prompt_missing_after_close: bool = True,
) -> dict[str, dict[int, float]]:
    """
    Open one window with 2 synchronized panels:
    - Interactive bullseye editor (Purkinje, PMJ, Thickness Bundles)
    - Mapped 3D values (with segment-id labels at barycenters)
    """
    try:
        import numpy as np
        import pyvista as pv
        from pyvistaqt import QtInteractor
    except Exception as exc:
        raise RuntimeError(
            "Unified GUI requires numpy + pyvista + pyvistaqt. "
            "Install with: pip install pyvistaqt"
        ) from exc

    try:
        from PySide6 import QtCore, QtWidgets
    except Exception:
        try:
            from PyQt5 import QtCore, QtWidgets  # type: ignore
        except Exception as exc:
            raise RuntimeError(
                "Unified GUI requires PySide6 or PyQt5."
            ) from exc

    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        from matplotlib.figure import Figure
        from matplotlib.patches import Circle, Wedge
    except Exception as exc:
        raise RuntimeError("Unified GUI requires matplotlib with Qt backend.") from exc

    purkinje_max_value = float(purkinje_max_value)
    pmj_max_value = float(pmj_max_value)
    thickness_bundles_max_value = float(thickness_bundles_max_value)
    if purkinje_max_value <= 0.0:
        raise ValueError("purkinje_max_value must be > 0.")
    if pmj_max_value <= 0.0:
        raise ValueError("pmj_max_value must be > 0.")
    if thickness_bundles_max_value <= 0.0:
        raise ValueError("thickness_bundles_max_value must be > 0.")

    if not vtk_path.exists():
        raise FileNotFoundError(f"VTK file not found: {vtk_path}")

    division = load_division_reference(
        reference_name=reference_name,
        config_source=config_source,
        config_file=config_file,
    )

    mesh_ctx = build_segment_mesh_context(
        vtk_path,
        segment_field=segment_field,
        segment_association=segment_association,
        endocardial_selector_field=endocardial_selector_field,
        endocardial_selector_value=endocardial_selector_value,
    )
    mesh = mesh_ctx.mesh
    selection = mesh_ctx.selection
    segment_ids = np.asarray(mesh_ctx.segment_ids).reshape(-1).astype(int)
    area_result = mesh_ctx.area_result
    endocardial_area_by_segment = dict(area_result.area_by_segment)
    endocardial_area_mapped = area_result.mapped_values
    endocardial_selector_source = str(area_result.selector_source)

    expected_ids = sorted(segment_ids_from_reference(division.reference_cfg))
    present_ids = sorted(int(x) for x in np.unique(segment_ids) if int(x) > 0)
    unknown_ids = sorted(set(present_ids) - set(expected_ids))
    if unknown_ids:
        print(
            "Warning: VTK contains segment IDs not in selected bullseye reference: "
            f"{unknown_ids}"
        )

    geometry = build_bullseye_segment_geometry(division.reference_cfg)
    ordered_ids = [g.segment_id for g in geometry]

    class UnifiedWorkbenchWindow(QtWidgets.QMainWindow):
        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle(f"Unified Bullseye Workbench ({division.reference_name})")
            self.resize(1480, 880)

            self.mesh = mesh
            self.out_mesh = mesh.copy(deep=True)
            self._ordered_ids = ordered_ids
            self._geometry = geometry
            self._cmap = plt.get_cmap("YlOrRd")
            self._default_face = "#d9d9d9"
            self._output_path = output_vtk_path
            self._closing_after_save = False

            self._tag_defs = {
                "purkinje": {
                    "display": "Purkinje Density",
                    "field": str(purkinje_field_name),
                    "max": float(purkinje_max_value),
                    "unit": "",
                },
                "pmj": {
                    "display": "PMJ /mm2",
                    "field": str(pmj_field_name),
                    "max": float(pmj_max_value),
                    "unit": "/mm2",
                },
                "thickness_bundles": {
                    "display": "Thickness Bundles",
                    "field": str(thickness_bundles_field_name),
                    "max": float(thickness_bundles_max_value),
                    "unit": "",
                },
            }
            self._tab_order = ["purkinje", "pmj", "thickness_bundles"]
            self._active_tag_key = self._tab_order[0]
            self._tag_states: dict[str, dict[str, Any]] = {
                key: {
                    "edit_values": {},
                    "preview_values": {},
                    "selected_seg_id": None,
                    "pending": False,
                }
                for key in self._tab_order
            }
            self._ui: dict[str, dict[str, Any]] = {}

            self._segment_label_points, self._segment_label_texts = self._compute_segment_label_data()
            self._build_ui()
            for key in self._tab_order:
                self._render_bullseye(key)
            self._on_tab_changed(0)

        def _compute_segment_label_data(self):
            points: list[Any] = []
            labels: list[str] = []
            if selection.association == "point":
                mesh_points = np.asarray(self.mesh.points)
                for seg_id in sorted(int(x) for x in np.unique(segment_ids) if int(x) > 0):
                    mask = segment_ids == int(seg_id)
                    if not np.any(mask):
                        continue
                    pts = mesh_points[mask]
                    if pts.size == 0:
                        continue
                    points.append(pts.mean(axis=0))
                    labels.append(str(int(seg_id)))
            else:
                centers = self.mesh.cell_centers().points
                for seg_id in sorted(int(x) for x in np.unique(segment_ids) if int(x) > 0):
                    mask = segment_ids == int(seg_id)
                    if not np.any(mask):
                        continue
                    pts = centers[mask]
                    if pts.size == 0:
                        continue
                    points.append(pts.mean(axis=0))
                    labels.append(str(int(seg_id)))
            if not points:
                return np.empty((0, 3)), []
            return np.asarray(points), labels

        def _build_ui(self) -> None:
            cw = QtWidgets.QWidget()
            root = QtWidgets.QHBoxLayout(cw)
            root.setContentsMargins(8, 8, 8, 8)
            root.setSpacing(10)
            self.setCentralWidget(cw)

            left_box = QtWidgets.QGroupBox("Bullseye Editors")
            left_layout = QtWidgets.QVBoxLayout(left_box)
            self.tab_widget = QtWidgets.QTabWidget()
            left_layout.addWidget(self.tab_widget)

            for key in self._tab_order:
                tab = QtWidgets.QWidget()
                tab_layout = QtWidgets.QVBoxLayout(tab)
                figure = Figure(figsize=(6.8, 6.8), dpi=110)
                canvas = FigureCanvas(figure)
                ax = figure.add_subplot(111)
                tab_layout.addWidget(canvas, 1)

                control_grid = QtWidgets.QGridLayout()
                selected_label = QtWidgets.QLabel("Selected: none")
                max_val = float(self._tag_defs[key]["max"])
                unit = str(self._tag_defs[key]["unit"])
                unit_text = "" if unit == "" else f" ({unit})"
                status_label = QtWidgets.QLabel(
                    f"Per-segment value{unit_text} in [0, {max_val:g}] | no global sum constraint"
                )
                value_edit = QtWidgets.QLineEdit("")
                value_edit.setPlaceholderText("Value")
                value_edit.returnPressed.connect(lambda k=key: self._set_value(k))
                set_btn = QtWidgets.QPushButton("Set Value")
                clear_btn = QtWidgets.QPushButton("Clear")
                apply3d_btn = QtWidgets.QPushButton("Apply to 3D")
                fill_btn = QtWidgets.QPushButton("Fill Missing With 0")
                save_exit_btn = QtWidgets.QPushButton("Save + Exit")

                set_btn.clicked.connect(lambda _=False, k=key: self._set_value(k))
                clear_btn.clicked.connect(lambda _=False, k=key: self._clear_value(k))
                apply3d_btn.clicked.connect(lambda _=False, k=key: self._apply_to_3d(k))
                fill_btn.clicked.connect(lambda _=False, k=key: self._fill_missing_with_zero(k))
                save_exit_btn.clicked.connect(self._save_and_exit)

                control_grid.addWidget(selected_label, 0, 0, 1, 2)
                control_grid.addWidget(status_label, 1, 0, 1, 2)
                control_grid.addWidget(value_edit, 2, 0, 1, 2)
                control_grid.addWidget(set_btn, 3, 0, 1, 1)
                control_grid.addWidget(clear_btn, 3, 1, 1, 1)
                control_grid.addWidget(apply3d_btn, 4, 0, 1, 1)
                control_grid.addWidget(fill_btn, 4, 1, 1, 1)
                control_grid.addWidget(save_exit_btn, 5, 0, 1, 2)
                tab_layout.addLayout(control_grid)

                canvas.mpl_connect(
                    "pick_event",
                    lambda event, k=key: self._on_pick_segment(k, event),
                )

                self._ui[key] = {
                    "figure": figure,
                    "canvas": canvas,
                    "ax": ax,
                    "selected_label": selected_label,
                    "status_label": status_label,
                    "value_edit": value_edit,
                    "patches": {},
                    "labels": {},
                    "normalizer": Normalize(vmin=0.0, vmax=max_val),
                }

                self.tab_widget.addTab(tab, self._tag_defs[key]["display"])

            self.tab_widget.currentChanged.connect(self._on_tab_changed)
            root.addWidget(left_box, 1)

            right_box = QtWidgets.QGroupBox("Mapped 3D Values + Segment IDs")
            right_layout = QtWidgets.QVBoxLayout(right_box)
            self.right_view = QtInteractor(right_box)
            right_layout.addWidget(self.right_view.interactor)
            root.addWidget(right_box, 1)

        def _state(self, tag_key: str) -> dict[str, Any]:
            return self._tag_states[tag_key]

        def _ui_for(self, tag_key: str) -> dict[str, Any]:
            return self._ui[tag_key]

        def _segment_text(self, tag_key: str, seg_id: int) -> str:
            vals = self._state(tag_key)["edit_values"]
            if seg_id in vals:
                return f"{seg_id}\n{vals[seg_id]:.3f}"
            return str(seg_id)

        def _render_bullseye(self, tag_key: str) -> None:
            ui = self._ui_for(tag_key)
            ax = ui["ax"]
            center_x = -0.24
            center_y = 0.0
            direction_radius_side = 1.30
            direction_radius_vertical = 1.18
            ax.clear()
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
                fontsize=9,
                fontweight="bold",
            )
            ax.text(
                center_x + direction_radius_side,
                center_y,
                "LATERAL",
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
            )
            ax.text(
                center_x,
                center_y + direction_radius_vertical,
                "POSTERIOR",
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
            )
            ax.text(
                center_x,
                center_y - direction_radius_vertical,
                "ANTERIOR",
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
            )

            ui["patches"] = {}
            ui["labels"] = {}

            for g in self._geometry:
                if g.is_apex_cap:
                    patch = Circle(
                        (center_x, center_y),
                        radius=g.r_outer,
                        facecolor=self._default_face,
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
                        facecolor=self._default_face,
                        edgecolor="black",
                        linewidth=1.2,
                    )
                    theta_mid = (g.theta1_deg + g.theta2_deg) * 0.5
                    r_mid = (g.r_inner + g.r_outer) * 0.5
                    x = center_x + r_mid * np.cos(np.deg2rad(theta_mid))
                    y = center_y + r_mid * np.sin(np.deg2rad(theta_mid))

                patch.set_picker(True)
                setattr(patch, "_seg_id", int(g.segment_id))
                ax.add_patch(patch)
                txt = ax.text(
                    x,
                    y,
                    self._segment_text(tag_key, g.segment_id),
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                )
                ui["patches"][g.segment_id] = patch
                ui["labels"][g.segment_id] = txt

            mappable = ScalarMappable(norm=ui["normalizer"], cmap=self._cmap)
            unit = str(self._tag_defs[tag_key]["unit"])
            bar_label = "Per-segment value" if unit == "" else f"Per-segment value ({unit})"
            ui["figure"].colorbar(
                mappable,
                ax=ax,
                fraction=0.046,
                pad=0.04,
                label=bar_label,
            )
            self._refresh_bullseye_values(tag_key)

        def _refresh_bullseye_values(self, tag_key: str) -> None:
            ui = self._ui_for(tag_key)
            vals = self._state(tag_key)["edit_values"]
            for seg_id in self._ordered_ids:
                patch = ui["patches"][seg_id]
                label = ui["labels"][seg_id]
                if seg_id in vals:
                    patch.set_facecolor(self._cmap(ui["normalizer"](vals[seg_id])))
                else:
                    patch.set_facecolor(self._default_face)
                label.set_text(self._segment_text(tag_key, seg_id))
            ui["canvas"].draw_idle()

        def _set_status(self, tag_key: str, msg: str) -> None:
            self._ui_for(tag_key)["status_label"].setText(msg)

        def _on_pick_segment(self, tag_key: str, event) -> None:
            patch = getattr(event, "artist", None)
            seg_id = getattr(patch, "_seg_id", None)
            if seg_id is None:
                return
            state = self._state(tag_key)
            ui = self._ui_for(tag_key)
            state["selected_seg_id"] = int(seg_id)
            ui["selected_label"].setText(f"Selected: {state['selected_seg_id']}")
            current = state["edit_values"].get(state["selected_seg_id"])
            ui["value_edit"].setText("" if current is None else f"{current:.3f}")
            max_val = float(self._tag_defs[tag_key]["max"])
            self._set_status(
                tag_key,
                f"Segment {state['selected_seg_id']} selected. "
                f"Enter value in [0, {max_val:g}] and click Set Value.",
            )

        def _set_value(self, tag_key: str) -> None:
            state = self._state(tag_key)
            ui = self._ui_for(tag_key)
            seg_id = state["selected_seg_id"]
            if seg_id is None:
                self._set_status(tag_key, "Select a segment first.")
                return
            raw = ui["value_edit"].text().strip()
            if raw == "":
                self._set_status(tag_key, "Empty value. Use Clear to remove current segment value.")
                return
            try:
                value = float(raw)
            except ValueError:
                self._set_status(tag_key, f"Invalid numeric value '{raw}'.")
                return
            max_val = float(self._tag_defs[tag_key]["max"])
            if value < 0.0 or value > max_val:
                self._set_status(tag_key, f"Out of range. Use [0, {max_val:g}].")
                return

            state["edit_values"][int(seg_id)] = float(value)
            state["pending"] = True
            self._refresh_bullseye_values(tag_key)
            self._set_status(
                tag_key,
                f"Segment {seg_id} = {value:.3f}. Click 'Apply to 3D' to refresh preview.",
            )

        def _clear_value(self, tag_key: str) -> None:
            state = self._state(tag_key)
            ui = self._ui_for(tag_key)
            seg_id = state["selected_seg_id"]
            if seg_id is None:
                self._set_status(tag_key, "Select a segment first.")
                return
            state["edit_values"].pop(int(seg_id), None)
            state["pending"] = True
            ui["value_edit"].setText("")
            self._refresh_bullseye_values(tag_key)
            self._set_status(
                tag_key,
                f"Segment {seg_id} cleared. Click 'Apply to 3D' to refresh preview.",
            )

        def _missing_segments(self, tag_key: str) -> list[int]:
            vals = self._state(tag_key)["edit_values"]
            return [seg_id for seg_id in self._ordered_ids if seg_id not in vals]

        def _fill_missing_with_zero(self, tag_key: str) -> None:
            state = self._state(tag_key)
            missing = self._missing_segments(tag_key)
            for seg_id in missing:
                state["edit_values"][seg_id] = 0.0
            state["pending"] = True
            self._refresh_bullseye_values(tag_key)
            if missing:
                self._set_status(
                    tag_key,
                    f"Filled missing segments with 0.0: {missing}. Click 'Apply to 3D' to refresh preview.",
                )
            else:
                self._set_status(tag_key, "No missing segments.")

        def _apply_to_3d(self, tag_key: str) -> None:
            state = self._state(tag_key)
            state["preview_values"] = dict(state["edit_values"])
            state["pending"] = False
            if tag_key == self._active_tag_key:
                self._update_mapped_view()
            self._set_status(tag_key, "3D preview updated for this tag.")

        def _apply_all_pending_to_3d(self) -> None:
            for key in self._tab_order:
                if self._state(key)["pending"]:
                    self._state(key)["preview_values"] = dict(self._state(key)["edit_values"])
                    self._state(key)["pending"] = False
            self._update_mapped_view()

        def _mapped_mesh(self):
            self.out_mesh = self.mesh.copy(deep=True)
            for key in self._tab_order:
                field_name = self._tag_defs[key]["field"]
                mapped = map_percentages_to_segments(
                    segment_ids=segment_ids,
                    segment_percentages=self._state(key)["preview_values"],
                    default_value=0.0,
                )
                if selection.association == "point":
                    self.out_mesh.point_data[field_name] = mapped
                else:
                    self.out_mesh.cell_data[field_name] = mapped
            area_field = str(endocardial_area_field_name)
            if selection.association == "point":
                self.out_mesh.point_data[area_field] = endocardial_area_mapped
            else:
                self.out_mesh.cell_data[area_field] = endocardial_area_mapped
            return self.out_mesh

        def _update_mapped_view(self) -> None:
            out_mesh = self._mapped_mesh()
            key = self._active_tag_key
            field_name = self._tag_defs[key]["field"]
            display = self._tag_defs[key]["display"]
            active_max = float(self._tag_defs[key]["max"])
            self.right_view.clear()
            self.right_view.add_mesh(
                out_mesh,
                scalars=field_name,
                preference=selection.association,
                cmap="YlOrRd",
                clim=(0.0, active_max),
                show_edges=False,
            )
            self.right_view.add_text(
                f"{display} | field={field_name} ({selection.association})",
                position="upper_left",
                font_size=10,
            )
            if self._segment_label_points.shape[0] > 0:
                self.right_view.add_point_labels(
                    self._segment_label_points,
                    self._segment_label_texts,
                    point_size=0,
                    font_size=12,
                    text_color="black",
                    shape_color="white",
                    shape_opacity=0.25,
                    always_visible=True,
                )
            self.right_view.reset_camera()

        def _on_tab_changed(self, index: int) -> None:
            if index < 0 or index >= len(self._tab_order):
                return
            self._active_tag_key = self._tab_order[index]
            self._update_mapped_view()
            state = self._state(self._active_tag_key)
            if state["pending"]:
                self._set_status(
                    self._active_tag_key,
                    "Pending edits exist. Click 'Apply to 3D' to refresh preview.",
                )

        def _format_missing(self, missing_by_tag: dict[str, list[int]]) -> str:
            lines = []
            for key in self._tab_order:
                vals = missing_by_tag.get(key, [])
                if vals:
                    lines.append(f"{self._tag_defs[key]['display']}: {vals}")
            return "\n".join(lines)

        def _save_vtk_dialog(self) -> bool:
            default = str(self._output_path) if self._output_path is not None else str(
                vtk_path.with_name(f"{vtk_path.stem}_densities_mapped{vtk_path.suffix}")
            )
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self,
                "Save mapped VTK",
                default,
                "VTK files (*.vtk *.vtu *.vtp);;All files (*)",
            )
            if not path:
                return False
            save_mesh_ascii(self._mapped_mesh(), path)
            self._set_status(self._active_tag_key, f"Saved mapped VTK: {path}")
            return True

        def _save_and_exit(self) -> None:
            missing_by_tag = {
                key: self._missing_segments(key)
                for key in self._tab_order
            }
            missing_by_tag = {k: v for k, v in missing_by_tag.items() if v}

            if missing_by_tag:
                missing_text = self._format_missing(missing_by_tag)
                if prompt_missing_after_close:
                    choice = QtWidgets.QMessageBox.question(
                        self,
                        "Incomplete Bullseye Input",
                        (
                            "Missing segment values detected:\n"
                            f"{missing_text}\n\n"
                            "Fill all missing values with 0.0 before Save + Exit?"
                        ),
                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel,
                        QtWidgets.QMessageBox.Yes,
                    )
                    if choice == QtWidgets.QMessageBox.Cancel:
                        return
                    if choice == QtWidgets.QMessageBox.Yes:
                        for key in self._tab_order:
                            for seg_id in self._missing_segments(key):
                                self._state(key)["edit_values"][seg_id] = 0.0
                            self._state(key)["pending"] = True
                            self._refresh_bullseye_values(key)
                else:
                    choice = QtWidgets.QMessageBox.question(
                        self,
                        "Incomplete Bullseye Input",
                        (
                            "Missing segment values detected:\n"
                            f"{missing_text}\n\n"
                            "Save + Exit anyway with partial values?"
                        ),
                        QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                        QtWidgets.QMessageBox.No,
                    )
                    if choice != QtWidgets.QMessageBox.Yes:
                        return

            self._apply_all_pending_to_3d()
            ok = self._save_vtk_dialog()
            if ok:
                self._closing_after_save = True
                self.close()

        def closeEvent(self, event) -> None:  # noqa: N802
            if self._closing_after_save:
                event.accept()
                return

            missing_by_tag = {
                key: self._missing_segments(key)
                for key in self._tab_order
            }
            missing_by_tag = {k: v for k, v in missing_by_tag.items() if v}
            if not missing_by_tag:
                event.accept()
                return

            missing_text = self._format_missing(missing_by_tag)
            if prompt_missing_after_close:
                choice = QtWidgets.QMessageBox.question(
                    self,
                    "Incomplete Bullseye Input",
                    (
                        "You are closing with incomplete segments:\n"
                        f"{missing_text}\n\n"
                        "Fill missing with 0.0 and close?"
                    ),
                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Cancel,
                    QtWidgets.QMessageBox.Yes,
                )
                if choice == QtWidgets.QMessageBox.Cancel:
                    event.ignore()
                    return
                if choice == QtWidgets.QMessageBox.Yes:
                    for key in self._tab_order:
                        for seg_id in self._missing_segments(key):
                            self._state(key)["edit_values"][seg_id] = 0.0
                            self._state(key)["pending"] = True
                        self._refresh_bullseye_values(key)
                    self._apply_all_pending_to_3d()
                event.accept()
                return

            choice = QtWidgets.QMessageBox.question(
                self,
                "Incomplete Bullseye Input",
                (
                    "You are closing with incomplete segments:\n"
                    f"{missing_text}\n\n"
                    "Close anyway and keep partial values?"
                ),
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No,
            )
            if choice == QtWidgets.QMessageBox.Yes:
                event.accept()
            else:
                event.ignore()

    app = QtWidgets.QApplication.instance()
    owns_app = app is None
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    window = UnifiedWorkbenchWindow()
    window.show()
    app.exec()

    if owns_app:
        app.quit()
    return {
        "purkinje": dict(window._state("purkinje")["edit_values"]),
        "pmj": dict(window._state("pmj")["edit_values"]),
        "thickness_bundles": dict(window._state("thickness_bundles")["edit_values"]),
        "endocardial_surface_area_by_tag": dict(endocardial_area_by_segment),
        "endocardial_selector_source": str(endocardial_selector_source),
    }


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            "Unified GUI: interactive bullseye + mapped 3D in one window."
        )
    )
    ap.add_argument("--vtk", type=Path, required=True, help="Input VTK/VTU/VTP mesh path.")
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output VTK path with all mapped density fields.",
    )
    ap.add_argument(
        "--reference",
        default=None,
        help="Division reference name (for example: aha17, apical4_mid5_basal5).",
    )
    ap.add_argument(
        "--config-source",
        choices=["anatomical", "local"],
        default="anatomical",
        help="Where to load bullseye config from.",
    )
    ap.add_argument(
        "--config-file",
        type=Path,
        default=None,
        help="Optional path to a lv_division_config.py file.",
    )
    ap.add_argument(
        "--segment-field",
        default=None,
        help="Anatomical segment field name in VTK. If omitted, auto-detected.",
    )
    ap.add_argument(
        "--segment-association",
        choices=["auto", "point", "cell"],
        default="auto",
        help="Whether segment IDs are in point_data or cell_data.",
    )
    ap.add_argument(
        "--purkinje-field",
        default="purkinje_density",
        help="Output scalar field name for Purkinje density.",
    )
    ap.add_argument(
        "--pmj-field",
        default="pmj_per_mm2",
        help="Output scalar field name for PMJ density in /mm2.",
    )
    ap.add_argument(
        "--thickness-bundles-field",
        default="thickness_bundles",
        help="Output scalar field name for Thickness Bundles.",
    )
    ap.add_argument(
        "--endocardial-area-field",
        default="endocardial_surface_area_by_tag",
        help="Output scalar field name for automatic endocardial area-by-tag values.",
    )
    ap.add_argument(
        "--endocardial-selector-field",
        default=None,
        help=(
            "Optional point/cell field to select endocardial cells for area "
            "computation. If omitted, auto-detects endo_surface tags and "
            "falls back to all cells."
        ),
    )
    ap.add_argument(
        "--endocardial-selector-value",
        default="1",
        help="Selector value used with --endocardial-selector-field (default: 1).",
    )
    ap.add_argument(
        "--purkinje-max-value",
        type=float,
        default=1.0,
        help="Maximum per-segment value for Purkinje density (default: 1).",
    )
    ap.add_argument(
        "--pmj-max-value",
        type=float,
        default=1.0,
        help="Maximum per-segment value for PMJ /mm2 (default: 1).",
    )
    ap.add_argument(
        "--thickness-bundles-max-value",
        type=float,
        default=200.0,
        help="Maximum per-segment value for Thickness Bundles (default: 200).",
    )
    ap.add_argument(
        "--allow-partial",
        action="store_true",
        help="Allow closing GUI with missing segments without auto-filling 0.0.",
    )
    return ap


def main() -> None:
    args = build_arg_parser().parse_args()
    values = launch_unified_workbench(
        vtk_path=args.vtk,
        output_vtk_path=args.output,
        reference_name=args.reference,
        config_source=args.config_source,
        config_file=args.config_file,
        segment_field=args.segment_field,
        segment_association=args.segment_association,
        purkinje_field_name=args.purkinje_field,
        pmj_field_name=args.pmj_field,
        thickness_bundles_field_name=args.thickness_bundles_field,
        endocardial_area_field_name=args.endocardial_area_field,
        endocardial_selector_field=args.endocardial_selector_field,
        endocardial_selector_value=args.endocardial_selector_value,
        purkinje_max_value=args.purkinje_max_value,
        pmj_max_value=args.pmj_max_value,
        thickness_bundles_max_value=args.thickness_bundles_max_value,
        prompt_missing_after_close=not args.allow_partial,
    )

    print("\nSolution list (copy/edit for quick adjustments):")
    print(f"Purkinje density (0..{args.purkinje_max_value:g}):")
    for seg_id in sorted(values["purkinje"]):
        print(f"  {seg_id}: {values['purkinje'][seg_id]:.3f}")
    print(f"PMJ density (/mm2, 0..{args.pmj_max_value:g}):")
    for seg_id in sorted(values["pmj"]):
        print(f"  {seg_id}: {values['pmj'][seg_id]:.3f}")
    print(f"Thickness Bundles (0..{args.thickness_bundles_max_value:g}):")
    for seg_id in sorted(values["thickness_bundles"]):
        print(f"  {seg_id}: {values['thickness_bundles'][seg_id]:.3f}")
    print("Endocardial surface area by tag:")
    for seg_id in sorted(values["endocardial_surface_area_by_tag"]):
        print(f"  {seg_id}: {values['endocardial_surface_area_by_tag'][seg_id]:.6f}")
    print(f"Endocardial selector source: {values['endocardial_selector_source']}")

    solution_payload = {
        "reference": args.reference,
        "purkinje_max_value": args.purkinje_max_value,
        "pmj_max_value": args.pmj_max_value,
        "thickness_bundles_max_value": args.thickness_bundles_max_value,
        "endocardial_area_field": args.endocardial_area_field,
        "endocardial_selector_field": args.endocardial_selector_field,
        "endocardial_selector_value": args.endocardial_selector_value,
        "endocardial_selector_source": values["endocardial_selector_source"],
        "purkinje": {str(k): float(v) for k, v in sorted(values["purkinje"].items())},
        "pmj": {str(k): float(v) for k, v in sorted(values["pmj"].items())},
        "thickness_bundles": {
            str(k): float(v) for k, v in sorted(values["thickness_bundles"].items())
        },
        "endocardial_surface_area_by_tag": {
            str(k): float(v)
            for k, v in sorted(values["endocardial_surface_area_by_tag"].items())
        },
    }
    print("\nSolution JSON:")
    print(json.dumps(solution_payload, indent=2, sort_keys=False))

    if args.output is not None:
        solution_path = Path(str(args.output) + ".solution.json")
        solution_path.parent.mkdir(parents=True, exist_ok=True)
        with open(solution_path, "w", encoding="utf-8") as f:
            json.dump(solution_payload, f, indent=2)
            f.write("\n")
        print(f"\nSaved solution list: {solution_path}")


if __name__ == "__main__":
    main()
