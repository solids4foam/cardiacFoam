# Postprocessing Protocol Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the postprocessing protocol — typed contract, `TableWriter` utility, `plot_builder` merge, schema v1.1, per-tutorial table summaries, and spec wiring.

**Architecture:** Protocol + Shared Base Utilities. `PostprocessingProtocol` stub enforces the `run_postprocessing(**kwargs) -> list[dict]` contract. `TableWriter` handles the standard metadata envelope for tabular outputs. `plot_builder.py` (from `driverFOAMagentPreparation`) is merged into the main package.

**Tech Stack:** Python 3.10+, pandas, plotly, numpy, unittest (stdlib)

**Working directory for all test commands:** `cardiacFoamsacred/applications/scripts/driverFoam/`

---

## File Map

| Action | File |
|--------|------|
| Create | `openfoam_driver/postprocessing/plot_builder.py` |
| Create | `openfoam_driver/postprocessing/table_writer.py` |
| Create | `openfoam_driver/tests/test_table_writer.py` |
| Create | `tutorials/NiedererEtAl2012/setupNiedererEtAl2012/postProcessing/table_summary.py` |
| Create | `tutorials/singleCell/setupSingleCell/postProcessing/table_summary.py` |
| Create | `tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing/table_summary.py` |
| Modify | `openfoam_driver/postprocessing/__init__.py` |
| Modify | `openfoam_driver/postprocessing/driver.py` |
| Modify | `openfoam_driver/tests/test_postprocessing_driver.py` |
| Modify | `tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/post_processing_manufactured.py` |
| Modify | `tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing_restCurves.py` |
| Modify | `openfoam_driver/core/defaults/niederer_2012.py` |
| Modify | `openfoam_driver/core/defaults/single_cell.py` |
| Modify | `openfoam_driver/core/defaults/restitution_curves.py` |
| Modify | `openfoam_driver/specs/tutorials/niederer_2012.py` |
| Modify | `openfoam_driver/specs/tutorials/single_cell.py` |
| Modify | `openfoam_driver/specs/tutorials/restitution_curves.py` |

All paths are relative to `cardiacFoamsacred/applications/scripts/driverFoam/` unless prefixed with `tutorials/`.

---

## Task 1: Copy `plot_builder.py` and update `__init__.py`

**Files:**
- Create: `openfoam_driver/postprocessing/plot_builder.py`
- Modify: `openfoam_driver/postprocessing/__init__.py`

- [ ] **Step 1: Copy `plot_builder.py` from `driverFOAMagentPreparation`**

Copy the file verbatim:

```bash
cp /Users/simaocastro/cardiacFoamEPsacred/driverFOAMagentPreparation/openfoam_driver/postprocessing/plot_builder.py \
   /Users/simaocastro/cardiacFoamEPsacred/cardiacFoamsacred/applications/scripts/driverFoam/openfoam_driver/postprocessing/plot_builder.py
```

- [ ] **Step 2: Verify import smoke test**

Run from `cardiacFoamsacred/applications/scripts/driverFoam/`:
```bash
python3 -c "from openfoam_driver.postprocessing.plot_builder import GroupShadedColors, PlotSpec, TraceSpec, build_line_traces, load_csv_folder, make_toggle_button, write_plot; print('OK')"
```
Expected: `OK`

- [ ] **Step 3: Update `__init__.py` to export `plot_builder` symbols and add `PostprocessingProtocol`**

Replace the full contents of `openfoam_driver/postprocessing/__init__.py` with:

```python
"""Shared utilities for tutorial post-processing scripts."""
from __future__ import annotations

from typing import Protocol, runtime_checkable

from .driver import PostprocessTask, run_postprocess_tasks
from .plot_builder import (
    DEFAULT_PALETTE,
    GroupShadedColors,
    PlotSpec,
    TraceSpec,
    build_line_traces,
    load_csv_folder,
    make_toggle_button,
    write_plot,
)
from .style import (
    apply_plotly_layout,
    configure_matplotlib_defaults,
    finalize_matplotlib_figure,
    style_matplotlib_axes,
    write_plotly_html,
)


@runtime_checkable
class PostprocessingProtocol(Protocol):
    """Typing stub for tutorial post-processing entry points.

    Every tutorial postprocessing script must expose a function named
    ``run_postprocessing`` that satisfies this signature.  Scripts are not
    required to subclass this Protocol — IDEs and mypy will flag mismatches
    when type-checking is enabled.
    """

    def __call__(
        self,
        *,
        output_dir: str,
        setup_root: str | None = None,
        **kwargs: object,
    ) -> list[dict]: ...


__all__ = [
    # contract
    "PostprocessingProtocol",
    # driver
    "PostprocessTask",
    "run_postprocess_tasks",
    # plot_builder — declarative Plotly helpers
    "DEFAULT_PALETTE",
    "GroupShadedColors",
    "PlotSpec",
    "TraceSpec",
    "build_line_traces",
    "load_csv_folder",
    "make_toggle_button",
    "write_plot",
    # style
    "apply_plotly_layout",
    "configure_matplotlib_defaults",
    "finalize_matplotlib_figure",
    "style_matplotlib_axes",
    "write_plotly_html",
]
```

- [ ] **Step 4: Verify `PostprocessingProtocol` is importable**

```bash
python3 -c "from openfoam_driver.postprocessing import PostprocessingProtocol; print('OK')"
```
Expected: `OK`

- [ ] **Step 5: Run existing tests to confirm nothing broke**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver -v 2>&1 | tail -5
```
Expected: `Ran 5 tests` … `OK`

- [ ] **Step 6: Commit**

```bash
git add openfoam_driver/postprocessing/plot_builder.py openfoam_driver/postprocessing/__init__.py
git commit -m "feat: add plot_builder and PostprocessingProtocol to postprocessing package"
```

---

## Task 2: Bump `schema_version` to "1.1" in `driver.py`

**Files:**
- Modify: `openfoam_driver/postprocessing/driver.py` (line 117)
- Modify: `openfoam_driver/tests/test_postprocessing_driver.py` (line 53)

- [ ] **Step 1: Update the failing test first**

In `openfoam_driver/tests/test_postprocessing_driver.py`, change line 53:
```python
# Before
self.assertEqual(manifest["schema_version"], "1.0")
# After
self.assertEqual(manifest["schema_version"], "1.1")
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_writes_schema_version_and_artifacts -v 2>&1 | tail -5
```
Expected: `FAIL` — `'1.0' != '1.1'`

- [ ] **Step 3: Bump schema_version in `driver.py`**

In `openfoam_driver/postprocessing/driver.py`, change line 117:
```python
# Before
        "schema_version": "1.0",
# After
        "schema_version": "1.1",
```

- [ ] **Step 4: Run all postprocessing tests**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver -v 2>&1 | tail -5
```
Expected: `Ran 5 tests` … `OK`

- [ ] **Step 5: Commit**

```bash
git add openfoam_driver/postprocessing/driver.py openfoam_driver/tests/test_postprocessing_driver.py
git commit -m "feat: bump plots.json schema_version to 1.1"
```

---

## Task 3: Write `table_writer.py` with tests

**Files:**
- Create: `openfoam_driver/postprocessing/table_writer.py`
- Create: `openfoam_driver/tests/test_table_writer.py`

- [ ] **Step 1: Write the failing tests**

Create `openfoam_driver/tests/test_table_writer.py`:

```python
from __future__ import annotations

import json
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter


class TestTableWriter(unittest.TestCase):
    def test_writes_csv_with_envelope(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            rows = [
                {"case_id": "case_A", "DX_mm": 0.1, "activation_ms": 42.3},
                {"case_id": "case_B", "DX_mm": 0.2, "activation_ms": 55.1},
            ]
            meta = TableMetadata(
                tutorial="TestTutorial",
                units={"activation_ms": "ms", "DX_mm": "mm"},
            )
            artifacts = TableWriter.write(rows, output_dir, "test_summary", "Test label", meta)

            csv_path = output_dir / "test_summary.csv"
            self.assertTrue(csv_path.exists())
            text = csv_path.read_text()
            self.assertIn("# tutorial: TestTutorial", text)
            self.assertIn("# generated_at:", text)
            self.assertIn('"activation_ms": "ms"', text)
            self.assertIn("case_id,DX_mm,activation_ms", text)
            self.assertIn("case_A,0.1,42.3", text)
            self.assertIn("case_B,0.2,55.1", text)

    def test_writes_html_with_metadata_and_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            rows = [{"col_a": "x", "col_b": 1}]
            meta = TableMetadata(tutorial="HtmlTest", units={"col_b": "ms"})
            TableWriter.write(rows, output_dir, "html_test", "HTML label", meta)

            html_text = (output_dir / "html_test.html").read_text()
            self.assertIn("HtmlTest", html_text)
            self.assertIn("col_a", html_text)
            self.assertIn("<td>x</td>", html_text)

    def test_returns_two_artifacts_with_correct_schema(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            meta = TableMetadata(tutorial="ArtifactTest", units={})
            artifacts = TableWriter.write(
                [{"x": 1}], output_dir, "art_stem", "Art label", meta
            )
            self.assertEqual(len(artifacts), 2)
            by_format = {a["format"]: a for a in artifacts}
            self.assertIn("csv", by_format)
            self.assertIn("html", by_format)
            for a in artifacts:
                self.assertEqual(a["kind"], "table")
                self.assertEqual(a["label"], "Art label")
                self.assertIn("path", a)

    def test_handles_empty_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            meta = TableMetadata(tutorial="EmptyTest", units={})
            artifacts = TableWriter.write([], output_dir, "empty_stem", "Empty", meta)
            self.assertEqual(len(artifacts), 2)
            csv_text = (output_dir / "empty_stem.csv").read_text()
            self.assertIn("# tutorial: EmptyTest", csv_text)

    def test_metadata_autofills_generated_at(self) -> None:
        meta = TableMetadata(tutorial="T", units={})
        self.assertNotEqual(meta.generated_at, "")
        # Must be a parseable ISO-8601 string
        datetime.fromisoformat(meta.generated_at)

    def test_metadata_preserves_explicit_generated_at(self) -> None:
        meta = TableMetadata(tutorial="T", units={}, generated_at="2026-01-01T00:00:00+00:00")
        self.assertEqual(meta.generated_at, "2026-01-01T00:00:00+00:00")

    def test_artifact_paths_are_relative_filenames(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            meta = TableMetadata(tutorial="T", units={})
            artifacts = TableWriter.write([{"a": 1}], output_dir, "stem", "L", meta)
            for a in artifacts:
                # path must be just the filename, not absolute
                self.assertFalse(Path(a["path"]).is_absolute())
                self.assertIn("stem", a["path"])


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
python3 -m unittest openfoam_driver.tests.test_table_writer -v 2>&1 | tail -10
```
Expected: `ImportError` — `table_writer` does not exist yet

- [ ] **Step 3: Implement `table_writer.py`**

Create `openfoam_driver/postprocessing/table_writer.py`:

```python
"""table_writer.py — Standard tabular output with metadata envelope."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


@dataclass
class TableMetadata:
    """Metadata envelope attached to every tutorial table output.

    Attributes
    ----------
    tutorial:
        Human-readable tutorial name, e.g. ``"NiedererEtAl2012"``.
    units:
        Mapping of column name → unit string, e.g.
        ``{"activationTime": "ms", "DX": "mm"}``.
    generated_at:
        UTC ISO-8601 timestamp string.  Auto-filled on construction if empty.
    """

    tutorial: str
    units: dict[str, str] = field(default_factory=dict)
    generated_at: str = ""

    def __post_init__(self) -> None:
        if not self.generated_at:
            self.generated_at = datetime.now(timezone.utc).isoformat()


class TableWriter:
    """Write tutorial summary tables as CSV (with comment envelope) and HTML."""

    @staticmethod
    def write(
        rows: list[dict[str, Any]],
        output_dir: str | Path,
        filename_stem: str,
        label: str,
        metadata: TableMetadata,
    ) -> list[dict[str, Any]]:
        """Write *rows* as ``<filename_stem>.csv`` and ``<filename_stem>.html``.

        Parameters
        ----------
        rows:
            List of dicts where every dict has the same keys (tutorial-defined
            columns).  An empty list is allowed — only the envelope is written.
        output_dir:
            Directory into which the files are written.
        filename_stem:
            Base name without extension, e.g. ``"NiedererEtAl2012_summary"``.
        label:
            Human-readable description used in the artifact entry.
        metadata:
            :class:`TableMetadata` providing tutorial name, units, and timestamp.

        Returns
        -------
        Two artifact dicts (CSV and HTML), each with keys ``path``, ``label``,
        ``kind``, ``format``.  Paths are relative filenames (not absolute).
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        csv_path = output_dir / f"{filename_stem}.csv"
        html_path = output_dir / f"{filename_stem}.html"

        fieldnames: list[str] = list(rows[0].keys()) if rows else []

        # ------------------------------------------------------------------
        # CSV: comment-line envelope then data
        # ------------------------------------------------------------------
        lines: list[str] = [
            f"# tutorial: {metadata.tutorial}",
            f"# generated_at: {metadata.generated_at}",
            f"# units: {json.dumps(metadata.units)}",
        ]
        if fieldnames:
            lines.append(",".join(fieldnames))
            for row in rows:
                lines.append(",".join(str(row.get(f, "")) for f in fieldnames))
        csv_path.write_text("\n".join(lines) + "\n")

        # ------------------------------------------------------------------
        # HTML: styled table with metadata header block
        # ------------------------------------------------------------------
        parts: list[str] = [
            "<!DOCTYPE html><html><head>",
            "<style>",
            "body{font-family:Arial,sans-serif;margin:20px}",
            "table{border-collapse:collapse;width:100%}",
            "th,td{border:1px solid #ccc;padding:6px 10px;text-align:left}",
            "th{background:#f0f0f0}",
            ".meta{color:#555;margin-bottom:14px;font-size:13px;line-height:1.6}",
            "</style></head><body>",
            "<div class='meta'>",
            f"<strong>tutorial:</strong> {metadata.tutorial}&nbsp;&nbsp;",
            f"<strong>generated_at:</strong> {metadata.generated_at}<br>",
            f"<strong>units:</strong> {json.dumps(metadata.units)}",
            "</div>",
            "<table><thead><tr>",
        ]
        for f in fieldnames:
            parts.append(f"<th>{f}</th>")
        parts.append("</tr></thead><tbody>")
        for row in rows:
            parts.append("<tr>")
            for f in fieldnames:
                parts.append(f"<td>{row.get(f, '')}</td>")
            parts.append("</tr>")
        parts.append("</tbody></table></body></html>")
        html_path.write_text("".join(parts))

        return [
            {
                "path": csv_path.name,
                "label": label,
                "kind": "table",
                "format": "csv",
            },
            {
                "path": html_path.name,
                "label": label,
                "kind": "table",
                "format": "html",
            },
        ]
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
python3 -m unittest openfoam_driver.tests.test_table_writer -v 2>&1 | tail -10
```
Expected: `Ran 7 tests` … `OK`

- [ ] **Step 5: Run all tests to confirm no regressions**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver openfoam_driver.tests.test_table_writer -v 2>&1 | tail -5
```
Expected: `Ran 12 tests` … `OK`

- [ ] **Step 6: Commit**

```bash
git add openfoam_driver/postprocessing/table_writer.py openfoam_driver/tests/test_table_writer.py
git commit -m "feat: add TableWriter utility for standard tabular output with metadata envelope"
```

---

## Task 4: Fix `post_processing_manufactured.py` type annotation

**Files:**
- Modify: `tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/post_processing_manufactured.py` (line 2745)

- [ ] **Step 1: Fix the return type annotation**

In `post_processing_manufactured.py` at line 2745, change:
```python
# Before
def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> None:
# After
def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **_: object) -> list[dict]:
```

- [ ] **Step 2: Run existing postprocessing driver tests (uses this script)**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver -v 2>&1 | tail -5
```
Expected: `Ran 5 tests` … `OK`

- [ ] **Step 3: Commit**

```bash
git add ../../tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/post_processing_manufactured.py
git commit -m "fix: correct run_postprocessing return type annotation in manufactured FDA script"
```

---

## Task 5: Fix `postProcessing_restCurves.py` to return artifacts

**Files:**
- Modify: `tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing_restCurves.py` (lines 509–537)

- [ ] **Step 1: Write a failing test in `test_postprocessing_driver.py`**

Add this test to `openfoam_driver/tests/test_postprocessing_driver.py` inside `TestPostprocessingDriver`:

```python
def test_restitution_run_postprocessing_returns_list(self) -> None:
    """run_postprocessing must return list[dict], not None."""
    repo_root = _repo_root_from_test()
    module_path = (
        repo_root
        / "tutorials"
        / "restitutionCurves_s1s2Protocol"
        / "setupRestitutionCurves_s1s2Protocol"
        / "postProcessing_restCurves.py"
    )
    import importlib.util

    spec = importlib.util.spec_from_file_location("restcurves", module_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    fn = mod.run_postprocessing
    import inspect

    hints = fn.__annotations__
    self.assertEqual(
        hints.get("return", type(None)),
        list,
        "run_postprocessing must annotate return as list[dict], not None",
    )
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_restitution_run_postprocessing_returns_list -v 2>&1 | tail -5
```
Expected: `FAIL`

- [ ] **Step 3: Update `run_postprocessing` in `postProcessing_restCurves.py`**

Replace lines 509–537 (the `run_postprocessing` function) with:

```python
def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    **kwargs,
) -> list[dict]:
    """PostprocessTask-compatible entry point for the openfoam_driver framework.

    Called by the driver engine after all S1–S2 simulations complete.
    Expected kwargs:
        ionic_models  (list[str])        - Models to post-process.
        tissue_map    (dict[str, list])  - Tissue types per ionic model.
        show_plots    (bool)             - Whether to display plots interactively.
    """
    ionic_models: list[str] = kwargs.get("ionic_models", [])
    tissue_map: dict[str, list[str]] = kwargs.get("tissue_map", {})
    show_plots: bool = kwargs.get("show_plots", False)

    output_dir_path = Path(output_dir)
    artifacts: list[dict] = []

    for model in ionic_models:
        tissues = list(tissue_map.get(model, []))
        postprocess_one_ionic_model(
            base_dir=output_dir_path.parent,
            output_folder=output_dir_path.name,
            ionic_model=model,
            tissues=tissues,
            show_plot=show_plots,
        )
        fig_path = output_dir_path / f"{model}_restitution.png"
        csv_path = output_dir_path / f"{model}_restitution.csv"
        if fig_path.exists():
            artifacts.append({
                "path": fig_path.name,
                "label": f"{model} restitution curve",
                "kind": "plot",
                "format": "png",
            })
        if csv_path.exists():
            artifacts.append({
                "path": csv_path.name,
                "label": f"{model} restitution data",
                "kind": "data",
                "format": "csv",
            })

    return artifacts
```

- [ ] **Step 4: Run all postprocessing tests**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver -v 2>&1 | tail -5
```
Expected: `Ran 6 tests` … `OK`

- [ ] **Step 5: Commit**

```bash
git add ../../tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing_restCurves.py \
        openfoam_driver/tests/test_postprocessing_driver.py
git commit -m "fix: make run_postprocessing in restCurves return artifacts list"
```

---

## Task 6: Write `NiedererEtAl2012/postProcessing/table_summary.py`

**Files:**
- Create: `tutorials/NiedererEtAl2012/setupNiedererEtAl2012/postProcessing/table_summary.py`

- [ ] **Step 1: Write the failing test**

Add to `openfoam_driver/tests/test_postprocessing_driver.py` inside `TestPostprocessingDriver`:

```python
def test_niederer_table_summary_produces_csv_and_html(self) -> None:
    repo_root = _repo_root_from_test()
    setup_root = (
        repo_root / "tutorials" / "NiedererEtAl2012" / "setupNiedererEtAl2012"
    )

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        # Write a fake points CSV matching the expected pattern
        csv_content = "activationTime,Points:0,Points:1,Points:2\n0.042,0.0,0.0,0.007\n0.055,0.02,0.003,0.007\n"
        (output_dir / "implicit_TNNP_epicardialCells_points_DT001_DX01.csv").write_text(csv_content)

        run_postprocess_tasks(
            setup_root=setup_root,
            output_dir=output_dir,
            tutorial_name="NiedererEtAl2012",
            tasks=[
                PostprocessTask(
                    module_relpath=Path("postProcessing/table_summary.py")
                )
            ],
        )

        self.assertTrue((output_dir / "NiedererEtAl2012_summary.csv").exists())
        self.assertTrue((output_dir / "NiedererEtAl2012_summary.html").exists())
        csv_text = (output_dir / "NiedererEtAl2012_summary.csv").read_text()
        self.assertIn("# tutorial: NiedererEtAl2012", csv_text)
        self.assertIn("case_id", csv_text)
        self.assertIn("implicit_TNNP_epicardialCells", csv_text)
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_niederer_table_summary_produces_csv_and_html -v 2>&1 | tail -5
```
Expected: `ERROR` — module `table_summary` not found

- [ ] **Step 3: Create the `postProcessing/` subfolder and write `table_summary.py`**

The file already exists at `tutorials/NiedererEtAl2012/setupNiedererEtAl2012/postProcessing/` (the folder exists). Create `table_summary.py`:

```python
"""table_summary.py — Activation time summary table for NiedererEtAl2012.

Reads *points_DT*_DX*.csv files from output_dir, extracts activation times
per probe point, and writes NiedererEtAl2012_summary.csv + .html via TableWriter.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

TUTORIALS_ROOT = Path(__file__).resolve().parents[4]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.plotting_common import extract_dx_dt
from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter


def _parse_filename(filename: str) -> tuple[str, float, float, str]:
    """Return (case_id, dx_mm, dt_ms, solver) from a points CSV filename.

    Expected pattern: ``{solver}_{model}_{tissue}_points_DT{tag}_DX{tag}.csv``
    """
    stem = filename.replace(".csv", "")
    dx, dt = extract_dx_dt(stem)
    case_id = stem.split("_points_DT")[0]
    solver = case_id.split("_")[0]
    return case_id, round(dx, 4), round(dt, 5), solver


def build_summary_rows(output_dir: Path) -> list[dict]:
    files = sorted(output_dir.glob("*points_DT*_DX*.csv"))
    if not files:
        print(f"[NiedererEtAl2012/table_summary] No points CSV files in {output_dir}")
        return []

    rows = []
    for fpath in files:
        case_id, dx, dt, solver = _parse_filename(fpath.name)
        df = pd.read_csv(fpath)
        activation_ms = df["activationTime"].values * 1000.0  # s → ms
        row: dict = {
            "case_id": case_id,
            "DX_mm": dx,
            "DT_ms": dt,
            "solver": solver,
        }
        for i, t in enumerate(activation_ms):
            row[f"point_{i}_activation_ms"] = round(float(t), 4)
        rows.append(row)

    return rows


def run_postprocessing(
    *, output_dir: str, setup_root: str | None = None, **_: object
) -> list[dict]:
    output_path = Path(output_dir)
    rows = build_summary_rows(output_path)
    if not rows:
        return []
    meta = TableMetadata(
        tutorial="NiedererEtAl2012",
        units={"activationTime": "ms", "DX": "mm", "DT": "ms"},
    )
    return TableWriter.write(
        rows,
        output_path,
        "NiedererEtAl2012_summary",
        "Niederer activation time summary",
        meta,
    )


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[1] / "NiedererFoam"
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
```

- [ ] **Step 4: Run the test**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_niederer_table_summary_produces_csv_and_html -v 2>&1 | tail -5
```
Expected: `OK`

- [ ] **Step 5: Run all tests**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver openfoam_driver.tests.test_table_writer -v 2>&1 | tail -5
```
Expected: `Ran 14 tests` … `OK`

- [ ] **Step 6: Commit**

```bash
git add ../../tutorials/NiedererEtAl2012/setupNiedererEtAl2012/postProcessing/table_summary.py \
        openfoam_driver/tests/test_postprocessing_driver.py
git commit -m "feat: add Niederer activation time table_summary.py"
```

---

## Task 7: Write `singleCell/postProcessing/table_summary.py`

**Files:**
- Create: `tutorials/singleCell/setupSingleCell/postProcessing/table_summary.py` (new subfolder)

- [ ] **Step 1: Write the failing test**

Add to `TestPostprocessingDriver` in `test_postprocessing_driver.py`:

```python
def test_singlecell_table_summary_produces_csv_and_html(self) -> None:
    repo_root = _repo_root_from_test()
    setup_root = repo_root / "tutorials" / "singleCell" / "setupSingleCell"

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        # Fake single-cell .txt output: time Vm (space-separated)
        # Resting ~-85 mV, peak ~40 mV, repolarises back to ~-85 mV
        import numpy as np
        t = np.linspace(0, 0.5, 500)
        vm = np.full_like(t, -85.0)
        vm[50:150] = np.linspace(-85, 40, 100)   # upstroke
        vm[150:350] = np.linspace(40, -85, 200)  # repolarisation
        txt_lines = ["time Vm"] + [f"{ti:.4f} {vi:.4f}" for ti, vi in zip(t, vm)]
        (output_dir / "TNNP_epicardialCells_run.txt").write_text("\n".join(txt_lines))

        run_postprocess_tasks(
            setup_root=setup_root,
            output_dir=output_dir,
            tutorial_name="singleCell",
            tasks=[
                PostprocessTask(
                    module_relpath=Path("postProcessing/table_summary.py")
                )
            ],
        )

        self.assertTrue((output_dir / "singleCell_summary.csv").exists())
        self.assertTrue((output_dir / "singleCell_summary.html").exists())
        csv_text = (output_dir / "singleCell_summary.csv").read_text()
        self.assertIn("# tutorial: singleCell", csv_text)
        self.assertIn("APD_ms", csv_text)
        self.assertIn("peak_voltage_mV", csv_text)
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_singlecell_table_summary_produces_csv_and_html -v 2>&1 | tail -5
```
Expected: `ERROR` — module `table_summary` not found

- [ ] **Step 3: Create `postProcessing/` subfolder and write `table_summary.py`**

```bash
mkdir -p /Users/simaocastro/cardiacFoamEPsacred/cardiacFoamsacred/tutorials/singleCell/setupSingleCell/postProcessing
```

Create `tutorials/singleCell/setupSingleCell/postProcessing/table_summary.py`:

```python
"""table_summary.py — Voltage and APD summary table for singleCell tutorial.

Reads .txt simulation output files from output_dir.  Each file contains a
space-separated time series with columns: time, Vm, [additional state vars...].
Extracts resting voltage, peak voltage, and APD at 90% repolarisation.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

TUTORIALS_ROOT = Path(__file__).resolve().parents[4]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.plotting_common import parse_model_and_cell
from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter

# 10% of excursion above resting remaining = 90% repolarisation
_APD_REPOL_FRACTION = 0.10


def _compute_apd90(time: np.ndarray, vm: np.ndarray) -> float | None:
    """Return APD90 in ms, or None if detection fails."""
    resting = float(vm[0])
    peak_idx = int(np.argmax(vm))
    peak = float(vm[peak_idx])
    if peak <= resting:
        return None
    threshold = resting + _APD_REPOL_FRACTION * (peak - resting)
    for i in range(peak_idx + 1, len(vm)):
        if vm[i] <= threshold:
            # Linear interpolation for sub-sample accuracy
            frac = (threshold - float(vm[i - 1])) / (float(vm[i]) - float(vm[i - 1]))
            t_repol = float(time[i - 1]) + frac * (float(time[i]) - float(time[i - 1]))
            return (t_repol - float(time[peak_idx])) * 1000.0  # s → ms
    return None


def _extract_row(fpath: Path) -> dict | None:
    try:
        df = pd.read_csv(fpath, sep=r"\s+", comment="#", engine="python")
    except Exception as exc:
        print(f"[singleCell/table_summary] Could not read {fpath.name}: {exc}")
        return None
    if df.empty or len(df.columns) < 2:
        return None

    time = df.iloc[:, 0].values.astype(float)
    vm = df.iloc[:, 1].values.astype(float)
    apd = _compute_apd90(time, vm)
    model, cell = parse_model_and_cell(fpath.name)

    return {
        "case_id": fpath.stem,
        "ionic_model": model or "unknown",
        "tissue": cell or "unknown",
        "APD_ms": round(apd, 3) if apd is not None else "N/A",
        "peak_voltage_mV": round(float(np.max(vm)), 3),
        "resting_voltage_mV": round(float(vm[0]), 3),
    }


def run_postprocessing(
    *, output_dir: str, setup_root: str | None = None, **_: object
) -> list[dict]:
    output_path = Path(output_dir)
    txt_files = sorted(output_path.glob("*.txt"))
    if not txt_files:
        print(f"[singleCell/table_summary] No .txt files found in {output_path}")
        return []

    rows = [r for f in txt_files if (r := _extract_row(f)) is not None]
    if not rows:
        return []

    meta = TableMetadata(
        tutorial="singleCell",
        units={"APD_ms": "ms", "peak_voltage_mV": "mV", "resting_voltage_mV": "mV"},
    )
    return TableWriter.write(
        rows,
        output_path,
        "singleCell_summary",
        "Single cell voltage and APD summary",
        meta,
    )


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[2]
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
```

- [ ] **Step 4: Run the test**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_singlecell_table_summary_produces_csv_and_html -v 2>&1 | tail -5
```
Expected: `OK`

- [ ] **Step 5: Run all tests**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver openfoam_driver.tests.test_table_writer -v 2>&1 | tail -5
```
Expected: `Ran 15 tests` … `OK`

- [ ] **Step 6: Commit**

```bash
git add ../../tutorials/singleCell/setupSingleCell/postProcessing/table_summary.py \
        openfoam_driver/tests/test_postprocessing_driver.py
git commit -m "feat: add singleCell APD/voltage table_summary.py"
```

---

## Task 8: Write `restitutionCurves/postProcessing/table_summary.py`

**Files:**
- Create: `tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing/table_summary.py` (new subfolder)

- [ ] **Step 1: Write the failing test**

Add to `TestPostprocessingDriver` in `test_postprocessing_driver.py`:

```python
def test_restitution_table_summary_consolidates_model_csvs(self) -> None:
    repo_root = _repo_root_from_test()
    setup_root = (
        repo_root
        / "tutorials"
        / "restitutionCurves_s1s2Protocol"
        / "setupRestitutionCurves_s1s2Protocol"
    )

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)
        # Write fake per-model restitution CSVs (written by postProcessing_restCurves.py)
        (output_dir / "TNNP_restitution.csv").write_text(
            "tissue,DI_ms,APD90_ms\nepicardiaCells,300.0,280.0\nmCells,250.0,240.0\n"
        )
        (output_dir / "BuenoOrovio_restitution.csv").write_text(
            "tissue,DI_ms,APD90_ms\nepicardiaCells,310.0,290.0\n"
        )

        run_postprocess_tasks(
            setup_root=setup_root,
            output_dir=output_dir,
            tutorial_name="restitutionCurves_s1s2Protocol",
            tasks=[
                PostprocessTask(
                    module_relpath=Path("postProcessing/table_summary.py")
                )
            ],
        )

        self.assertTrue((output_dir / "restitutionCurves_summary.csv").exists())
        self.assertTrue((output_dir / "restitutionCurves_summary.html").exists())
        csv_text = (output_dir / "restitutionCurves_summary.csv").read_text()
        self.assertIn("# tutorial: restitutionCurves_s1s2Protocol", csv_text)
        self.assertIn("ionic_model", csv_text)
        self.assertIn("TNNP", csv_text)
        self.assertIn("BuenoOrovio", csv_text)
```

- [ ] **Step 2: Run test to confirm it fails**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_restitution_table_summary_consolidates_model_csvs -v 2>&1 | tail -5
```
Expected: `ERROR` — module `table_summary` not found

- [ ] **Step 3: Create `postProcessing/` subfolder and write `table_summary.py`**

```bash
mkdir -p /Users/simaocastro/cardiacFoamEPsacred/cardiacFoamsacred/tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing
```

Create `tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing/table_summary.py`:

```python
"""table_summary.py — Consolidated restitution table for restitutionCurves tutorial.

Reads per-model ``*_restitution.csv`` files produced by postProcessing_restCurves.py,
merges them into a single table with an ``ionic_model`` column, and writes
restitutionCurves_summary.csv + .html via TableWriter.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

TUTORIALS_ROOT = Path(__file__).resolve().parents[4]
if str(TUTORIALS_ROOT) not in sys.path:
    sys.path.insert(0, str(TUTORIALS_ROOT))

from openfoam_driver.postprocessing.table_writer import TableMetadata, TableWriter


def run_postprocessing(
    *, output_dir: str, setup_root: str | None = None, **_: object
) -> list[dict]:
    output_path = Path(output_dir)
    csv_files = sorted(output_path.glob("*_restitution.csv"))
    if not csv_files:
        print(f"[restitutionCurves/table_summary] No *_restitution.csv in {output_path}")
        return []

    rows: list[dict] = []
    for fpath in csv_files:
        ionic_model = fpath.stem.replace("_restitution", "")
        try:
            df = pd.read_csv(fpath)
        except Exception as exc:
            print(f"[restitutionCurves/table_summary] Could not read {fpath.name}: {exc}")
            continue
        for _, row in df.iterrows():
            rows.append({
                "ionic_model": ionic_model,
                "tissue": str(row.get("tissue", "unknown")),
                "DI_ms": round(float(row.get("DI_ms", 0.0)), 3),
                "APD_ms": round(float(row.get("APD90_ms", 0.0)), 3),
            })

    if not rows:
        return []

    meta = TableMetadata(
        tutorial="restitutionCurves_s1s2Protocol",
        units={"DI_ms": "ms", "APD_ms": "ms"},
    )
    return TableWriter.write(
        rows,
        output_path,
        "restitutionCurves_summary",
        "Restitution curves consolidated summary",
        meta,
    )


if __name__ == "__main__":
    folder = Path(__file__).resolve().parents[2]
    print(f"[table_summary] Default folder = {folder}")
    run_postprocessing(output_dir=str(folder))
```

- [ ] **Step 4: Run the test**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver.TestPostprocessingDriver.test_restitution_table_summary_consolidates_model_csvs -v 2>&1 | tail -5
```
Expected: `OK`

- [ ] **Step 5: Run all tests**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver openfoam_driver.tests.test_table_writer -v 2>&1 | tail -5
```
Expected: `Ran 16 tests` … `OK`

- [ ] **Step 6: Commit**

```bash
git add ../../tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing/table_summary.py \
        openfoam_driver/tests/test_postprocessing_driver.py
git commit -m "feat: add restitutionCurves consolidated table_summary.py"
```

---

## Task 9: Wire `table_summary` into tutorial specs via defaults

**Files:**
- Modify: `openfoam_driver/core/defaults/niederer_2012.py`
- Modify: `openfoam_driver/core/defaults/single_cell.py`
- Modify: `openfoam_driver/core/defaults/restitution_curves.py`
- Modify: `openfoam_driver/specs/tutorials/niederer_2012.py`
- Modify: `openfoam_driver/specs/tutorials/single_cell.py`
- Modify: `openfoam_driver/specs/tutorials/restitution_curves.py`

### 9a — Add `TABLE_SUMMARY_RELPATH` to defaults

- [ ] **Step 1: Update `core/defaults/niederer_2012.py`**

After line 55 (`CACHE_POSTPROCESS_FUNCTION = "run_postprocessing"`), add:
```python
TABLE_SUMMARY_RELPATH = Path("postProcessing/table_summary.py")
```

In the `__all__` list (around line 62), add `"TABLE_SUMMARY_RELPATH"` after the existing `POSTPROCESS` entries.

- [ ] **Step 2: Update `core/defaults/single_cell.py`**

After line 30 (`POSTPROCESS_FUNCTION_NAME = "run_postprocessing"`), add:
```python
TABLE_SUMMARY_RELPATH = Path("postProcessing/table_summary.py")
```

- [ ] **Step 3: Update `core/defaults/restitution_curves.py`**

After line 49 (`POSTPROCESS_FUNCTION_NAME = "run_postprocessing"`), add:
```python
TABLE_SUMMARY_RELPATH = Path("postProcessing/table_summary.py")
```

In the `__all__` list, add `"TABLE_SUMMARY_RELPATH"`.

- [ ] **Step 4: Verify defaults import**

```bash
python3 -c "
from openfoam_driver.core.defaults import niederer_2012, single_cell, restitution_curves
print(niederer_2012.TABLE_SUMMARY_RELPATH)
print(single_cell.TABLE_SUMMARY_RELPATH)
print(restitution_curves.TABLE_SUMMARY_RELPATH)
"
```
Expected:
```
postProcessing/table_summary.py
postProcessing/table_summary.py
postProcessing/table_summary.py
```

### 9b — Update `niederer_2012.py` spec

- [ ] **Step 5: Add `table_summary_relpath` to `_postprocess()` in `specs/tutorials/niederer_2012.py`**

In `_postprocess()` (starting at line 430), add parameter after `case_postprocess_cache_dirname`:
```python
    table_summary_relpath: Path = defaults.TABLE_SUMMARY_RELPATH,
```

Add a new `PostprocessTask` at the end of the `tasks` list (after the `points_postprocess_relpath` task):
```python
            PostprocessTask(
                module_relpath=table_summary_relpath,
            ),
```

- [ ] **Step 6: Add `table_summary_relpath` to `make_spec()` in `specs/tutorials/niederer_2012.py`**

In `make_spec()` (around line 499–507), add after `case_postprocess_cache_dirname`:
```python
    table_summary_relpath: str | Path = defaults.TABLE_SUMMARY_RELPATH,
```

In the body of `make_spec()`, after `cache_postprocess_path = Path(cache_postprocess_relpath)`, add:
```python
    table_summary_path = Path(table_summary_relpath)
```

In the `partial(_postprocess, ...)` call (around line 595), add:
```python
            table_summary_relpath=table_summary_path,
```

### 9c — Update `single_cell.py` spec

- [ ] **Step 7: Add `table_summary_relpath` to `_postprocess()` in `specs/tutorials/single_cell.py`**

In `_postprocess()` (line 93), add parameter after `strict_artifacts`:
```python
    table_summary_relpath: Path = defaults.TABLE_SUMMARY_RELPATH,
```

Add a new `PostprocessTask` in the `tasks` list after the existing postprocess task:
```python
            PostprocessTask(
                module_relpath=table_summary_relpath,
            ),
```

- [ ] **Step 8: Add `table_summary_relpath` to `make_spec()` in `specs/tutorials/single_cell.py`**

In `make_spec()` (around line 131), add after `postprocess_function_name`:
```python
    table_summary_relpath: str | Path = defaults.TABLE_SUMMARY_RELPATH,
```

In the body of `make_spec()`, after resolving `postprocess_script_relpath`, add:
```python
    table_summary_path = Path(table_summary_relpath)
```

Pass `table_summary_relpath=table_summary_path` to `partial(_postprocess, ...)`.

### 9d — Update `restitution_curves.py` spec

- [ ] **Step 9: Add `table_summary_relpath` to `_postprocess()` in `specs/tutorials/restitution_curves.py`**

In `_postprocess()` (line 160), add parameter after `strict_artifacts`:
```python
    table_summary_relpath: Path = defaults.TABLE_SUMMARY_RELPATH,
```

Add a new `PostprocessTask` after the existing one:
```python
            PostprocessTask(
                module_relpath=table_summary_relpath,
            ),
```

- [ ] **Step 10: Add `table_summary_relpath` to `make_spec()` in `specs/tutorials/restitution_curves.py`**

In `make_spec()` (around line 212), add after `postprocess_function_name`:
```python
    table_summary_relpath: str | Path = defaults.TABLE_SUMMARY_RELPATH,
```

In body of `make_spec()`, add:
```python
    table_summary_path = Path(table_summary_relpath)
```

Pass `table_summary_relpath=table_summary_path` to `partial(_postprocess, ...)`.

### 9e — Verify

- [ ] **Step 11: Verify specs import cleanly**

```bash
python3 -c "
from openfoam_driver.specs.tutorials.niederer_2012 import make_spec as n
from openfoam_driver.specs.tutorials.single_cell import make_spec as s
from openfoam_driver.specs.tutorials.restitution_curves import make_spec as r
print('All specs import OK')
"
```
Expected: `All specs import OK`

- [ ] **Step 12: Run all tests**

```bash
python3 -m unittest openfoam_driver.tests.test_postprocessing_driver openfoam_driver.tests.test_table_writer -v 2>&1 | tail -5
```
Expected: `Ran 16 tests` … `OK`

- [ ] **Step 13: Commit**

```bash
git add openfoam_driver/core/defaults/niederer_2012.py \
        openfoam_driver/core/defaults/single_cell.py \
        openfoam_driver/core/defaults/restitution_curves.py \
        openfoam_driver/specs/tutorials/niederer_2012.py \
        openfoam_driver/specs/tutorials/single_cell.py \
        openfoam_driver/specs/tutorials/restitution_curves.py
git commit -m "feat: wire table_summary PostprocessTask into NiedererEtAl2012, singleCell, restitutionCurves specs"
```
