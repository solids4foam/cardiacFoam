# Postprocessing Protocol Design

**Date:** 2026-04-08
**Status:** Approved
**Scope:** `applications/scripts/driverFoam/openfoam_driver/postprocessing/` + all tutorial `setup<Tutorial>/postProcessing/` scripts

---

## Problem

Each tutorial (`NiedererEtAl2012`, `singleCell`, `manufacturedFDA`, `restitutionCurves`) has ad-hoc postprocessing scripts that return inconsistent artifact shapes, lack tabular summary outputs, and do not integrate cleanly with the `plots.json` manifest consumed by the downstream app. A prior draft in `driverFOAMagentPreparation/` introduced good patterns (`plot_builder.py`, updated `driver.py`) but was never merged.

---

## Goals

1. A typed protocol/contract that every tutorial postprocessing script must conform to
2. Consistent `plots.json` manifest (v1.1) readable by the existing app
3. Structured tabular summaries (tutorial-specific columns, standard metadata envelope) exported as `.csv` + `.html`
4. Merge `driverFOAMagentPreparation` improvements into the main package
5. Minimal migration burden — existing plotting logic is not rewritten

## Non-Goals

- Static PNG/PDF export (handled by the downstream app via kaleido)
- Mandatory use of `plot_builder` for custom visualizations (3D subplots, complex updatemenus)
- A heart-simulation-level table schema (future work, user-defined)

---

## Approach

**Protocol + Shared Base Utilities.** The contract is enforced through a typed `PostprocessingProtocol` stub and a new `TableWriter` utility. Tutorial scripts conform by signature and return schema — no mandatory base class inheritance.

---

## Section 1 — The Contract

### Entry point

Every tutorial postprocessing script exposes exactly one public function:

```python
def run_postprocessing(*, output_dir: str, setup_root: str | None = None, **kwargs) -> list[dict]:
    ...
```

**Rules:**
- All parameters are keyword-only (`*` in signature)
- `output_dir` — always injected by the driver; directory where simulation outputs live and where artifacts are written
- `setup_root` — optional; path to `setup<Tutorial>/` folder, used when scripts reference static assets (e.g. Niederer reference Excel file)
- `**kwargs` — absorbs future driver-injected context without breaking existing scripts
- Return value is always a `list[dict]` — empty list `[]` when nothing to report, never `None`

### Artifact dict schema

Each item in the returned list must conform to:

```python
{
    "path":   str,   # file path relative to output_dir
    "label":  str,   # human-readable title shown in the app
    "kind":   Literal["plot", "table", "data", "report"],
    "format": Literal["html", "png", "csv", "json", "md", "dir"],
}
```

`kind` semantics:
- `"plot"` — interactive or static visualization
- `"table"` — structured tabular data (key metrics, summaries)
- `"data"` — raw simulation output (CSV samples, VTK files)
- `"report"` — human-readable narrative (markdown, PDF)

### Protocol stub

A `PostprocessingProtocol` typing stub is defined in `openfoam_driver/postprocessing/__init__.py`:

```python
from typing import Protocol, runtime_checkable

@runtime_checkable
class PostprocessingProtocol(Protocol):
    def __call__(
        self,
        *,
        output_dir: str,
        setup_root: str | None = None,
        **kwargs: object,
    ) -> list[dict]: ...
```

Scripts are not required to subclass it. IDEs and `mypy` will flag signature mismatches when type-checking is enabled.

---

## Section 2 — Shared Utilities

All modules live in `applications/scripts/driverFoam/openfoam_driver/postprocessing/`.

### `table_writer.py` (new)

Handles the standard table envelope so tutorials do not re-implement it.

```python
@dataclass
class TableMetadata:
    tutorial: str
    units: dict[str, str]          # e.g. {"activationTime": "ms", "DX": "mm"}
    generated_at: str = ""         # auto-filled with UTC ISO-8601 if empty

class TableWriter:
    @staticmethod
    def write(
        rows: list[dict],
        output_dir: str | Path,
        filename_stem: str,
        label: str,
        metadata: TableMetadata,
    ) -> list[dict]:               # returns artifact dicts (csv + html)
        ...
```

**CSV output** — envelope injected as comment lines before the header:

```
# tutorial: NiedererEtAl2012
# generated_at: 2026-04-08T12:00:00Z
# units: {"activationTime": "ms", "DX": "mm", "DT": "ms"}
case_id,DX,DT,point_0_ms,point_1_ms,...
implicit_TNNP_epi_DT001_DX01,0.1,0.01,42.3,55.1,...
```

**HTML output** — styled table with metadata displayed in a header block above the data rows.

Returns two artifact dicts:
```python
[
    {"path": "NiedererEtAl2012_summary.csv", "label": label, "kind": "table", "format": "csv"},
    {"path": "NiedererEtAl2012_summary.html", "label": label, "kind": "table", "format": "html"},
]
```

### `plot_builder.py` (from `driverFOAMagentPreparation`, merged in)

Declarative Plotly helpers — `PlotSpec`, `TraceSpec`, `GroupShadedColors`, `build_line_traces`, `write_plot`, `make_toggle_button`. Recommended for standard line/scatter plots. Not required for custom visualizations (3D subplots, complex updatemenus remain imperative).

### `plotting_common.py` + `style.py` (existing, replaced with `driverFOAMagentPreparation` versions)

`extract_dx_dt`, `lighten_hex_color`, `parse_model_and_cell`, `ordered_unique`, `build_visibility_mask`, `rename_cardiacfoam_trace`, `apply_plotly_layout`, `write_plotly_html` — no API changes, improved implementations from the draft.

---

## Section 3 — Driver Integration

### `driver.py` (replaced with `driverFOAMagentPreparation` version)

`PostprocessTask` dataclass and `run_postprocess_tasks()` are unchanged from the draft. The driver:

1. Loads each task module dynamically from `module_relpath`
2. Resolves `$OUTPUT_DIR` and `$SETUP_ROOT` placeholders in `kwargs`
3. Calls `run_postprocessing(**resolved_kwargs)`
4. Normalizes returned artifact dicts (adds `absolute_path`, `exists`, `task_index`, `module`, `function`)
5. Writes `plots.json` to `output_dir`

### `plots.json` schema — v1.1

Only change from v1.0: `schema_version` bumped to `"1.1"`, `kind: "table"` is now a valid value. The app's existing reader requires no breaking changes — new `kind` values are additive.

```json
{
    "schema_version": "1.1",
    "tutorial": "NiedererEtAl2012",
    "output_dir": "/path/to/output",
    "generated_at_utc": "2026-04-08T12:00:00Z",
    "artifact_count": 5,
    "artifacts": [
        {
            "path": "Niederer_vs_cardiacFoam.html",
            "absolute_path": "/full/path/Niederer_vs_cardiacFoam.html",
            "exists": true,
            "kind": "plot",
            "format": "html",
            "label": "Niederer vs cardiacFoam activation times",
            "task_index": 1,
            "module": "line_postProcessing.py",
            "function": "run_postprocessing"
        },
        {
            "path": "NiedererEtAl2012_summary.csv",
            "absolute_path": "/full/path/NiedererEtAl2012_summary.csv",
            "exists": true,
            "kind": "table",
            "format": "csv",
            "label": "Activation time summary table",
            "task_index": 3,
            "module": "table_summary.py",
            "function": "run_postprocessing"
        }
    ]
}
```

### `PostprocessTask` — unchanged

```python
@dataclass
class PostprocessTask:
    module_relpath: str
    function_name: str = "run_postprocessing"
    kwargs: dict = field(default_factory=dict)
```

### Tutorial spec wiring pattern

```python
tasks = [
    PostprocessTask("postProcessing/cache_postProcessing.py"),
    PostprocessTask(
        "postProcessing/line_postProcessing.py",
        kwargs={"excel_path": "$SETUP_ROOT/postProcessing/Niederer_graphs.../WebPlotDigitilizerdata.xlsx"},
    ),
    PostprocessTask("postProcessing/points_postProcessing.py"),
    PostprocessTask("postProcessing/table_summary.py"),   # new per tutorial
]
```

---

## Section 4 — Per-Tutorial Migration

### Existing scripts — minimal changes only

| Script | Required change |
|--------|----------------|
| `cache_postProcessing.py` | Verify returned dict has `kind` and `format` fields |
| `line_postProcessing.py` | Add `"kind": "plot", "format": "html"` to each returned dict |
| `points_postProcessing.py` | Add `"kind": "plot", "format": "html"` to returned dict |
| `singleCellinteractivePlots.py` | Add `"kind": "plot", "format": "html"` to returned dict |
| `post_processing_manufactured.py` | Add `kind`/`format` to all returned dicts |
| `postProcessing_restCurves.py` | Add `kind`/`format` to returned dicts |

All existing plotting logic, Plotly figures, and color schemes are untouched.

### New `table_summary.py` per tutorial

Location: `setup<Tutorial>/postProcessing/table_summary.py`

```
tutorials/NiedererEtAl2012/setupNiedererEtAl2012/postProcessing/table_summary.py
tutorials/singleCell/setupSingleCell/postProcessing/table_summary.py
tutorials/manufacturedSolutions/monodomainPseudoECG/setupManufacturedFDA/postProcessing/table_summary.py
tutorials/restitutionCurves_s1s2Protocol/setupRestitutionCurves_s1s2Protocol/postProcessing/table_summary.py
```

Each script:
1. Reads its tutorial's CSV outputs from `output_dir`
2. Extracts key metrics (tutorial-specific columns)
3. Calls `TableWriter.write()` with a `TableMetadata` envelope
4. Returns the artifact list from `TableWriter`

**Tutorial-specific columns:**

| Tutorial | Key metrics | Source |
|----------|------------|--------|
| NiedererEtAl2012 | `case_id`, `DX_mm`, `DT_ms`, `solver`, `point_{n}_activation_ms` (one col per probe point) | `*points_DT*_DX*.csv` files in `output_dir` |
| singleCell | `case_id`, `ionic_model`, `tissue`, `APD_ms`, `peak_voltage_mV`, `resting_voltage_mV` | `.txt` time-series files in `output_dir`; APD extracted by 90% repolarization threshold crossing on voltage column |
| manufacturedFDA | `case_id`, `dimension`, `solver`, `n_cells`, `L2_error`, `convergence_rate` | L2 errors read from `.dat` files in `output_dir`; convergence_rate computed inside the postprocessing script as log-ratio across consecutive refinement levels |
| restitutionCurves | `case_id`, `ionic_model`, `BCL_ms`, `DI_ms`, `APD_ms` | `.txt` trace files; APD/DI detected by threshold crossings as in existing `postProcessing_restCurves.py` |

### File naming conventions

| Output type | Pattern | Example |
|-------------|---------|---------|
| Interactive plot | `<descriptive_name>.html` | `Niederer_vs_cardiacFoam.html` |
| Table CSV | `<TutorialName>_summary.csv` | `NiedererEtAl2012_summary.csv` |
| Table HTML | `<TutorialName>_summary.html` | `NiedererEtAl2012_summary.html` |
| Manifest | `plots.json` | `plots.json` (fixed name, always) |

---

## Section 5 — Files Merged from `driverFOAMagentPreparation`

| Source | Destination | Action |
|--------|-------------|--------|
| `postprocessing/driver.py` | `openfoam_driver/postprocessing/driver.py` | Replace |
| `postprocessing/plot_builder.py` | `openfoam_driver/postprocessing/plot_builder.py` | New file |
| `postprocessing/style.py` | `openfoam_driver/postprocessing/style.py` | Replace |
| `postprocessing/plotting_common.py` | `openfoam_driver/postprocessing/plotting_common.py` | Replace |
| `postprocessing/POSTPROCESSING_ARCHITECTURE.md` | `docs/superpowers/specs/` | Reference doc |

`table_writer.py` is new — written from scratch per this design.

---

## Future Extension Point

A heart-simulation-level `HeartSimulationTableTemplate` will be added later by the user. It will follow the same `TableWriter` + `TableMetadata` pattern with a fixed column schema defined at that time. The `table_summary.py` scripts in each tutorial will remain unchanged — the heart template is an additional layer, not a replacement.
