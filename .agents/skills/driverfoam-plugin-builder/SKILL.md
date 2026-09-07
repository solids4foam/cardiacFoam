---
name: driverfoam-plugin-builder
description: >
  Use this skill when a user wants to integrate a NEW OpenFOAM solver into
  driverFOAM. Trigger when a user says "add a new solver", "create a plugin",
  "port myFoam to driverFOAM", "how do I support my solver in driverFOAM",
  "make driverFOAM work with X", or "extend driverFOAM".
---

# driverFOAM Plugin Builder Skill

As an AI agent, your goal is to help users build a working driverFOAM solver
plugin from scratch. driverFOAM is solver-agnostic at its core: any OpenFOAM
solver can be orchestrated by implementing the `SolverPlugin` contract and
registering an entry-point.

---

## What is a driverFOAM Plugin?

A **plugin** is a Python class that adapts a specific OpenFOAM solver
(e.g. `shallowWaterFoam`, `simpleFoam`, `rhoCentralFoam`) for orchestration
by driverFOAM. It has three parts:

1. **A Python class** implementing the `SolverPlugin` (v1) + `SolverPluginV2`
   (v2) contracts defined in `core/plugin_interface.py`.
2. **A `plugin.yaml` manifest** declaring the case files, C++ source roots,
   and identity (id, api version).
3. **An entry-point registration** in `pyproject.toml` so driverFOAM can
   discover the plugin by name at runtime.

The plugin draws a **clean boundary** between the generic execution engine
(DAG runner, RunDocument, strict planner) and all solver-specific knowledge
(which binaries to run, which dictionaries to write, which artifacts to expect).

---

## The 6 Mandatory Artifacts

Before you call `validate_plugin()`, you must have:

| # | Artifact | Location |
|---|---|---|
| 1 | Plugin Python class file | `my_package/my_solver_plugin.py` |
| 2 | `plugin.yaml` manifest | Same directory as the class file |
| 3 | Entry-point in `pyproject.toml` | `[project.entry-points."driverfoam.plugins"]` |
| 4 | `get_dict_entries()` catalog | Inside the class |
| 5 | At least one tutorial in `get_tutorial_catalog()` | (or return the empty struct) |
| 6 | Pass `validate_plugin()` | Confirmed by `load_plugin_context('mysolver')` |

---

## Step-by-Step Build Workflow

### Step 1 — Copy the scaffold

Copy `openfoam_driver/core/generic_plugin.py` to your package as
`my_solver_plugin.py`. Fill in the four identity properties:

```python
@property
def plugin_name(self) -> str:
    return "shallowWaterFoam"          # Human display name

@property
def plugin_id(self) -> str:
    return "org.myproject.shallowwater" # Reverse-DNS; must match plugin.yaml

@property
def plugin_version(self) -> str:
    return "0.1.0"                     # Plugin semantics version

@property
def plugin_api_version(self) -> str:
    return "2"                         # Use "2"; "1" is legacy
```

> **Rule**: `plugin_id` must use only lowercase letters, digits, dots, and
> hyphens, and must not start or end with punctuation.

### Step 2 — Author `plugin.yaml`

Create a `plugin.yaml` alongside your class file. Use the annotated template
in `openfoam_driver/core/generic-plugin.yaml` for field reference. A minimal
example for a solver that reads `constant/transportProperties`:

```yaml
schema_version: 1
plugin:
  id: org.myproject.shallowwater  # Must match plugin_id
  api_version: "2"

case_profile:
  dictionaries:
    - path: system/controlDict
      kind: openfoam_dictionary
      role: openfoam.control_dict
      required: always
    - path: system/fvSchemes
      kind: openfoam_dictionary
      role: openfoam.discretisation
      required: always
    - path: system/fvSolution
      kind: openfoam_dictionary
      role: openfoam.solver_settings
      required: always
    - path: constant/transportProperties
      kind: openfoam_dictionary
      role: plugin.configuration
      required: always
    - path: Allrun
      kind: case_script
      role: openfoam.entrypoint
      required: conditional

# cxx_mapping: omit if you have no C++ sources to audit
```

Load it in `get_profile()`:

```python
from functools import lru_cache
from pathlib import Path
from openfoam_driver.core.plugin_profile import load_plugin_profile

@staticmethod
@lru_cache(maxsize=1)
def get_profile():
    return load_plugin_profile(Path(__file__).with_name("plugin.yaml"))
```

> **Rule**: `profile.plugin_id` and `profile.api_version` must match your
> class's `plugin_id` and `plugin_api_version` or `driver_context()` raises.

### Step 3 — Register the entry-point

In your `pyproject.toml` add:

```toml
[project.entry-points."driverfoam.plugins"]
shallowwater = "my_package.my_solver_plugin:ShallowWaterPlugin"
```

The name (`shallowwater`) is what users pass to `--plugin`. It must be unique
across all installed packages.

Also declare your `plugin.yaml` in `package-data`:

```toml
[tool.setuptools.package-data]
"my_package" = ["plugin.yaml"]
```

Install in editable mode so the entry-point is live:

```bash
pip install -e .
```

Verify discovery:

```bash
python -c "from importlib.metadata import entry_points; print(list(entry_points(group='driverfoam.plugins')))"
```

### Step 4 — Implement the `DictEntry` catalog

Create one `DictEntry` for each key your solver reads from its dictionaries.
Import from `openfoam_driver.core.contracts.dictionary`:

```python
from openfoam_driver.core.contracts.dictionary import DictEntry
from openfoam_driver.core.contracts.dictionary_catalog import DictionaryCatalog

_TRANSPORT_ENTRIES = (
    DictEntry(
        driver_path="transportProperties.nu",
        description="Kinematic viscosity [m2/s]",
        value_kind="openfoam_literal",
        required=True,
        unit="m2 s-1",
        typical_value="1e-6",
    ),
    DictEntry(
        driver_path="transportProperties.model",
        description="Transport model name",
        value_kind="openfoam_literal",
        enum_values=("Newtonian", "CrossPowerLaw"),
        required=True,
    ),
)

def get_dict_entries(self):
    return _TRANSPORT_ENTRIES

def get_dictionary_catalog(self):
    return DictionaryCatalog({"transportProperties": _TRANSPORT_ENTRIES})

def get_dict_groups(self):
    return {"transport": _TRANSPORT_ENTRIES}

def get_dict_entry_catalog(self):
    return {"transportProperties": _TRANSPORT_ENTRIES}
```

> **Rule**: `driver_path` strings must be globally unique across all entries.
> Duplicates cause `driver_context()` to raise immediately.

### Step 5 — Register tutorials

At minimum return the empty catalog structure so the driver knows where to
look if tutorials are added later:

```python
def get_tutorial_catalog(self):
    return {
        "registered_tutorials": (),
        "spec_factories": {},
        "make_generic_case_spec": None,
    }

def get_tutorial_displays(self):
    return ()
```

To register a real tutorial, create a `make_spec` factory and add it:

```python
from openfoam_driver.core.runtime.models import TutorialSpec, CaseConfig

def _make_dam_break_spec(*, case_root, **_):
    return TutorialSpec(
        name="damBreak",
        case_root=case_root,
        config=CaseConfig(solver="shallowWaterFoam"),
    )

def get_tutorial_catalog(self):
    return {
        "registered_tutorials": ("damBreak",),
        "spec_factories": {"damBreak": _make_dam_break_spec},
        "make_generic_case_spec": None,
    }
```

### Step 6 — Declare solver commands (v2)

```python
def get_solver_commands(self) -> frozenset[str]:
    return frozenset({"shallowWaterFoam"})

def get_auxiliary_commands(self) -> frozenset[str]:
    return frozenset({"blockMesh", "decomposePar", "reconstructPar", "checkMesh"})

def get_solve_step_commands(self) -> frozenset[str]:
    return frozenset({"shallowWaterFoam"})
```

Leave the remaining v2 stubs from the scaffold in place (they return empty
values and are valid for most solvers).

### Step 7 — Validate end-to-end

```bash
# Confirm the plugin loads
python -c "
from openfoam_driver.core.plugin_interface import load_plugin_context
ctx = load_plugin_context('shallowwater')
print('Loaded:', ctx.identity.id, 'api_version:', ctx.identity.api_version)
"

# Run strict plan against a tutorial or case folder
driverFoam --plugin shallowwater plan --strict --entry damBreak
```

Fix all diagnostics reported by `plan --strict` before proceeding to run.

### Step 8 (Optional) — Enable sweeps

Without these two hooks, `driverFoam sweep-run` is **refused by name**:

```python
def route_sweep_case_values(self, *, base, resolved_axis_values, driver_context):
    """Map one sweep axis combination onto your solver's case vocabulary.
    Must be pure — no filesystem writes; materialize_sweep_case does those.
    """
    routed = dict(base)
    routed.update(resolved_axis_values)  # e.g. {"nu": 1e-5}
    return routed

def materialize_sweep_case(self, *, case_dir, routed):
    """Write one sweep case to disk from the routed values."""
    # Example: patch transportProperties.nu
    from foamlib import FoamFile
    props = FoamFile(case_dir / "constant" / "transportProperties")
    props["nu"] = routed["nu"]
```

### Step 9 (Optional) — Add post-processing

Create a `postprocess/` directory alongside your plugin. Each script must
expose `run_postprocessing` matching `PostprocessingProtocol`:

```python
from openfoam_driver.postprocessing import PostprocessingProtocol

def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    **kwargs: object,
) -> list[dict]:
    """Loads depth-field CSVs and plots free-surface evolution."""
    # Implementation here
    return []
```

The docstring is used by the agent to decide which script matches the task.

---

## v1 Required Members Cheat-Sheet

All 14 must be present or `validate_plugin()` raises:

```python
class MySolverPlugin:
    @property
    def plugin_name(self) -> str: ...              # "shallowWaterFoam"
    @property
    def plugin_id(self) -> str: ...                # "org.myproject.shallowwater"
    @property
    def plugin_version(self) -> str: ...           # "0.1.0"
    @property
    def plugin_api_version(self) -> str: ...       # "2"
    def get_profile(self): ...                     # load_plugin_profile(...)
    def get_dict_entries(self): ...                # tuple[DictEntry, ...]
    def get_dictionary_catalog(self): ...          # DictionaryCatalog(...)
    def get_dict_groups(self): ...                 # dict[str, tuple[DictEntry, ...]]
    def get_capabilities(self): ...                # build_capability_manifest(...)
    def get_tutorial_catalog(self): ...            # dict with spec_factories
    def get_tutorial_displays(self): ...           # tuple[TutorialDisplay, ...]
    def validate_configuration(self, spec): ...    # tuple[StrictDiagnostic, ...]
    def validate_run_semantics(self, context): ... # tuple[...]
    def predict_data_artifacts(self, case_root, spec): ... # tuple[DataArtifact, ...]
```

## v2 Additional Required Members

All 13 must be **callable** when `plugin_api_version == "2"`:

```python
    def get_solver_commands(self) -> frozenset[str]: ...
    def get_auxiliary_commands(self) -> frozenset[str]: ...
    def get_utility_manifests(self) -> dict: ...
    def get_utility_roots(self) -> tuple: ...
    def resolve_case_models(self, case_root) -> dict: ...
    def get_samplable_fields(self, resolved) -> dict: ...
    def get_override_schema(self, tutorial_name, make_spec_info) -> dict: ...
    def get_run_document_config_schema(self) -> dict: ...
    def get_dict_entry_catalog(self) -> dict: ...
    def get_solve_step_commands(self) -> frozenset[str]: ...
    def get_telemetry_source_globs(self, command) -> tuple: ...
    def get_extra_provenance_paths(self, case_root) -> tuple: ...
    def get_artifact_value_reader(self, artifact_format) -> None: ...
```

---

## Key Optional Hooks

These are probed with `getattr` — implement only what you need:

| Hook | Fallback when absent | Why implement it |
|---|---|---|
| `route_sweep_case_values(...)` | **Sweeps refused** | Required for `driverFoam sweep-run` |
| `materialize_sweep_case(...)` | **Sweeps refused** | Required for `driverFoam sweep-run` |
| `has_case_marker(case_root)` | `False` | Auto-detect if a folder belongs to your plugin |
| `is_nondimensional_case(spec)` | `False` (SI checks on) | Skip mesh-scale diagnostics for non-dimensional cases |
| `get_mesh_geometry_diagnostics(case_root)` | `()` | Add geometry checks for your non-polyMesh point sets |
| `build_run_document_config(spec)` | `({}, ())` | Build a typed `RunDocument.config` object |
| `get_override_scopes()` | `()` | Enable `--apply` patch overrides |
| `get_regeneration_scopes()` | `()` | Enable `--apply` regenerating overrides |
| `get_report_catalog()` | `()` | Register post-run reports |
| `get_named_catalogs()` | `{}` | Expose catalogs in `driverFoam describe` output |
| `get_required_inputs(...)` | `()` (all unknowns = inputs) | Tighten provenance tracking |
| `get_generated_output_globs(...)` | `()` | Improve artifact vs input classification |

---

## Worked Example — `ShallowWaterPlugin`

A complete minimal plugin that passes `validate_plugin()`:

```python
# my_package/shallow_water_plugin.py
"""driverFOAM plugin for the built-in shallowWaterFoam solver."""

from __future__ import annotations
from functools import lru_cache
from pathlib import Path

from openfoam_driver.core.plugin_profile import load_plugin_profile
from openfoam_driver.core.contracts.dictionary import DictEntry
from openfoam_driver.core.contracts.dictionary_catalog import DictionaryCatalog
from openfoam_driver.core.capability_manifest import build_capability_manifest

_TRANSPORT_ENTRIES = (
    DictEntry(
        driver_path="transportProperties.nu",
        description="Kinematic viscosity",
        value_kind="openfoam_literal",
        required=True,
        unit="m2 s-1",
        typical_value="1e-6",
    ),
)


class ShallowWaterPlugin:
    """driverFOAM plugin for shallowWaterFoam."""

    # --- Identity -----------------------------------------------------------
    @property
    def plugin_name(self) -> str:
        return "shallowWaterFoam"

    @property
    def plugin_id(self) -> str:
        return "org.myproject.shallowwater"

    @property
    def plugin_version(self) -> str:
        return "0.1.0"

    @property
    def plugin_api_version(self) -> str:
        return "2"

    # --- Profile ------------------------------------------------------------
    @staticmethod
    @lru_cache(maxsize=1)
    def get_profile():
        return load_plugin_profile(Path(__file__).with_name("plugin.yaml"))

    # --- Dictionary catalog -------------------------------------------------
    def get_dict_entries(self):
        return _TRANSPORT_ENTRIES

    def get_dictionary_catalog(self):
        return DictionaryCatalog({"transportProperties": _TRANSPORT_ENTRIES})

    def get_dict_groups(self):
        return {"transport": _TRANSPORT_ENTRIES}

    def get_dict_entry_catalog(self):
        return {"transportProperties": _TRANSPORT_ENTRIES}

    # --- Capabilities -------------------------------------------------------
    def get_capabilities(self):
        return build_capability_manifest(
            plugin_commands=self.get_solver_commands() | self.get_auxiliary_commands(),
            utility_manifests=self.get_utility_manifests(),
            samplable_fields=self.get_samplable_fields({}),
        )

    # --- Tutorials ----------------------------------------------------------
    def get_tutorial_catalog(self):
        return {"registered_tutorials": (), "spec_factories": {}, "make_generic_case_spec": None}

    def get_tutorial_displays(self):
        return ()

    # --- Validation ---------------------------------------------------------
    def validate_configuration(self, spec):
        return ()

    def validate_run_semantics(self, context):
        return ()

    def predict_data_artifacts(self, case_root, spec):
        return ()

    # --- v2: Commands -------------------------------------------------------
    def get_solver_commands(self) -> frozenset[str]:
        return frozenset({"shallowWaterFoam"})

    def get_auxiliary_commands(self) -> frozenset[str]:
        return frozenset({"blockMesh", "decomposePar", "reconstructPar"})

    def get_utility_manifests(self) -> dict:
        return {}

    def get_utility_roots(self) -> tuple:
        return ()

    # --- v2: Case introspection ---------------------------------------------
    def resolve_case_models(self, case_root) -> dict:
        del case_root
        return {}

    def get_samplable_fields(self, resolved) -> dict:
        del resolved
        return {}

    # --- v2: Configuration vocabulary ---------------------------------------
    def get_override_schema(self, tutorial_name, make_spec_info) -> dict:
        del tutorial_name, make_spec_info
        return {}

    def get_run_document_config_schema(self) -> dict:
        return {"type": "object", "additionalProperties": True}

    # --- v2: Runtime evidence -----------------------------------------------
    def get_solve_step_commands(self) -> frozenset[str]:
        return frozenset({"shallowWaterFoam"})

    def get_telemetry_source_globs(self, command: str) -> tuple:
        del command
        return ()

    def get_extra_provenance_paths(self, case_root) -> tuple:
        del case_root
        return ()

    def get_artifact_value_reader(self, artifact_format: str):
        del artifact_format
        return None
```

Matching `plugin.yaml`:

```yaml
schema_version: 1
plugin:
  id: org.myproject.shallowwater
  api_version: "2"

case_profile:
  dictionaries:
    - path: system/controlDict
      kind: openfoam_dictionary
      role: openfoam.control_dict
      required: always
    - path: constant/transportProperties
      kind: openfoam_dictionary
      role: plugin.configuration
      required: always
    - path: Allrun
      kind: case_script
      role: openfoam.entrypoint
      required: conditional
```

---

## Anti-Patterns to Avoid

1. **Cardiac key leak** — Do not return `{"anatomy": {}, "physics": {}, ...}`
   from `build_run_document_config()` unless your solver actually uses those
   keys. Use `{}` or your own domain-specific keys.

2. **Duplicate `driver_path`** — Every `DictEntry.driver_path` must be unique
   across all entries returned by `get_dict_entries()`. Duplicates raise at
   `driver_context()` construction time.

3. **Wrong entry-point group name** — The group must be exactly
   `driverfoam.plugins` (note: no capitalisation, no underscores).
   A typo produces a silent discovery failure — use `python -c
   "from importlib.metadata import entry_points; print(list(entry_points(group='driverfoam.plugins')))"`.

4. **Profile id / api_version mismatch** — `plugin.yaml → plugin.id` must
   exactly match `plugin_id`, and `plugin.yaml → plugin.api_version` must
   match `plugin_api_version`. A mismatch raises
   `TypeError: SolverPlugin profile id does not match plugin_id`.

5. **Calling `default_driver_context()`** — This returns a cardiac-shaped
   context (for backward compatibility). In your own plugin always use
   `load_plugin_context('mysolver')` or `driver_context(plugin_instance, source=...)`.

---

## Troubleshooting

| Error | Cause | Fix |
|---|---|---|
| `KeyError: 'mysolver'` | Plugin not discovered | Check entry-point group name; reinstall with `pip install -e .` |
| `TypeError: SolverPlugin is missing required members: ...` | Missing v1 methods | Implement the listed methods |
| `TypeError: SolverPlugin declares plugin_api_version '2' but does not implement the v2 contract; missing: ...` | Missing v2 methods | Implement the listed callables |
| `TypeError: SolverPlugin profile id does not match plugin_id` | YAML `plugin.id` ≠ class `plugin_id` | Make them identical strings |
| `TypeError: SolverPlugin.plugin_api_version '2' is not supported` | Wrong api_version string | Use `"1"` or `"2"` only |
| `TypeError: SolverPlugin dictionary catalog has duplicate paths: ...` | Two DictEntry with same `driver_path` | Rename one |
| Sweeps refused: `plugin ... does not implement route_sweep_case_values` | Optional sweep hook absent | Implement Steps 8 |
