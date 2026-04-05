from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Final


@dataclass(frozen=True)
class RouteSpec:
    id: str
    path: str
    title: str
    purpose: str
    primary_data: tuple[str, ...]
    actions: tuple[str, ...] = ()


@dataclass(frozen=True)
class ViewModelSpec:
    id: str
    description: str
    fields: tuple[str, ...]
    source: str


ROUTES: Final[tuple[RouteSpec, ...]] = (
    RouteSpec(
        id="tutorial_catalog",
        path="/tutorials",
        title="Tutorial Catalog",
        purpose="List registered tutorials, discovered case folders, and entry points into configuration.",
        primary_data=("registered_tutorials", "available_tutorials", "case_directories"),
        actions=("select tutorial", "open config", "open generic case"),
    ),
    RouteSpec(
        id="tutorial_overview",
        path="/tutorials/:tutorialId",
        title="Tutorial Overview",
        purpose="Show tutorial resolution, resolved paths, make_spec schema, and summary metadata.",
        primary_data=("requested_tutorial", "resolution", "resolved_name", "make_spec", "spec.metadata"),
        actions=("edit config", "preview cases", "run dry preview"),
    ),
    RouteSpec(
        id="tutorial_config",
        path="/tutorials/:tutorialId/config",
        title="Tutorial Configuration",
        purpose="Edit tutorial-level make_spec parameters and dict overrides before execution.",
        primary_data=("make_spec.parameters", "common_override_keys", "dict_entries"),
        actions=("save config", "reset defaults", "preview cases"),
    ),
    RouteSpec(
        id="tutorial_cases",
        path="/tutorials/:tutorialId/cases",
        title="Planned Cases",
        purpose="Preview the expanded sweep or generic-case plan produced by the current configuration.",
        primary_data=("spec.cases.count", "spec.cases.items"),
        actions=("filter cases", "inspect params", "launch run"),
    ),
    RouteSpec(
        id="runs_index",
        path="/runs",
        title="Runs",
        purpose="Show prior and in-progress driver runs discovered from run manifests.",
        primary_data=("run_manifest.json",),
        actions=("open run", "refresh"),
    ),
    RouteSpec(
        id="run_detail",
        path="/runs/:runId",
        title="Run Detail",
        purpose="Inspect a single run manifest, case statuses, output paths, and errors.",
        primary_data=("run_manifest.json", "plots.json"),
        actions=("open outputs", "inspect case error", "rerun with config"),
    ),
)


VIEW_MODELS: Final[tuple[ViewModelSpec, ...]] = (
    ViewModelSpec(
        id="TutorialCatalogViewModel",
        description="Home screen model for all runnable tutorials and discovered case folders.",
        fields=(
            "registered_tutorials",
            "special_tutorial_aliases",
            "available_tutorials",
            "case_directories",
        ),
        source="describe_tutorial(...), route-independent catalog data",
    ),
    ViewModelSpec(
        id="TutorialOverviewViewModel",
        description="Resolved tutorial metadata and top-level configuration summary.",
        fields=(
            "requested_tutorial",
            "resolution",
            "resolved_name",
            "factory_overrides",
            "make_spec",
            "spec.name",
            "spec.case_root",
            "spec.setup_root",
            "spec.output_dir",
            "spec.metadata",
        ),
        source="describe_tutorial(...)",
    ),
    ViewModelSpec(
        id="SpecParameterViewModel",
        description="One tutorial-level make_spec parameter rendered as a form field.",
        fields=("name", "kind", "required", "annotation", "default", "current_value"),
        source="describe_tutorial(...).make_spec.parameters",
    ),
    ViewModelSpec(
        id="DictEntryViewModel",
        description="One editable dict override entry with GUI typing hints.",
        fields=(
            "driver_path",
            "description",
            "notes",
            "value_kind",
            "ui_control",
            "enum_values",
            "examples",
            "dynamic_path",
            "source_refs",
        ),
        source="describe_tutorial(...).dict_entries",
    ),
    ViewModelSpec(
        id="PlannedCaseViewModel",
        description="One expanded case from the current tutorial configuration.",
        fields=("case_id", "params"),
        source="describe_tutorial(...).spec.cases.items",
    ),
    ViewModelSpec(
        id="LaunchRunRequest",
        description="Payload persisted by the GUI and convertible to driver --config JSON.",
        fields=(
            "tutorial",
            "tutorials_root",
            "make_spec_overrides",
            "electro_property_overrides",
            "physics_property_overrides",
        ),
        source="GUI-generated from form state",
    ),
    ViewModelSpec(
        id="RunSummaryViewModel",
        description="High-level run status extracted from run_manifest.json.",
        fields=("tutorial", "case_root", "output_dir", "dry_run", "results"),
        source="run_manifest.json",
    ),
    ViewModelSpec(
        id="RunCaseResultViewModel",
        description="Per-case execution state inside a run manifest.",
        fields=("case_id", "status", "duration_s", "params", "error"),
        source="run_manifest.json.results[]",
    ),
)


def describe_gui_schema() -> dict[str, object]:
    return {
        "routes": [asdict(route) for route in ROUTES],
        "view_models": [asdict(view_model) for view_model in VIEW_MODELS],
    }
