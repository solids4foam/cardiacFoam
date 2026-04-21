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
        id="entry_catalog",
        path="/entries",
        title="Entry Catalog",
        purpose="List runnable tutorials, workflow templates, workflow cases, and discovered case folders.",
        primary_data=("entry_catalog", "workflow_catalog", "registered_tutorials", "case_directories"),
        actions=("select entry", "filter by kind", "open config"),
    ),
    RouteSpec(
        id="entry_overview",
        path="/entries/:entryId",
        title="Entry Overview",
        purpose="Show entry classification, workflow metadata, resolved paths, and make_spec schema.",
        primary_data=("entry", "workflow", "resolution", "resolved_name", "make_spec", "spec.metadata"),
        actions=("edit config", "preview cases", "run dry preview"),
    ),
    RouteSpec(
        id="entry_config",
        path="/entries/:entryId/config",
        title="Entry Configuration",
        purpose="Edit entry-level make_spec parameters and dict overrides before execution.",
        primary_data=(
            "make_spec.parameters",
            "common_override_keys",
            "dict_entries",
            "ionic_model_catalog",
            "active_tension_catalog",
        ),
        actions=("save config", "reset defaults", "preview cases"),
    ),
    RouteSpec(
        id="entry_cases",
        path="/entries/:entryId/cases",
        title="Planned Cases",
        purpose="Preview the expanded sweep or generic-case plan produced by the current entry configuration.",
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
        id="EntryCatalogViewModel",
        description="Home screen model for tutorials, workflow templates, workflow cases, and generic case folders.",
        fields=(
            "entry_catalog",
            "workflow_catalog",
            "registered_tutorials",
            "special_tutorial_aliases",
            "available_tutorials",
            "case_directories",
        ),
        source="describe_entry(...), route-independent catalog data",
    ),
    ViewModelSpec(
        id="WorkflowFamilyViewModel",
        description="Workflow family metadata linking symbolic templates to runnable reference cases.",
        fields=(
            "workflow_family",
            "template_entry",
            "reference_cases",
            "workflow_templates",
        ),
        source="describe_entry(...).workflow_catalog[]",
    ),
    ViewModelSpec(
        id="EntryOverviewViewModel",
        description="Resolved entry metadata, workflow context, and top-level configuration summary.",
        fields=(
            "requested_entry",
            "entry",
            "entry_kind",
            "is_runnable",
            "workflow",
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
        source="describe_entry(...)",
    ),
    ViewModelSpec(
        id="SpecParameterViewModel",
        description="One entry-level make_spec parameter rendered as a form field.",
        fields=("name", "kind", "required", "annotation", "default", "current_value"),
        source="describe_entry(...).make_spec.parameters",
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
        source="describe_entry(...).dict_entries",
    ),
    ViewModelSpec(
        id="PlannedCaseViewModel",
        description="One expanded case from the current entry configuration.",
        fields=("case_id", "params"),
        source="describe_entry(...).spec.cases.items",
    ),
    ViewModelSpec(
        id="LaunchRunRequest",
        description="Payload persisted by the GUI and convertible to driver --config JSON.",
        fields=(
            "entry",
            "entry_kind",
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
        fields=("entry", "entry_kind", "case_root", "output_dir", "dry_run", "results"),
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
