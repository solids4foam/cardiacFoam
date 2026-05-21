"""Scratch-construct OpenFOAM dictionary files from agent intent.

`build_electro_properties` synthesises a complete `constant/electroProperties`
text from selectors + overrides. The pipeline reuses the existing dict-entry
catalog (`dict_entries.py`), the structured-constraint validator
(`validation.py`), and the path conventions encoded in `slot_key`.

Convergence (plan §9.1): no new catalog, no new validator. The builder
composes existing primitives. Every output passes through `validate_run`
before being returned; an agent that gets a string back is guaranteed it
is validator-clean.
"""
from __future__ import annotations

from typing import Any

from openfoam_driver.dict_entries import (
    DictEntry,
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
)
from openfoam_driver.core.runtime.run_model import RunDocument
from openfoam_driver.specs.validation import (
    _entry_is_applicable,
    _predicate_matches,
    primary_phase,
    slot_key,
    validate_run,
)


def _all_electro_entries() -> list[DictEntry]:
    out: list[DictEntry] = []
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        out.extend(group)
    return out


# Block-presence virtual keys auto-inferred from override prefixes.
# When an agent declares any override under `bathPotentialDomain.*`,
# `ecgDomains.*`, or `conductionNetworkDomains.*`, the dict_builder sets
# the corresponding `$..._present`/`$..._configured` virtual key on the
# resolved context. Entries gated by `applicable_when={"$..._...": True}`
# then become visible and their typical_value fallbacks fire. Without an
# override under one of these prefixes the block stays off — plain
# bidomain does not get bath leaves leaking in.
_VIRTUAL_PRESENCE_TRIGGERS: tuple[tuple[str, str], ...] = (
    ("bathPotentialDomain.", "$bathPotentialDomain_configured"),
    ("ecgDomains.", "$ecgDomains_present"),
    ("conductionNetworkDomains.", "$conductionNetworkDomains_present"),
)


def _infer_virtual_presence(ctx: dict[str, Any]) -> None:
    """Set block-presence virtual keys in-place on `ctx` whenever any
    real slot_key starts with one of the documented block prefixes.
    Idempotent and safe to call after the override merge."""
    for prefix, virtual_key in _VIRTUAL_PRESENCE_TRIGGERS:
        if virtual_key in ctx:
            continue
        for existing_key in ctx:
            if existing_key.startswith(prefix):
                ctx[virtual_key] = True
                break


def resolve_context(
    selectors: dict[str, str],
    *,
    overrides: dict[str, str] | None = None,
) -> dict[str, Any]:
    """Collapse selectors + overrides into a single `{slot_key: value}` dict.

    Selectors enter at their raw key (`myocardiumSolver`, `ionicModel`, ...);
    overrides go through `slot_key` so the `$ELECTRO_MODEL_COEFFS.` prefix
    is stripped. This is the same context shape `validation._flatten_context`
    produces from a Run document — so the validator can be reused unchanged.

    Also infers virtual presence keys (`$bathPotentialDomain_configured`,
    `$ecgDomains_present`, `$conductionNetworkDomains_present`) so the
    matching `applicable_when` predicates fire only when the agent has
    actually declared overrides under the corresponding block.
    """
    ctx: dict[str, Any] = dict(selectors)
    if overrides:
        for driver_path, value in overrides.items():
            ctx[slot_key(driver_path)] = value
    _infer_virtual_presence(ctx)
    return ctx


def select_applicable_entries(
    context: dict[str, Any],
    *,
    entries: list[DictEntry] | None = None,
) -> list[DictEntry]:
    """Return only entries whose `applicable_when` predicate matches the
    context. Entries with no `applicable_when` are always included. The
    `entries` kwarg lets callers (and tests) curate the input pool;
    default is every electro entry."""
    pool = entries if entries is not None else _all_electro_entries()
    return [e for e in pool if _entry_is_applicable(e, context)]


def _is_required_in_context(
    entry: DictEntry,
    context: dict[str, Any],
) -> bool:
    """Decide whether an entry is required *for this context*.

    Three cases:
    - `required=True` AND `required_when` empty → always required.
    - `required=True` AND `required_when` non-empty → required ONLY when at
      least one `required_when` predicate matches. This reads the two fields
      together as the entry author's intent ("required, but only under
      these conditions"). Without this rule, entries like `phiERefPoint`
      (required=True, required_when={"myocardiumSolver":"bidomainSolver"})
      would fire missing-required errors on singleCell cases.
    - `required=False` AND `required_when` non-empty → required only when a
      predicate matches (the validator's existing semantics).
    """
    if entry.required_when:
        return any(
            _predicate_matches(context, key, expected)
            for key, expected in entry.required_when.items()
        )
    return entry.required


def check_required(
    entries: list[DictEntry],
    populated: dict[str, str],
    *,
    context: dict[str, Any] | None = None,
) -> None:
    """Raise `ValueError` if any required entry in `entries` is missing from
    `populated`. Inapplicable entries are assumed already filtered out by
    `select_applicable_entries`; optional entries are silently ignored.

    `dynamic_path=True` entries are skipped — they describe template paths
    (e.g. ``domainCouplings.<name>.electroDomainCoupler``) rather than
    concrete required leaves. The user's overrides supply concrete paths
    when those blocks are actually configured; the P5e validators catch
    dangling references at run-construction time.

    The optional `context` enables `_is_required_in_context` to honour
    `required_when` predicates; when omitted, `required=True` is treated
    unconditionally for backward compat with simple callers.
    """
    missing: list[str] = []
    ctx = context if context is not None else {}
    for entry in entries:
        if entry.dynamic_path:
            continue
        if not _is_required_in_context(entry, ctx):
            continue
        key = slot_key(entry.driver_path)
        if key not in populated or populated[key] in (None, ""):
            missing.append(entry.driver_path)
    if missing:
        raise ValueError(
            "build_electro_properties: required entries have no value and "
            "no typical_value fallback was applicable:\n  - "
            + "\n  - ".join(missing)
        )


def populate_values(
    entries: list[DictEntry],
    context: dict[str, Any],
    *,
    typical_value_fallback: bool = True,
) -> dict[str, str]:
    """For each entry, resolve the final value to write into the dict.

    Precedence per entry:
      1. Explicit value already in `context` (from selectors or overrides).
      2. `entry.typical_value` if non-empty AND `typical_value_fallback`.
      3. Omit — caller's downstream `check_required` decides whether that's
         a problem for required entries.
    """
    populated: dict[str, str] = {}
    for entry in entries:
        key = slot_key(entry.driver_path)
        if key in context and context[key] not in (None, ""):
            populated[key] = str(context[key])
            continue
        if typical_value_fallback and entry.typical_value:
            populated[key] = entry.typical_value
            continue
        # Omit. Downstream check_required will flag if it was required.
    return populated


def _foamfile_preamble(object_name: str) -> str:
    """Standard FoamFile preamble for a given `object` value.

    The object name distinguishes the two dict files this module writes —
    `electroProperties` vs `physicsProperties`. Future dict targets follow
    the same pattern.
    """
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | cardiacFoam dict_builder synthesis              |\n"
        "|   \\\\  /    A nd           |                                                 |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        "    class       dictionary;\n"
        "    location    \"constant\";\n"
        f"    object      {object_name};\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    )


# Backwards-compat alias — pre-P9b callers and tests reference this name.
_FOAMFILE_PREAMBLE = _foamfile_preamble("electroProperties")


def _populated_to_run(
    populated: dict[str, str],
    entries: list[DictEntry],
) -> RunDocument:
    """Distribute populated values into a Run document keyed by each
    entry's primary phase. Selector keys (which may not correspond to any
    entry, but always do here for the dict-builder entry pool) are placed
    in the physics slice as a sensible default."""
    config: dict[str, dict[str, str]] = {
        "anatomy": {}, "physics": {}, "stimulus": {}, "solver": {},
    }
    placed: set[str] = set()
    for entry in entries:
        key = slot_key(entry.driver_path)
        if key not in populated:
            continue
        ph = primary_phase(entry) or "physics"
        config[ph][key] = populated[key]
        placed.add(key)
    # Any populated keys without a matching entry land in physics. This
    # only triggers for selector keys that don't correspond to DictEntry —
    # uncommon, but safe.
    for key, val in populated.items():
        if key not in placed:
            config["physics"][key] = val
    return RunDocument(id="dict_builder", name="dict_builder",
                       status="draft", config=config)


def build_electro_properties(
    selectors: dict[str, str],
    *,
    overrides: dict[str, str] | None = None,
    typical_value_fallback: bool = True,
) -> str:
    """Synthesise a complete `electroProperties` dict from intent.

    Args:
        selectors: top-level discriminators (myocardiumSolver, ionicModel,
            tissue, ...). Required keys depend on the chosen solver.
        overrides: full driver_path → value mappings for entries whose
            `typical_value` is not appropriate.
        typical_value_fallback: when True (default), applicable entries
            with no override fall back to `DictEntry.typical_value` if
            declared. When False, only explicit overrides count.

    Returns:
        OpenFOAM-format text including the standard `FoamFile` preamble.

    Raises:
        ValueError: required+applicable entry has no value, mutex violation,
            or any structured-constraint violation from `validate_run`.
    """
    context = resolve_context(selectors, overrides=overrides)
    entries = select_applicable_entries(context)
    populated = populate_values(
        entries, context, typical_value_fallback=typical_value_fallback,
    )

    # Safety net: run the full validator scoped to electro entries only.
    # The validator now subsumes required-field checks (its section 1
    # honours `required_when` + `dynamic_path` the same way `check_required`
    # does), so we don't pre-call `check_required` from the public builder
    # entry-point. `check_required` stays exported for callers that want
    # just the required-field subset.
    run = _populated_to_run(populated, entries)
    errors = [e for e in validate_run(run, entries=entries) if e.level == "error"]
    if errors:
        raise ValueError(
            "build_electro_properties: validator rejected synthesised dict:\n  - "
            + "\n  - ".join(e.message for e in errors)
        )

    body = _serialize(populated, entries, selectors["myocardiumSolver"])
    return _FOAMFILE_PREAMBLE + "\n" + body


_COEFFS_PREFIX = "$ELECTRO_MODEL_COEFFS."


def _set_nested(node: dict, path: list[str], value: str) -> None:
    """Insert `value` at `path` inside the nested dict `node`, creating
    intermediate sub-dicts as needed. A leaf already present is
    overwritten — the populated dict has unique slot_keys so this is safe."""
    cursor = node
    for segment in path[:-1]:
        cursor = cursor.setdefault(segment, {})
        # If a prior leaf collided with a sub-block name, replace the leaf
        # with a sub-block — should not happen with the current catalog but
        # is defensive.
        if not isinstance(cursor, dict):
            raise ValueError(
                f"Path collision in serialiser at segment {segment!r}: a leaf "
                f"value exists where a sub-block is needed."
            )
    cursor[path[-1]] = value


def _serialize_block(tree: dict, indent: int) -> str:
    """Emit nested OpenFOAM block syntax. Leaves are `key value;`,
    sub-blocks are `key\\n{\\n  ...\\n}` recursively."""
    lines: list[str] = []
    pad = " " * indent
    for key, value in tree.items():
        if isinstance(value, dict):
            lines.append(f"{pad}{key}")
            lines.append(f"{pad}{{")
            inner = _serialize_block(value, indent + 4)
            if inner:
                lines.append(inner)
            lines.append(f"{pad}}}")
        else:
            lines.append(f"{pad}{key} {value};")
    return "\n".join(lines)


def _serialize(
    populated: dict[str, str],
    entries: list[DictEntry],
    myocardium_solver: str,
) -> str:
    """Group populated values by scope and emit the OpenFOAM dict body.

    Top-level keys (entries whose `driver_path` does not start with the
    `$ELECTRO_MODEL_COEFFS.` prefix) are emitted at the root. Everything
    else nests under the resolved `<solver>Coeffs` block.
    """
    top_level: dict[str, str] = {}
    coeffs: dict = {}
    for entry in entries:
        key = slot_key(entry.driver_path)
        if key not in populated:
            continue
        value = populated[key]
        if entry.driver_path.startswith(_COEFFS_PREFIX):
            segments = entry.driver_path[len(_COEFFS_PREFIX):].split(".")
            _set_nested(coeffs, segments, value)
        else:
            top_level[entry.driver_path] = value

    parts: list[str] = []
    for key, value in top_level.items():
        parts.append(f"{key} {value};")
    if coeffs:
        coeffs_scope = f"{myocardium_solver}Coeffs"
        parts.append("")
        parts.append(coeffs_scope)
        parts.append("{")
        parts.append(_serialize_block(coeffs, indent=4))
        parts.append("}")
    parts.append("")
    return "\n".join(parts)


def build_physics_properties(
    selectors: dict[str, str],
    *,
    overrides: dict[str, str] | None = None,
    typical_value_fallback: bool = True,
) -> str:
    """Synthesise a complete `physicsProperties` dict from intent.

    Mirrors :func:`build_electro_properties` but against
    :data:`PHYSICS_PROPERTY_ENTRIES`. There is no ``<solver>Coeffs``
    wrapper — every physics key lives at the dict root. Today the only
    entry is ``type``; future physics-level selectors slot in unchanged.

    Args:
        selectors: top-level physics keys (e.g. ``{"type": "electroModel"}``).
        overrides: full driver_path → value mappings for future expansion.
        typical_value_fallback: when True, applicable entries with no
            override fall back to ``DictEntry.typical_value`` if declared.

    Returns:
        OpenFOAM-format text with a ``physicsProperties``-typed FoamFile
        preamble.

    Raises:
        ValueError: required entry has no value, or any structured
            constraint violation from `validate_run`.
    """
    context = resolve_context(selectors, overrides=overrides)
    # Scope to physics entries — electro entries don't belong here.
    entries = select_applicable_entries(context, entries=list(PHYSICS_PROPERTY_ENTRIES))
    populated = populate_values(
        entries, context, typical_value_fallback=typical_value_fallback,
    )

    run = _populated_to_run(populated, entries)
    errors = [e for e in validate_run(run, entries=entries) if e.level == "error"]
    if errors:
        raise ValueError(
            "build_physics_properties: validator rejected synthesised dict:\n  - "
            + "\n  - ".join(e.message for e in errors)
        )

    # Physics keys are root-level — no <solver>Coeffs wrapper.
    body_lines: list[str] = []
    for entry in entries:
        key = slot_key(entry.driver_path)
        if key in populated:
            body_lines.append(f"{key} {populated[key]};")
    body = "\n".join(body_lines) + "\n"
    return _foamfile_preamble("physicsProperties") + "\n" + body
