"""Run-document validator.

The validator reports three kinds of issue:

1. **Required-field omissions.** Each ``DictEntry`` flagged ``required=True``
   must have a value in the Run document slice owned by its *primary* phase.
2. **Enum violations.** Entries with ``value_kind="enum"`` whose value is not
   one of the declared ``enum_values``.
3. **Structured constraints (P5c/P5b).** Each ``DictEntry`` may declare
   ``applicable_when``, ``forbidden_when``, ``required_when``, and
   ``mutually_exclusive_with``. The validator evaluates these against a
   flattened view of the run config; rule families and semantics are
   documented in plan §5.1. The legacy hardcoded ``eikonalSolver``/
   ``ionicModel`` cross-field check was removed in P5b — the
   ``ionicModel`` entry now carries ``forbidden_when={"myocardiumSolver":
   "eikonalSolver"}`` which is evaluated programmatically here.

The *primary* phase is the first phase in workflow order
(``anatomy → physics → stimulus → solver``) that the entry claims.
Multi-phase entries are validated there; the other phases do not duplicate
validation errors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

from openfoam_driver.dict_entries import (
    DictEntry,
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
    Phase,
)
from openfoam_driver.solver_coupling import SOLVER_COMPATIBILITY_RULES

_PHASE_ORDER: tuple[Phase, ...] = (
    "anatomy", "physics", "stimulus", "solver",
)


@dataclass(frozen=True)
class ValidationError:
    phase: Phase
    field: str
    message: str
    level: str  # "error" | "warning"


def _all_entries():
    yield from PHYSICS_PROPERTY_ENTRIES
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        yield from group


def primary_phase(entry) -> Phase | None:
    """Return the editing phase for a (possibly multi-phase) entry.

    Walks ``_PHASE_ORDER`` and returns the first phase the entry claims;
    every other declared phase is a read-only mirror.
    """
    for ph in _PHASE_ORDER:
        if ph in entry.phases:
            return ph
    return None


_COEFFS_PREFIX = "$ELECTRO_MODEL_COEFFS."


def slot_key(driver_path: str) -> str:
    """Map a driver_path to its slot key inside a phase slice.

    Strips the ``$ELECTRO_MODEL_COEFFS.`` prefix when present; otherwise
    returns the path as-is. Multi-segment unprefixed paths are kept intact
    so that nested-group leaves don't collide with top-level keys of the
    same name (e.g. ``$ELECTRO_MODEL_COEFFS.bathPotentialDomain.phiEReferenceValue`` must not overwrite the
    top-level ``type`` entry inside the physics slice).
    """
    if driver_path.startswith(_COEFFS_PREFIX):
        return driver_path[len(_COEFFS_PREFIX):]
    return driver_path


def _slice_value(run, phase: Phase, driver_path: str):
    """Look up the slot value for a driver_path inside a phase slice."""
    slice_ = run.config.get(phase, {}) or {}
    return slice_.get(slot_key(driver_path))


def _flatten_context(run) -> dict[str, Any]:
    """Build a flat predicate-key → value view across every phase slice.

    Structured-constraint predicates reference slot-keys (the post-
    ``slot_key`` form: ``myocardiumSolver``, ``ionicModel``,
    ``singleCellStimulus.stim_amplitude``, ...). Multiple phases never
    write the same slot-key in practice; if they do, the last wins —
    document the convention rather than silently merging.
    """
    context: dict[str, Any] = {}
    for slice_ in run.config.values():
        if not slice_:
            continue
        for key, val in slice_.items():
            if val not in (None, ""):
                context[key] = val
    return context


def _predicate_matches(
    context: dict[str, Any],
    key: str,
    expected: str | tuple[str, ...],
) -> bool:
    """Return True iff ``context[key]`` matches ``expected``.

    Scalar ``expected`` → equality. Tuple ``expected`` → membership.
    A missing key is treated as not-matching (the predicate's
    precondition is absent).
    """
    if key not in context:
        return False
    actual = context[key]
    if isinstance(expected, tuple):
        return actual in expected
    return actual == expected


def _entry_is_applicable(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Evaluate ``applicable_when`` — every predicate must match (AND)."""
    if not entry.applicable_when:
        return True
    return all(
        _predicate_matches(context, key, expected)
        for key, expected in entry.applicable_when.items()
    )


def is_required_in_context(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Whether an entry's ``required`` semantics fire under this context.

    Many entries declare BOTH ``required=True`` and ``required_when={...}``
    — the author's intent is "required, but only when the predicate
    matches". This helper reads the two fields together:

    - ``required_when`` non-empty → required iff any predicate matches.
    - ``required_when`` empty     → ``entry.required`` is taken at face value.
    """
    if entry.required_when:
        return any(
            _predicate_matches(context, key, expected)
            for key, expected in entry.required_when.items()
        )
    return entry.required


def _entry_value_present(entry: DictEntry, context: dict[str, Any]) -> bool:
    """Is the entry's own slot set in the flattened context?"""
    key = slot_key(entry.driver_path)
    return key in context and context[key] not in (None, "")


def _format_predicate(predicate: dict[str, Any]) -> str:
    parts: list[str] = []
    for key, expected in predicate.items():
        if isinstance(expected, tuple):
            parts.append(f"{key} ∈ {{{', '.join(expected)}}}")
        else:
            parts.append(f"{key}={expected}")
    return " and ".join(parts)


def _all_entries_list():
    return list(_all_entries())


def validate_run(
    run,
    *,
    entries: Iterable[DictEntry] | None = None,
) -> list[ValidationError]:
    """Validate ``run`` against the dict-entry catalog.

    ``entries`` overrides the live catalog for testability and for callers
    that want to validate against a curated subset (e.g., dict_builder).
    When omitted, the full live catalog is used.
    """
    entry_list: list[DictEntry] = (
        list(entries) if entries is not None else _all_entries_list()
    )
    context = _flatten_context(run)
    errors: list[ValidationError] = []

    # 1) Required-field checks. Skip entries whose applicable_when fails —
    #    requiredness is conditional on applicability. When an entry has
    #    both ``required=True`` and a non-empty ``required_when``, the
    #    latter narrows the former: requiredness fires only when at least
    #    one ``required_when`` predicate matches.
    for e in entry_list:
        if not _entry_is_applicable(e, context):
            continue
        if not is_required_in_context(e, context):
            continue
        if e.dynamic_path:
            # Dynamic-path entries describe templates; concrete required
            # leaves are the user's responsibility when those blocks are
            # actually configured.
            continue
        ph = primary_phase(e)
        if ph is None:
            continue
        val = _slice_value(run, ph, e.driver_path)
        if val in (None, ""):
            errors.append(ValidationError(
                phase=ph,
                field=e.driver_path,
                message=f"{e.driver_path} is required.",
                level="error",
            ))

    # 2) Enum checks. Skip inapplicable entries for the same reason.
    for e in entry_list:
        if e.value_kind != "enum" or not e.enum_values:
            continue
        if not _entry_is_applicable(e, context):
            continue
        ph = primary_phase(e)
        if ph is None:
            continue
        val = _slice_value(run, ph, e.driver_path)
        if val is None or val == "":
            continue
        if val not in e.enum_values:
            errors.append(ValidationError(
                phase=ph,
                field=e.driver_path,
                message=f"{val!r} is not one of {list(e.enum_values)}.",
                level="error",
            ))

    # 3) Structured constraints (P5c/P5b).
    # (Legacy hardcoded eikonalSolver+ionicModel check removed in P5b — the
    # ionicModel entry now carries forbidden_when={"myocardiumSolver": "eikonalSolver"}
    # which the section below evaluates programmatically.)
    errors.extend(_evaluate_structured(entry_list, context))

    # 4) Solver-coupling consistency (P5e). Closes the
    # conductionSystemSolver / electroDomainCoupler gap that the four-family
    # structured constraints could not express (cross-domain pairing).
    errors.extend(_evaluate_solver_coupling(context))

    # 5) Block-reference integrity (P5e). Closes the
    # domainCouplings.<name>.conductionNetworkDomain gap (referential
    # integrity to a sibling block).
    errors.extend(_evaluate_block_references(context))

    # 6) Tissue heterogeneity (Phase 2). Cross-field rules the four predicate
    # families cannot express: model-capability gate + interface ordering.
    errors.extend(_evaluate_heterogeneity(context))

    # 7) Tissue / ionic-model compatibility (warning-level). The tissue
    # selector should be one of the chosen model's compatible_tissues.
    errors.extend(_evaluate_tissue_compatibility(context))

    return errors


def _evaluate_structured(
    entries: list[DictEntry],
    context: dict[str, Any],
) -> list[ValidationError]:
    """Evaluate the four structured-constraint families per entry."""
    errors: list[ValidationError] = []
    paths_set = {slot_key(e.driver_path) for e in entries
                 if _entry_value_present(e, context)}

    for e in entries:
        # Skip entries whose applicable_when precondition fails entirely.
        if not _entry_is_applicable(e, context):
            continue
        ph = primary_phase(e) or "physics"

        # forbidden_when: fires when ANY predicate matches AND the entry's
        # own slot has a value. Each matching predicate emits its own
        # ValidationError so the reason text stays specific.
        if _entry_value_present(e, context):
            for key, expected in e.forbidden_when.items():
                if _predicate_matches(context, key, expected):
                    errors.append(ValidationError(
                        phase=ph,
                        field=e.driver_path,
                        message=(
                            f"{e.driver_path} is forbidden when "
                            f"{_format_predicate({key: expected})}."
                        ),
                        level="error",
                    ))

        # required_when: handled by section 1 of validate_run via
        # is_required_in_context. Section 3 does NOT re-emit a violation
        # to avoid double-firing on the same entry. The structured
        # required_when field is still consumed — its predicates feed
        # is_required_in_context which gates section 1's required check.

        # mutually_exclusive_with: fires when BOTH this entry's slot is set
        # AND any of the listed sibling slots is also set. To avoid
        # double-reporting symmetric relations we only flag the side that
        # *declares* the relation.
        if _entry_value_present(e, context):
            for sibling_path in e.mutually_exclusive_with:
                sibling_slot = slot_key(sibling_path)
                if sibling_slot in paths_set:
                    errors.append(ValidationError(
                        phase=ph,
                        field=e.driver_path,
                        message=(
                            f"{e.driver_path} is mutually exclusive with "
                            f"{sibling_path}."
                        ),
                        level="error",
                    ))

    return errors


# -------- P5e: Cross-block and pairing validators --------
#
# The three constraints listed in plan §5 as prose-only
# (conductionSystemSolver, electroDomainCoupler, conductionNetworkDomain)
# share two non-DictEntry features:
#   1. their predicates need a wildcard scan over slot_keys with dynamic
#      <name> segments, and
#   2. their rules are already captured elsewhere (SOLVER_COMPATIBILITY_RULES
#      for the first two; the conductionNetworkDomains block declarations
#      for the third).
#
# We honour the convergence principle: both evaluators consult existing
# tables/context rather than introducing a fifth DictEntry family.


_CONDUCTION_SOLVER_SUFFIX = ".purkinjeGraphModelCoeffs.conductionSystemSolver"
_COUPLER_SUFFIX = ".electroDomainCoupler"
_NETWORK_REF_SUFFIX = ".conductionNetworkDomain"
_CONDUCTION_NET_PREFIX = "conductionNetworkDomains."
_DOMAIN_COUPLINGS_PREFIX = "domainCouplings."


def _is_template_slot_key(key: str) -> bool:
    """Slot keys carrying an un-substituted dynamic-path placeholder
    (e.g. ``domainCouplings.<name>.conductionNetworkDomain``) are template
    forms that ``_filled_run`` synthesises for required-field coverage but
    do not represent a real run-time coupling. Both P5e evaluators skip
    them so they do not generate false-positive dangling-reference or
    coupler-mismatch errors on stub-filled fixtures."""
    return "<" in key or ">" in key


def _find_conduction_system_solver(context: dict[str, Any]) -> str | None:
    """Locate the conductionSystemSolver value if any conductionNetworkDomain
    block declares one. Returns the solver name or None if no Purkinje
    coupling is declared."""
    for key, val in context.items():
        if _is_template_slot_key(key):
            continue
        if (
            key.startswith(_CONDUCTION_NET_PREFIX)
            and key.endswith(_CONDUCTION_SOLVER_SUFFIX)
        ):
            return str(val)
    return None


def _find_declared_couplers(context: dict[str, Any]) -> list[str]:
    """Return every electroDomainCoupler value declared under
    domainCouplings.<name>.electroDomainCoupler."""
    return [
        str(val) for key, val in context.items()
        if not _is_template_slot_key(key)
        and key.startswith(_DOMAIN_COUPLINGS_PREFIX)
        and key.endswith(_COUPLER_SUFFIX)
    ]


def _evaluate_solver_coupling(context: dict[str, Any]) -> list[ValidationError]:
    """Enforce SOLVER_COMPATIBILITY_RULES on the (myocardium, purkinje)
    pair declared in context, plus the required coupler.

    Silent when no Purkinje coupling is declared (the common case for
    single-cell / pure-myocardium runs).
    """
    errors: list[ValidationError] = []
    myocardium = context.get("myocardiumSolver")
    if myocardium is None:
        return errors

    purkinje = _find_conduction_system_solver(context)
    declared_couplers = _find_declared_couplers(context)

    # Build a "no Purkinje declared" placeholder so the rules table still
    # works for the bidomain-without-Purkinje invalid case (where the rule
    # uses "*" as the wildcard for the Purkinje side).
    purkinje_for_rule_match = purkinje if purkinje is not None else None

    for rule in SOLVER_COMPATIBILITY_RULES:
        if rule["myocardium_solver"] != myocardium:
            continue
        # Wildcard '*' on the rule side matches any non-None Purkinje value;
        # otherwise the values must match exactly.
        rule_purkinje = rule["purkinje_solver"]
        if rule_purkinje == "*":
            # The wildcard rules only fire when SOMETHING declares Purkinje.
            if purkinje_for_rule_match is None:
                continue
        elif rule_purkinje != purkinje_for_rule_match:
            continue

        if not rule["valid"]:
            errors.append(ValidationError(
                phase="physics",
                field="myocardiumSolver/conductionSystemSolver",
                message=(
                    f"Incompatible solver pair: myocardiumSolver={myocardium} "
                    f"with conductionSystemSolver={purkinje_for_rule_match}. "
                    f"{rule.get('reason', '')}"
                ).strip(),
                level="error",
            ))
            continue

        # Valid pair: check the required coupler is the one declared.
        required = rule.get("required_coupler")
        if required is None:
            continue
        if not declared_couplers:
            errors.append(ValidationError(
                phase="physics",
                field="electroDomainCoupler",
                message=(
                    f"electroDomainCoupler is required for myocardiumSolver="
                    f"{myocardium} + conductionSystemSolver={purkinje}; "
                    f"expected {required}."
                ),
                level="error",
            ))
        else:
            for actual in declared_couplers:
                if actual != required:
                    errors.append(ValidationError(
                        phase="physics",
                        field="electroDomainCoupler",
                        message=(
                            f"electroDomainCoupler={actual!r} is incompatible "
                            f"with myocardiumSolver={myocardium} + "
                            f"conductionSystemSolver={purkinje}; "
                            f"expected {required}."
                        ),
                        level="error",
                    ))

    return errors


def _evaluate_block_references(
    context: dict[str, Any],
) -> list[ValidationError]:
    """Every domainCouplings.<name>.conductionNetworkDomain must point at a
    network name that has at least one declared sub-key under
    conductionNetworkDomains.<that-name>.*.

    Empty context, no domainCouplings, or no conductionNetworkDomain
    references → no errors emitted.
    """
    errors: list[ValidationError] = []
    declared_networks: set[str] = set()
    for key in context:
        if _is_template_slot_key(key):
            continue
        if not key.startswith(_CONDUCTION_NET_PREFIX):
            continue
        rest = key[len(_CONDUCTION_NET_PREFIX):]
        if "." not in rest:
            continue  # malformed; only count fully-qualified declarations
        declared_networks.add(rest.split(".", 1)[0])

    for key, val in context.items():
        if _is_template_slot_key(key):
            continue
        if not (
            key.startswith(_DOMAIN_COUPLINGS_PREFIX)
            and key.endswith(_NETWORK_REF_SUFFIX)
        ):
            continue
        referenced = str(val)
        if referenced not in declared_networks:
            errors.append(ValidationError(
                phase="physics",
                field=key,
                message=(
                    f"conductionNetworkDomain references {referenced!r} but "
                    f"no matching block is declared under "
                    f"conductionNetworkDomains.{referenced}.*"
                ),
                level="error",
            ))

    return errors


# -------- Phase 2: ionic heterogeneity + tissue compatibility --------
#
# These cross-field checks cannot be expressed by the four DictEntry
# predicate families (numeric ordering, catalog cross-reference), so they
# are evaluated here against the flattened context, mirroring the
# solver-coupling and block-reference evaluators above.

_HETEROGENEITY_PREFIX = "ionicHeterogeneity."


def _evaluate_heterogeneity(context: dict[str, Any]) -> list[ValidationError]:
    """Validate an ``ionicHeterogeneity`` block when one is present.

    Fires only when at least one ``ionicHeterogeneity.*`` slot is set, so
    runs without heterogeneity are unaffected. Two rules:

    1. The selected ``ionicModel`` must support transmural heterogeneity.
    2. ``endoMInterface`` must be strictly less than ``mEpiInterface``.
    """
    errors: list[ValidationError] = []
    het_keys = [k for k in context if k.startswith(_HETEROGENEITY_PREFIX)]
    if not het_keys:
        return errors

    from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG

    model = context.get("ionicModel")
    if model is not None:
        entry = IONIC_MODEL_CATALOG.get(model)
        if entry is not None and not getattr(entry, "supports_heterogeneity", False):
            errors.append(ValidationError(
                phase="physics",
                field="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity",
                message=(
                    f"ionicHeterogeneity is configured but ionicModel "
                    f"{model!r} does not support transmural heterogeneity. "
                    f"Supported models: BuenoOrovio, TNNP, TWorld, "
                    f"ToRORd_dynCl (and their compactBatched variants)."
                ),
                level="error",
            ))

    endo = context.get("ionicHeterogeneity.endoMInterface")
    mepi = context.get("ionicHeterogeneity.mEpiInterface")
    if endo is not None and mepi is not None:
        try:
            if float(endo) >= float(mepi):
                errors.append(ValidationError(
                    phase="physics",
                    field="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.endoMInterface",
                    message=(
                        f"endoMInterface ({endo}) must be strictly less than "
                        f"mEpiInterface ({mepi})."
                    ),
                    level="error",
                ))
        except (TypeError, ValueError):
            # Non-numeric values are caught by value-kind handling elsewhere.
            pass

    return errors


def _evaluate_tissue_compatibility(context: dict[str, Any]) -> list[ValidationError]:
    """Reject a ``tissue`` selector not in the model's compatible set.

    Error-level: the C++ ``ionicSelector`` hard-fails when the tissue is not
    among the model's ``supportedTissueTypes()``, so an incompatible pairing
    cannot run. We reject it pre-launch rather than let the solver abort.

    Skipped when the model is unknown (not in the catalogue) or declares no
    compatible tissues, to avoid false positives. Manufactured models never
    carry a ``tissue`` selector (their dict entry's ``applicable_when``
    excludes them), so they are not affected.
    """
    errors: list[ValidationError] = []
    model = context.get("ionicModel")
    tissue = context.get("tissue")
    if model is None or tissue is None:
        return errors

    from openfoam_driver.ionic_model_catalog import IONIC_MODEL_CATALOG

    entry = IONIC_MODEL_CATALOG.get(model)
    if entry is None or not entry.compatible_tissues:
        return errors
    # Manufactured/verification models use the `dimension` selector, not the
    # biological tissue selector; their `tissue` value is inert, so do not
    # police it (the "manufactured" sentinel is not a real tissue name).
    if "manufactured" in entry.compatible_tissues:
        return errors
    if tissue not in entry.compatible_tissues:
        errors.append(ValidationError(
            phase="physics",
            field="$ELECTRO_MODEL_COEFFS.tissue",
            message=(
                f"tissue {tissue!r} is not in the compatible tissues for "
                f"ionicModel {model!r}: {list(entry.compatible_tissues)}."
            ),
            level="error",
        ))

    return errors
