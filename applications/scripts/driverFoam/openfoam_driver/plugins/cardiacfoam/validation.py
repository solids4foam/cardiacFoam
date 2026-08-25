from typing import TYPE_CHECKING, Any

from openfoam_driver.plugins.cardiacfoam.solver_coupling import SOLVER_COMPATIBILITY_RULES
from openfoam_driver.core.specs.validation_types import ValidationError

if TYPE_CHECKING:
    from openfoam_driver.core.planning_types import StrictDiagnostic


_CONDUCTION_SOLVER_SUFFIX = ".purkinjeGraphModelCoeffs.conductionSystemSolver"
_COUPLER_SUFFIX = ".electroDomainCoupler"
_NETWORK_REF_SUFFIX = ".conductionNetworkDomain"
import re as _re

# A dynamic-path placeholder segment, e.g. <name> / <region_name> / <patch>.
_PLACEHOLDER_RE = _re.compile(r"<[^>]+>")

_CONDUCTION_NET_PREFIX = "conductionNetworkDomains."
_DOMAIN_COUPLINGS_PREFIX = "domainCouplings."


def _is_template_slot_key(key: str) -> bool:
    """Slot keys carrying an un-substituted dynamic-path placeholder
    (e.g. ``domainCouplings.<name>.conductionNetworkDomain``) are template
    forms that ``_filled_run`` synthesises for required-field coverage but
    do not represent a real run-time coupling."""
    return "<" in key or ">" in key


def _find_conduction_system_solver(context: dict[str, Any]) -> str | None:
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
    return [
        str(val) for key, val in context.items()
        if not _is_template_slot_key(key)
        and key.startswith(_DOMAIN_COUPLINGS_PREFIX)
        and key.endswith(_COUPLER_SUFFIX)
    ]


def _evaluate_solver_coupling(context: dict[str, Any]) -> list[ValidationError]:
    errors: list[ValidationError] = []
    myocardium = context.get("myocardiumSolver")
    if myocardium is None:
        return errors

    purkinje = _find_conduction_system_solver(context)
    declared_couplers = _find_declared_couplers(context)

    purkinje_for_rule_match = purkinje if purkinje is not None else None

    for rule in SOLVER_COMPATIBILITY_RULES:
        if rule["myocardium_solver"] != myocardium:
            continue
        rule_purkinje = rule["purkinje_solver"]
        if rule_purkinje == "*":
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

        break

    return errors


def _declared_conduction_networks(context: dict[str, Any]) -> set[str]:
    """Names of every conductionNetworkDomains.<name> block that has at
    least one sub-key present in context (i.e. is actually configured, not
    just a template slot)."""
    declared_networks: set[str] = set()
    for key in context:
        if _is_template_slot_key(key):
            continue
        if not key.startswith(_CONDUCTION_NET_PREFIX):
            continue
        rest = key[len(_CONDUCTION_NET_PREFIX):]
        if "." not in rest:
            continue
        declared_networks.add(rest.split(".", 1)[0])
    return declared_networks


def _evaluate_block_references(
    context: dict[str, Any],
) -> list[ValidationError]:
    errors: list[ValidationError] = []
    declared_networks = _declared_conduction_networks(context)

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


def _value_matches_scoped(actual: Any, expected: str | tuple) -> bool:
    """Scalar equality / tuple membership for a single instance's own value.

    Deliberately does not reuse ``specs.validation._predicate_matches``:
    that helper treats an un-substituted ``<name>`` placeholder as a
    wildcard matching ANY configured instance, which is correct for its
    once-per-catalog-entry call site but wrong here -- this function is
    called once per declared network specifically so one network's value
    can never satisfy another network's requirement.
    """
    if isinstance(expected, tuple):
        return actual in expected
    return actual == expected


def _dynamic_block_templates() -> dict[str, list[Any]]:
    """Group every dynamic catalogue entry by its block template.

    A template is the slot key truncated through its FIRST placeholder --
    ``domainCouplings.<name>``, ``ecgDomains.<name>``,
    ``ionicHeterogeneity.regions.<region_name>``. Entries sharing a template
    are the leaves that one configured instance of that block may carry.
    """
    from openfoam_driver.dict_entries import get_electro_property_entry_groups
    from openfoam_driver.core.specs.validation import slot_key

    templates: dict[str, list[Any]] = {}
    for group in get_electro_property_entry_groups().values():
        for entry in group:
            if not entry.dynamic_path:
                continue
            key = slot_key(entry.driver_path)
            match = _PLACEHOLDER_RE.search(key)
            if not match:
                continue
            templates.setdefault(key[: match.end()], []).append(entry)
    return templates


def _literal_segments_at(position: int, templates: dict[str, list[Any]]) -> set[str]:
    """Literal (non-placeholder) segments any template uses at ``position``.

    ``ecgDomains.<name>`` and ``ecgDomains.electrodePositions.<electrode>``
    share a prefix, so a naive instance scan would read the literal
    ``electrodePositions`` as an ECG domain called "electrodePositions".
    Excluding segments the catalogue itself uses literally at that depth
    keeps the two apart without hardcoding either name.
    """
    literals: set[str] = set()
    for template in templates:
        segments = template.split(".")
        if len(segments) > position and not _PLACEHOLDER_RE.fullmatch(segments[position]):
            literals.add(segments[position])
    return literals


def _declared_instances(
    context: dict[str, Any], template: str, templates: dict[str, list[Any]]
) -> set[str]:
    """Names of every ``<template>`` block with at least one key in context."""
    prefix_segments = template.split(".")[:-1]
    prefix = ".".join(prefix_segments) + "."
    position = len(prefix_segments)
    reserved = _literal_segments_at(position, templates)

    instances: set[str] = set()
    for key in context:
        if _is_template_slot_key(key) or not key.startswith(prefix):
            continue
        remainder = key[len(prefix):].split(".")
        if len(remainder) < 2:
            # The block's own leaf must sit BELOW the instance name.
            continue
        name = remainder[0]
        if name and name not in reserved:
            instances.add(name)
    return instances


def _instance_applicable(
    entry: Any,
    context: dict[str, Any],
    template_prefix: str,
    instance_prefix: str,
) -> bool:
    """Is ``entry`` applicable for this concrete instance?

    Scoped exactly like ``required_when``: a predicate naming a sibling leaf
    inside the same block is resolved against THIS instance's keys, so one
    instance's configuration can never make another instance's leaf
    applicable. Predicates pointing outside the block are left to the generic
    pass and treated as satisfied here.
    """
    from openfoam_driver.core.specs.validation import slot_key

    if not entry.applicable_when:
        return True
    for pred_template, expected in entry.applicable_when.items():
        pred_slot = slot_key(pred_template)
        if not pred_slot.startswith(template_prefix):
            continue
        actual = context.get(instance_prefix + pred_slot[len(template_prefix):])
        if actual is None or not _value_matches_scoped(actual, expected):
            return False
    return True


def _evaluate_dynamic_required_fields(context: dict[str, Any]) -> list[ValidationError]:
    """Required-field checks for every configured dynamic block.

    ``validate_run``'s generic required-field pass (``specs/validation.py``)
    skips every ``dynamic_path`` entry outright: "concrete required leaves
    are the user's responsibility when those blocks are actually
    configured." This is where that responsibility is discharged.

    Each declared instance gets its own scoped substitution of the catalogue
    template, built from only that instance's own keys, so a sibling's value
    can never satisfy this instance's requirement (and vice versa).

    This was previously hardcoded to ``conductionNetworkDomains.<name>.*``,
    which left the other dynamic blocks unenforced -- notably
    ``domainCouplings.<name>``, whose leaves all carry an empty
    ``typical_value`` and so are not silently filled by the builder the way
    ``ecgDomains.<name>.sampling.*`` are. The blocks are now discovered from
    the catalogue instead of named here, so a new dynamic block is covered
    the moment it is catalogued.
    """
    from openfoam_driver.core.specs.validation import slot_key

    errors: list[ValidationError] = []
    templates = _dynamic_block_templates()

    for template, entries in sorted(templates.items()):
        template_prefix = f"{template}."
        for instance in sorted(_declared_instances(context, template, templates)):
            instance_prefix = ".".join(template.split(".")[:-1] + [instance]) + "."
            for entry in entries:
                suffix = slot_key(entry.driver_path)[len(template_prefix):]
                if not suffix:
                    continue
                concrete_key = instance_prefix + suffix

                required = entry.required
                if entry.required_when:
                    required = False
                    for pred_template, expected in entry.required_when.items():
                        pred_slot = slot_key(pred_template)
                        if not pred_slot.startswith(template_prefix):
                            # Predicate references something outside this
                            # instance's own block (e.g. a top-level
                            # selector); out of scope for this per-instance
                            # pass.
                            continue
                        pred_suffix = pred_slot[len(template_prefix):]
                        actual = context.get(instance_prefix + pred_suffix)
                        if actual is None:
                            continue
                        if _value_matches_scoped(actual, expected):
                            required = True
                            break

                if not required:
                    continue
                if not _instance_applicable(entry, context, template_prefix, instance_prefix):
                    # Required only WHERE APPLICABLE. The three
                    # ecgDomains.<name>.sampling.* leaves are gated on
                    # ecgSolver == eikonalECG; an ECG domain using any other
                    # solver never emits them, so demanding them there is a
                    # false error. Every other dynamic required entry has an
                    # empty applicable_when, so this narrows nothing else.
                    continue
                if concrete_key in context and context[concrete_key] not in (None, ""):
                    continue
                block = template.split(".")[:-1]
                errors.append(ValidationError(
                    phase="physics",
                    field=concrete_key,
                    message=(
                        f"{concrete_key} is required for "
                        f"{'.'.join(block)}.{instance} but has no value."
                    ),
                    level="error",
                ))

    return errors


_HETEROGENEITY_PREFIX = "ionicHeterogeneity."
_APEX_BASE_PREFIX = "ionicHeterogeneity.apexBaseBands."


def _evaluate_heterogeneity(context: dict[str, Any]) -> list[ValidationError]:
    errors: list[ValidationError] = []
    het_keys = [k for k in context if k.startswith(_HETEROGENEITY_PREFIX)]
    if not het_keys:
        return errors

    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG

    transmural_keys = [k for k in het_keys if not k.startswith(_APEX_BASE_PREFIX)]
    ab_keys = [k for k in het_keys if k.startswith(_APEX_BASE_PREFIX)]

    model = context.get("ionicModel")
    entry = IONIC_MODEL_CATALOG.get(model) if model is not None else None

    if transmural_keys and entry is not None and not getattr(entry, "supports_heterogeneity", False):
        capable_models = sorted(
            n for n, e in IONIC_MODEL_CATALOG.items()
            if getattr(e, "supports_heterogeneity", False)
            and not n.endswith("compactBatched")
        )
        errors.append(ValidationError(
            phase="physics",
            field="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity",
            message=(
                f"ionicHeterogeneity is configured but ionicModel "
                f"{model!r} does not support transmural heterogeneity. "
                f"Supported models: {', '.join(capable_models)} "
                f"(and their compactBatched variants where available)."
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
            pass

    mode = context.get("ionicHeterogeneity.mode", "transmuralBands")
    if mode == "namedRegions":
        ranges = []
        for k, v in context.items():
            if k.startswith("ionicHeterogeneity.regions.") and k.endswith(".range"):
                region_name = k.split(".")[2]
                try:
                    if isinstance(v, str):
                        clean_v = v.strip("()[] ")
                        parts = clean_v.split()
                        min_v, max_v = float(parts[0]), float(parts[1])
                    else:
                        min_v, max_v = float(v[0]), float(v[1])
                    ranges.append((min_v, max_v, region_name, k))
                except (ValueError, TypeError, IndexError):
                    pass

        ranges.sort(key=lambda x: x[0])

        for i in range(len(ranges)):
            min_v, max_v, name, k = ranges[i]
            if min_v >= max_v:
                errors.append(ValidationError(
                    phase="physics",
                    field=f"$ELECTRO_MODEL_COEFFS.{k}",
                    message=f"Region '{name}' range [{min_v}, {max_v}] must be strictly increasing.",
                    level="error",
                ))
            if min_v < 0.0 or max_v > 1.0:
                errors.append(ValidationError(
                    phase="physics",
                    field=f"$ELECTRO_MODEL_COEFFS.{k}",
                    message=f"Region '{name}' range [{min_v}, {max_v}] must be within [0, 1].",
                    level="error",
                ))
            if i > 0:
                prev_min, prev_max, prev_name, prev_k = ranges[i - 1]
                if min_v < prev_max:
                    errors.append(ValidationError(
                        phase="physics",
                        field=f"$ELECTRO_MODEL_COEFFS.{k}",
                        message=f"Region '{name}' range [{min_v}, {max_v}] overlaps with region '{prev_name}' [{prev_min}, {prev_max}].",
                        level="error",
                    ))

    if ab_keys:
        if entry is not None and not getattr(entry, "supports_apex_base_heterogeneity", False):
            capable_models = sorted(
                n for n, e in IONIC_MODEL_CATALOG.items()
                if getattr(e, "supports_apex_base_heterogeneity", False)
                and not n.endswith("compactBatched")
            )
            errors.append(ValidationError(
                phase="physics",
                field="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.apexBaseBands",
                message=(
                    f"ionicHeterogeneity.apexBaseBands is configured but ionicModel "
                    f"{model!r} does not support apex-to-base heterogeneity. "
                    f"Supported models: {', '.join(capable_models)} "
                    f"(and their compactBatched variants where available)."
                ),
                level="error",
            ))

        beta = context.get("ionicHeterogeneity.apexBaseBands.beta")
        if beta is not None:
            try:
                if float(beta) <= 0:
                    errors.append(ValidationError(
                        phase="physics",
                        field="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.apexBaseBands.beta",
                        message=f"apexBaseBands.beta ({beta}) must be > 0.",
                        level="error",
                    ))
            except (TypeError, ValueError):
                pass

        scaling_min = context.get("ionicHeterogeneity.apexBaseBands.scalingMin")
        scaling_max = context.get("ionicHeterogeneity.apexBaseBands.scalingMax")
        if scaling_min is not None:
            try:
                if float(scaling_min) <= 0:
                    errors.append(ValidationError(
                        phase="physics",
                        field="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.apexBaseBands.scalingMin",
                        message=f"apexBaseBands.scalingMin ({scaling_min}) must be > 0.",
                        level="error",
                    ))
            except (TypeError, ValueError):
                pass
        if scaling_min is not None and scaling_max is not None:
            try:
                if float(scaling_min) > float(scaling_max):
                    errors.append(ValidationError(
                        phase="physics",
                        field="$ELECTRO_MODEL_COEFFS.ionicHeterogeneity.apexBaseBands.scalingMin",
                        message=(
                            f"apexBaseBands.scalingMin ({scaling_min}) must be <= "
                            f"scalingMax ({scaling_max})."
                        ),
                        level="error",
                    ))
            except (TypeError, ValueError):
                pass

    return errors


def _evaluate_tissue_compatibility(context: dict[str, Any]) -> list[ValidationError]:
    errors: list[ValidationError] = []
    model = context.get("ionicModel")
    tissue = context.get("tissue")
    if model is None or tissue is None:
        return errors

    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import IONIC_MODEL_CATALOG

    entry = IONIC_MODEL_CATALOG.get(model)
    if entry is None or not entry.compatible_tissues:
        return errors
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


_RPVJ_COUPLER = "reactionDiffusionPvjCoupler"


def _graph_has_terminal_resistances(graph_path: Any) -> bool:
    """Read-only, structural check of a materialized Purkinje graph file for
    a non-empty top-level ``pvjResistances`` list.

    Mirrors ``conductionGraph::readFromDict`` (src/electroModels/
    electroDomains/conductionSystemDomain/conductionGraph.H): that C++ side
    reads ``pvjResistances`` as an optional top-level scalar list and, if
    present, ``conductionSystemDomain::terminalResistances()`` (conductionSystemDomain.H)
    returns it whenever non-empty. Uses foamlib rather than the line-based
    scanner in core/runtime/mutators.py because this is a pure existence/
    shape check on a file this module did not write and never needs
    verbatim values from -- foamlib parses in-process without evaluating
    ``#calc``/``#codeStream`` (see core/runtime/foam_backend.py's header
    comment), so it is safe to point at an arbitrary materialized file.
    """
    from foamlib import FoamFile

    try:
        resistances = FoamFile(graph_path).get("pvjResistances")
    except (OSError, ValueError):
        return False
    if resistances is None:
        return False
    try:
        return len(resistances) > 0
    except TypeError:
        return bool(resistances)


def _evaluate_pvj_resistance_requirement(
    case_root: Any, electro_path: Any,
) -> tuple["StrictDiagnostic", ...]:
    """Graph-aware requiredness check for ``reactionDiffusionPvjCoupler``'s
    ``rPvj``.

    ``reactionDiffusionPvjCoupler.C`` (src/electroModels/electroCouplers/
    pvjCoupler/reactionDiffusion/reactionDiffusionPvjCoupler.C:120-134)
    only reads ``dict.get<scalar>("rPvj")`` -- a hard ``FatalError`` if
    absent -- when the graph's own ``terminalResistances()`` is null (i.e.
    the graph file provides no ``pvjResistances``). The generic catalog has
    no ``required_when`` predicate that can express "required unless a
    FILE says otherwise", so this is a plugin-specific semantic check
    rather than a ``DictEntry`` field, consumed by ``validate_configuration``
    (the strict pre-flight check gating ``driverFoam run --strict``, which has
    filesystem access to the materialized case) rather than
    ``validate_run_semantics`` (which only ever sees an abstract run
    document, never a real graph file on disk).

    Three-way outcome per ``reactionDiffusionPvjCoupler``-coupled network:

    - The graph file is not yet materialized on disk: DEFER (no
      diagnostic). cardiacCore may still generate it in a later step; a
      later ``--strict`` re-check (e.g. the one immediately before launch)
      will see the materialized file and correctly resolve this.
    - The graph IS materialized and provides ``pvjResistances``, OR
      ``rPvj`` is set directly: silent (either source is sufficient, same
      as the C++ precedence).
    - The graph IS materialized, provides no ``pvjResistances``, and
      ``rPvj`` is absent: ERROR -- neither source exists, which is exactly
      what the C++ side would hard-FatalError on at runtime.
    """
    from pathlib import Path as _Path

    from openfoam_driver.core.planning_types import diagnostic as _diagnostic
    from openfoam_driver.plugins.cardiacfoam.dict_builder import parse_electro_properties
    from openfoam_driver.core.specs.validation import slot_key

    case_root = _Path(case_root)
    electro_path = _Path(electro_path)
    if not electro_path.exists():
        return ()

    overrides = parse_electro_properties(electro_path)["overrides"]

    couplings: dict[str, dict[str, str]] = {}
    for k, v in overrides.items():
        sk = slot_key(k)
        if not sk.startswith(_DOMAIN_COUPLINGS_PREFIX):
            continue
        rest = sk[len(_DOMAIN_COUPLINGS_PREFIX):]
        if rest.endswith(_COUPLER_SUFFIX):
            name = rest[: -len(_COUPLER_SUFFIX)]
            couplings.setdefault(name, {})["coupler"] = v
        elif rest.endswith(_NETWORK_REF_SUFFIX):
            name = rest[: -len(_NETWORK_REF_SUFFIX)]
            couplings.setdefault(name, {})["network"] = v

    diagnostics: list["StrictDiagnostic"] = []
    for coupling_name, info in couplings.items():
        if info.get("coupler") != _RPVJ_COUPLER:
            continue
        network = info.get("network")
        if network is None:
            # Dangling/absent reference is _evaluate_block_references's
            # concern, not this function's.
            continue

        # rPvj lives on the COUPLER's own dict block (domainCouplings.<name>),
        # not on the network's purkinjeGraphModelCoeffs -- it's the argument
        # reactionDiffusionPvjCoupler's own constructor reads via
        # dict.get<scalar>("rPvj") on the dictionary it was constructed
        # with, which is the coupler dict, confirmed against the catalog
        # entry ($ELECTRO_MODEL_COEFFS.domainCouplings.<name>.rPvj) and the
        # real purkinjeNiedererEtAl2011/monodomain1D3D tutorial fixtures.
        rpvj_key = f"$ELECTRO_MODEL_COEFFS.{_DOMAIN_COUPLINGS_PREFIX}{coupling_name}.rPvj"
        if rpvj_key in overrides:
            continue

        graph_key = (
            f"$ELECTRO_MODEL_COEFFS.{_CONDUCTION_NET_PREFIX}{network}"
            f".purkinjeGraphModelCoeffs.graphFile"
        )
        graph_name = overrides.get(graph_key)
        if graph_name is None:
            # No graph reference at all -- required-field checks own
            # flagging a missing graphFile; not this function's concern.
            continue

        graph_path = case_root / "constant" / graph_name
        if not graph_path.exists():
            continue  # DEFER: not yet materialized.

        if _graph_has_terminal_resistances(graph_path):
            continue

        diagnostics.append(_diagnostic(
            "error",
            "missing_rpvj",
            (
                f"conductionNetworkDomains.{network} is coupled via "
                f"{_RPVJ_COUPLER} (domainCouplings.{coupling_name}) but "
                f"neither rPvj nor a graph-provided pvjResistances list is "
                f"available -- the materialized graph file {graph_name!r} "
                f"has no pvjResistances, and rPvj is not set. "
                f"reactionDiffusionPvjCoupler.C will FatalError on "
                f"dict.get<scalar>(\"rPvj\") at runtime."
            ),
            source=str(electro_path),
            field=rpvj_key,
        ))

    return tuple(diagnostics)
