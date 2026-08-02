from typing import Any

from openfoam_driver.plugins.cardiacfoam.solver_coupling import SOLVER_COMPATIBILITY_RULES
from .validation_types import ValidationError


_CONDUCTION_SOLVER_SUFFIX = ".purkinjeGraphModelCoeffs.conductionSystemSolver"
_COUPLER_SUFFIX = ".electroDomainCoupler"
_NETWORK_REF_SUFFIX = ".conductionNetworkDomain"
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


def _evaluate_block_references(
    context: dict[str, Any],
) -> list[ValidationError]:
    errors: list[ValidationError] = []
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
