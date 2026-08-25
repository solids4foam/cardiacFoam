from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from ...core.runtime.mutators import (
    ensure_foam_dict,
    remove_foam_dict,
    remove_foam_entry,
    update_foam_entry,
)
from .detection import detect_electro_coeffs_scope


def _resolve_scope_tokens(
    path: str,
    *,
    electro_properties_path: Path | None = None,
) -> tuple[str, ...]:
    resolved_parts: list[str] = []
    for token in path.split("."):
        if token == "$ELECTRO_MODEL_COEFFS":
            if electro_properties_path is None:
                raise ValueError(
                    "Scope token '$ELECTRO_MODEL_COEFFS' requires electro_properties_path"
                )
            resolved_parts.append(detect_electro_coeffs_scope(electro_properties_path))
            continue
        resolved_parts.append(token)
    return tuple(part for part in resolved_parts if part)


def normalize_entry_overrides(
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
    *,
    electro_properties_path: Path | None = None,
) -> list[dict[str, Any]]:
    if overrides is None:
        return []

    def normalize_from_key_value(key_path: str, value: Any) -> dict[str, Any]:
        parts = _resolve_scope_tokens(
            str(key_path),
            electro_properties_path=electro_properties_path,
        )
        if not parts:
            raise ValueError("Override path cannot be empty")
        if len(parts) == 1:
            return {"key": parts[0], "value": value, "scope": None}
        return {"key": parts[-1], "value": value, "scope": parts[:-1]}

    if isinstance(overrides, Mapping):
        return [normalize_from_key_value(key_path, value) for key_path, value in overrides.items()]

    normalized: list[dict[str, Any]] = []
    for item in overrides:
        if not isinstance(item, Mapping):
            raise TypeError("Entry overrides must be a mapping or sequence of mappings")
        if "key" not in item or "value" not in item:
            raise KeyError("Override items must define 'key' and 'value'")

        key = str(item["key"])
        value = item["value"]
        if "scope" in item:
            raw_scope = item["scope"]
            if raw_scope is None:
                scope = None
            elif isinstance(raw_scope, str):
                scope = _resolve_scope_tokens(
                    raw_scope,
                    electro_properties_path=electro_properties_path,
                )
            else:
                scope = tuple(
                    part
                    for token in raw_scope
                    for part in _resolve_scope_tokens(
                        str(token),
                        electro_properties_path=electro_properties_path,
                    )
                )
        else:
            normalized_item = normalize_from_key_value(key, value)
            normalized.append(normalized_item)
            continue

        normalized.append(
            {
                "key": key,
                "value": value,
                "scope": scope,
            }
        )

    return normalized


def apply_entry_overrides(
    file_path: Path,
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
    *,
    electro_properties_path: Path | None = None,
) -> None:
    for item in normalize_entry_overrides(
        overrides,
        electro_properties_path=electro_properties_path,
    ):
        update_foam_entry(
            file_path,
            item["key"],
            item["value"],
            scope=item["scope"],
        )


def apply_electro_property_overrides(
    electro_properties_path: Path,
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
) -> None:
    apply_entry_overrides(
        electro_properties_path,
        overrides,
        electro_properties_path=electro_properties_path,
    )


def apply_physics_property_overrides(
    physics_properties_path: Path,
    overrides: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None,
) -> None:
    apply_entry_overrides(physics_properties_path, overrides)


def remove_electro_property_dict(
    electro_properties_path: Path,
    dict_name: str,
    *,
    scope: str | Sequence[str] | None = None,
    missing_ok: bool = False,
) -> None:
    resolved_scope = None
    if scope is not None:
        raw_scope = (scope,) if isinstance(scope, str) else tuple(scope)
        resolved_scope = tuple(
            part
            for token in raw_scope
            for part in _resolve_scope_tokens(
                str(token),
                electro_properties_path=electro_properties_path,
            )
        )

    remove_foam_dict(
        electro_properties_path,
        dict_name,
        scope=resolved_scope,
        missing_ok=missing_ok,
    )


def ensure_electro_property_entry(
    electro_properties_path: Path,
    entry_name: str,
    value: Any,
    *,
    scope: str | Sequence[str] | None = None,
) -> None:
    """Set a scalar entry, adding it if the key is absent.

    The general override path deliberately requires a key to exist already, so
    that a typo fails loudly instead of silently growing a new entry. That is
    the right default, but it cannot express a key whose *presence* legitimately
    varies -- the bath-bidomain patch entries, where which of
    ``groundPatches``/``surfaceCurrentPatches`` holds a patch depends on the
    boundary variant the case was last written for.

    Use this only for such entries. Everything else should stay strict.
    """
    resolved_scope = None
    if scope is not None:
        raw_scope = (scope,) if isinstance(scope, str) else tuple(scope)
        resolved_scope = tuple(
            part
            for token in raw_scope
            for part in _resolve_scope_tokens(
                str(token),
                electro_properties_path=electro_properties_path,
            )
        )

    update_foam_entry(
        electro_properties_path,
        entry_name,
        value,
        scope=resolved_scope,
        add_if_missing=True,
    )


def remove_electro_property_entry(
    electro_properties_path: Path,
    entry_name: str,
    *,
    scope: str | Sequence[str] | None = None,
    missing_ok: bool = False,
) -> None:
    """Scalar counterpart of :func:`remove_electro_property_dict`.

    Same ``$TOKEN`` scope resolution; delegates to
    :func:`~openfoam_driver.core.runtime.mutators.remove_foam_entry` so a
    ``name value;`` entry can be removed without the block remover rejecting
    it for having no opening brace.
    """
    resolved_scope = None
    if scope is not None:
        raw_scope = (scope,) if isinstance(scope, str) else tuple(scope)
        resolved_scope = tuple(
            part
            for token in raw_scope
            for part in _resolve_scope_tokens(
                str(token),
                electro_properties_path=electro_properties_path,
            )
        )

    remove_foam_entry(
        electro_properties_path,
        entry_name,
        scope=resolved_scope,
        missing_ok=missing_ok,
    )


def ensure_electro_property_dict(
    electro_properties_path: Path,
    dict_name: str,
    block_text: str,
    *,
    scope: str | Sequence[str] | None = None,
) -> bool:
    resolved_scope = None
    if scope is not None:
        raw_scope = (scope,) if isinstance(scope, str) else tuple(scope)
        resolved_scope = tuple(
            part
            for token in raw_scope
            for part in _resolve_scope_tokens(
                str(token),
                electro_properties_path=electro_properties_path,
            )
        )

    return ensure_foam_dict(
        electro_properties_path,
        dict_name,
        block_text,
        scope=resolved_scope,
    )


def _resolve_electro_model_coeffs_entry(
    driver_path: str, case_root: Path,
) -> tuple[list[str] | None, str]:
    from .detection import detect_myocardium_solver_name
    from .dict_builder import _entry_scope_and_key

    electro_path = case_root / "constant" / "electroProperties"
    coeffs_scope = f"{detect_myocardium_solver_name(electro_path)}Coeffs"
    return _entry_scope_and_key(driver_path, coeffs_scope)


def electro_model_coeffs_scope() -> "OverrideScope":
    """The cardiac plugin's one `step --strict --apply` override scope:
    $ELECTRO_MODEL_COEFFS -> constant/electroProperties, addressed against
    the "electroProperties" dictionary-catalog group, with the active
    solver's <solver>Coeffs block resolved per-case at apply time."""
    from openfoam_driver.core.specs.apply_overrides import OverrideScope

    return OverrideScope(
        token="ELECTRO_MODEL_COEFFS",
        file_relpath="constant/electroProperties",
        catalog_group="electroProperties",
        resolve_entry=_resolve_electro_model_coeffs_entry,
    )


def electro_properties_regeneration_scope() -> "RegenerationScope":
    """The cardiac plugin's one `step --strict --apply` regeneration scope:
    the bare ``myocardiumSolver`` selector -> constant/electroProperties,
    rebuilt (not key-patched) via
    :func:`openfoam_driver.plugins.cardiacfoam.dict_builder.regenerate_electro_properties`
    because switching it renames the active ``<solver>Coeffs`` sub-block and
    changes which sibling keys the catalog allows -- something a single
    key/value/scope patch cannot express. Only ``myocardiumSolver`` is
    wired in: the other three ``_SELECTOR_KEYS`` (``ionicModel``,
    ``tissue``, ``conductivitySource``) are ``$ELECTRO_MODEL_COEFFS.``-scoped
    leaves that change a value in place without renaming anything, so they
    stay on the ordinary key-patch route above."""
    from openfoam_driver.core.specs.apply_overrides import RegenerationScope
    from .dict_builder import regenerate_electro_properties

    return RegenerationScope(
        selector_keys=frozenset({"myocardiumSolver"}),
        file_relpath="constant/electroProperties",
        catalog_group="electroProperties",
        regenerate=regenerate_electro_properties,
    )
