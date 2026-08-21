from pathlib import Path

from foamlib import FoamFile


def _read_top_level(electro_properties_path: Path, key: str) -> str:
    try:
        value = FoamFile(electro_properties_path)[key]
    except (KeyError, ValueError) as exc:
        raise KeyError(
            f"Could not determine {key} from {electro_properties_path}"
        ) from exc
    return str(value)


def detect_myocardium_solver_name(electro_properties_path: Path) -> str:
    return _read_top_level(electro_properties_path, "myocardiumSolver")


def detect_electro_coeffs_scope(electro_properties_path: Path) -> str:
    return f"{detect_myocardium_solver_name(electro_properties_path)}Coeffs"


def _coeffs_scope(electro_properties_path: Path):
    """Return the FoamFile.SubDict for the active <solver>Coeffs block.

    Raises the same KeyError shape the pre-migration scanner did if the
    scope itself is missing or the file cannot be parsed at all.
    """
    scope = detect_electro_coeffs_scope(electro_properties_path)
    try:
        return FoamFile(electro_properties_path)[scope]
    except (KeyError, ValueError) as exc:
        raise KeyError(
            f"Could not determine {scope} from {electro_properties_path}"
        ) from exc


def detect_ionic_model_name(electro_properties_path: Path) -> str:
    """Return the ionicModel value from the active <solver>Coeffs block."""
    coeffs = _coeffs_scope(electro_properties_path)
    scope = detect_electro_coeffs_scope(electro_properties_path)
    try:
        return str(coeffs["ionicModel"])
    except KeyError as exc:
        raise KeyError(
            f"Could not determine ionicModel from {electro_properties_path} "
            f"(scope {scope!r})"
        ) from exc


def detect_ionic_export_list(
    electro_properties_path: Path,
) -> tuple[str, ...] | None:
    """Return the names declared in <solver>Coeffs.outputVariables.ionic.export ( ... ).

    Scoped to the active solver's Coeffs block -- the pre-migration regex
    searched the whole file, so an export list belonging to an inactive
    solver's Coeffs block could leak in. This is a deliberate behaviour
    correction alongside the parser migration, not an accidental change.
    """
    try:
        coeffs = _coeffs_scope(electro_properties_path)
        export = coeffs["outputVariables"]["ionic"]["export"]
    except (KeyError, ValueError):
        return None
    if isinstance(export, str):
        return (export,)
    return tuple(str(t) for t in export)


def electro_properties_has_block(
    electro_properties_path: Path,
    block_name: str,
) -> bool:
    """Return True if ``electro_properties_path`` declares a top-level OpenFOAM block."""
    try:
        value = FoamFile(electro_properties_path)[block_name]
    except (KeyError, ValueError):
        return False
    return hasattr(value, "keys")


def detect_verification_model_type(
    electro_properties_path: Path,
) -> str | None:
    """Return the value of ``verificationModel.type`` inside the active ``<solver>Coeffs`` block."""
    try:
        coeffs = _coeffs_scope(electro_properties_path)
        return str(coeffs["verificationModel"]["type"])
    except (KeyError, ValueError):
        return None


def detect_active_tension_model_name(
    electro_properties_path: Path,
) -> str | None:
    """Return the ``activeTensionModel`` value from inside ``<solver>Coeffs``."""
    try:
        coeffs = _coeffs_scope(electro_properties_path)
        return str(coeffs["activeTensionModel"])
    except (KeyError, ValueError):
        return None


def detect_active_tension_export_list(
    electro_properties_path: Path,
) -> tuple[str, ...] | None:
    """Return the names declared in <solver>Coeffs.outputVariables.activeTension.export ( ... )."""
    try:
        coeffs = _coeffs_scope(electro_properties_path)
        export = coeffs["outputVariables"]["activeTension"]["export"]
    except (KeyError, ValueError):
        return None
    if isinstance(export, str):
        return (export,)
    return tuple(str(t) for t in export)
