import re
from pathlib import Path

def detect_myocardium_solver_name(electro_properties_path: Path) -> str:
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not stripped.startswith("myocardiumSolver"):
            continue
        tokens = stripped.rstrip(";").split()
        if len(tokens) < 2:
            break
        return tokens[1]
    raise KeyError(f"Could not determine myocardiumSolver from {electro_properties_path}")


def detect_electro_coeffs_scope(electro_properties_path: Path) -> str:
    return f"{detect_myocardium_solver_name(electro_properties_path)}Coeffs"


def detect_ionic_model_name(electro_properties_path: Path) -> str:
    """Return the ionicModel value from the active <solver>Coeffs block."""
    scope = detect_electro_coeffs_scope(electro_properties_path)
    in_scope = False
    depth = 0
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not in_scope:
            if stripped == scope or stripped.startswith(f"{scope} ") or stripped.startswith(f"{scope}{{"):
                in_scope = True
            continue
        if "{" in stripped:
            depth += stripped.count("{")
        if stripped.startswith("ionicModel") and depth == 1:
            tokens = stripped.rstrip(";").split()
            if len(tokens) >= 2:
                return tokens[1]
        if "}" in stripped:
            depth -= stripped.count("}")
            if depth <= 0:
                break
    raise KeyError(
        f"Could not determine ionicModel from {electro_properties_path} "
        f"(scope {scope!r})"
    )


_IONIC_EXPORT_RE = re.compile(r"\bionic\s*\{[^}]*\bexport\s*\(([^)]*)\)", re.DOTALL)


def detect_ionic_export_list(
    electro_properties_path: Path,
) -> tuple[str, ...] | None:
    """Return the names declared in ``outputVariables.ionic.export ( ... )``."""
    text = electro_properties_path.read_text()
    cleaned = "\n".join(line.split("//", 1)[0] for line in text.splitlines())
    match = _IONIC_EXPORT_RE.search(cleaned)
    if match is None:
        return None
    tokens = tuple(t for t in match.group(1).split() if t)
    return tokens if tokens else None


_BLOCK_DECL_RE = re.compile(
    r"^\s*(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*(?:\{|$)",
)


def electro_properties_has_block(
    electro_properties_path: Path,
    block_name: str,
) -> bool:
    """Return True if ``electro_properties_path`` declares a top-level OpenFOAM block."""
    text = electro_properties_path.read_text()
    cleaned_lines = [line.split("//", 1)[0] for line in text.splitlines()]
    for i, line in enumerate(cleaned_lines):
        match = _BLOCK_DECL_RE.match(line)
        if not match or match.group("name") != block_name:
            continue
        if "{" in match.group(0):
            return True
        rest = line[match.end():].lstrip()
        if rest.startswith("{"):
            return True
        for following in cleaned_lines[i + 1:]:
            stripped = following.strip()
            if not stripped:
                continue
            if stripped.startswith("{"):
                return True
            break
    return False


def detect_verification_model_type(
    electro_properties_path: Path,
) -> str | None:
    """Return the value of ``verificationModel.type`` inside the active ``<solver>Coeffs`` block."""
    scope = detect_electro_coeffs_scope(electro_properties_path)
    in_scope = False
    in_verification = False
    depth = 0
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not in_scope:
            if (
                stripped == scope
                or stripped.startswith(f"{scope} ")
                or stripped.startswith(f"{scope}{{")
            ):
                in_scope = True
            continue
        if "{" in stripped:
            depth += stripped.count("{")
        if (
            not in_verification
            and stripped.startswith("verificationModel")
            and depth == 1
        ):
            in_verification = True
        if in_verification and stripped.startswith("type") and depth == 2:
            tokens = stripped.rstrip(";").split()
            if len(tokens) >= 2:
                return tokens[1]
        if "}" in stripped:
            depth -= stripped.count("}")
            if in_verification and depth <= 1:
                in_verification = False
            if depth <= 0:
                break
    return None


_AT_EXPORT_RE = re.compile(
    r"\bactiveTension\s*\{[^}]*\bexport\s*\(([^)]*)\)", re.DOTALL
)


def detect_active_tension_model_name(
    electro_properties_path: Path,
) -> str | None:
    """Return the ``activeTensionModel`` value from inside ``<solver>Coeffs``."""
    scope = detect_electro_coeffs_scope(electro_properties_path)
    in_scope = False
    depth = 0
    for line in electro_properties_path.read_text().splitlines():
        stripped = line.split("//", 1)[0].strip()
        if not in_scope:
            if (
                stripped == scope
                or stripped.startswith(f"{scope} ")
                or stripped.startswith(f"{scope}{{")
            ):
                in_scope = True
            continue
        if depth == 1 and stripped.startswith("activeTensionModel"):
            tokens = stripped.rstrip(";").split()
            if len(tokens) >= 2:
                return tokens[1]
        if "{" in stripped:
            depth += stripped.count("{")
        if "}" in stripped:
            depth -= stripped.count("}")
            if depth <= 0:
                break
    return None


def detect_active_tension_export_list(
    electro_properties_path: Path,
) -> tuple[str, ...] | None:
    """Return the names declared in ``outputVariables.activeTension.export ( ... )``."""
    text = electro_properties_path.read_text()
    cleaned = "\n".join(line.split("//", 1)[0] for line in text.splitlines())
    match = _AT_EXPORT_RE.search(cleaned)
    if match is None:
        return None
    tokens = tuple(t for t in match.group(1).split() if t)
    return tokens if tokens else None
