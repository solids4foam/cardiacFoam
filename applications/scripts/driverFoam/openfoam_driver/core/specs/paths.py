from pathlib import Path


def repo_root_default() -> Path:
    """Locate the repository / package root using a three-tier fallback.

    Tier 1 (monorepo): ancestor directory that has both ``tutorials/`` and
        ``src/`` siblings — the full cardiacFoam checkout.
    Tier 2 (standalone-with-tutorials): ancestor directory that has a
        ``tutorials/`` sibling but no ``src/`` — driverFOAM cloned with a
        companion tutorials tree.
    Tier 3 (fully standalone): the ``driverFoam/`` directory that contains
        ``openfoam_driver/`` — i.e. the package root itself.  This is the
        fallback when neither a monorepo nor an external tutorials tree is
        present, e.g. in a temp-folder clone or CI.
    """
    current = Path(__file__).resolve()
    tier2_candidate: Path | None = None
    for parent in current.parents:
        has_tutorials = (parent / "tutorials").exists()
        has_src = (parent / "src").exists()
        # Tier 1: full monorepo layout
        if has_tutorials and has_src:
            return parent
        # Remember first ancestor with tutorials/ only (Tier 2)
        if has_tutorials and tier2_candidate is None:
            tier2_candidate = parent
    if tier2_candidate is not None:
        return tier2_candidate
    # Tier 3: the directory containing openfoam_driver/ (driverFoam/ package root)
    # current = …/openfoam_driver/core/specs/paths.py → parents[2] = openfoam_driver/
    # parents[3] = driverFoam/ (the package root with pyproject.toml).
    # (Verified empirically, not just by counting path segments by eye —
    # this exact line was off by one at the pre-move depth, silently dead
    # code because Tier 1/2 always win inside this monorepo checkout.)
    return current.parents[3]


def tutorials_root_default() -> Path:
    repo_root = repo_root_default()
    tutorials_root = repo_root / "tutorials"
    if tutorials_root.exists():
        return tutorials_root
    return repo_root


def driverfoam_scratch_root() -> Path:
    """Return the repository-local root for disposable driverFOAM data."""
    return repo_root_default() / ".tmp" / "driverfoam"


def default_sweep_output_dir(spec_path: str | Path) -> Path:
    """Return the standard output location for a sweep specification."""
    return driverfoam_scratch_root() / "sweeps" / Path(spec_path).stem


def default_setup_dir_name(case_dir_name: str) -> str:
    normalized_case_dir = case_dir_name.strip()
    if not normalized_case_dir:
        raise ValueError("case_dir_name cannot be empty")
    case_leaf = Path(normalized_case_dir).name
    return f"setup{case_leaf[:1].upper()}{case_leaf[1:]}"


def resolve_spec_paths(
    *,
    tutorials_root: Path | None,
    case_dir_name: str,
    setup_dir_name: str | None = None,
    output_dir_name: str | Path | None = None,
    default_output_dir_name: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    resolved_tutorials_root = (
        Path(tutorials_root) if tutorials_root is not None else tutorials_root_default()
    )
    resolved_case_dir = case_dir_name.strip()
    if not resolved_case_dir:
        raise ValueError("case_dir_name cannot be empty")
    resolved_setup_dir = setup_dir_name or default_setup_dir_name(resolved_case_dir)
    resolved_output_dir = output_dir_name or default_output_dir_name
    if resolved_output_dir is None:
        raise ValueError("output_dir_name or default_output_dir_name must be provided")

    case_root = resolved_tutorials_root / resolved_case_dir
    setup_root = case_root / resolved_setup_dir
    output_dir = case_root / Path(resolved_output_dir)
    return case_root, setup_root, output_dir


def resolve_run_script_path(
    *,
    tutorials_root: Path | None,
    run_script_relpath: Path,
) -> Path:
    if run_script_relpath.is_absolute():
        return run_script_relpath

    candidate_roots: list[Path] = []
    if tutorials_root is not None:
        candidate_roots.append(Path(tutorials_root))
    candidate_roots.append(repo_root_default())
    candidate_roots.append(tutorials_root_default())

    checked_paths: list[Path] = []
    seen: set[Path] = set()
    for root in candidate_roots:
        resolved_root = root.resolve()
        if resolved_root in seen:
            continue
        seen.add(resolved_root)
        candidate = resolved_root / run_script_relpath
        checked_paths.append(candidate)
        if candidate.exists():
            return candidate

    checked_str = ", ".join(str(path) for path in checked_paths)
    raise FileNotFoundError(
        f"Run script not found for '{run_script_relpath}'. Checked: {checked_str}. "
        "Use an absolute path or a path relative to repository root."
    )
