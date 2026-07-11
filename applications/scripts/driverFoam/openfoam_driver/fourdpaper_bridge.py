from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path


DEFAULT_FIELDS = ("Vm", "activationTime")
DEFAULT_RENDER_DOC = "main.qmd"
MANIFEST_NAME = "4dpaper_manifest.json"


@dataclass(frozen=True)
class FourDPaperAsset:
    path: str
    kind: str
    source: str


@dataclass(frozen=True)
class FourDPaperBridgeResult:
    status: str
    case_root: str
    paper_root: str
    data_path: str
    qmd_fragment: str
    manifest_path: str
    source_foam: str
    fields: list[str]
    purkinje_series: str | None = None
    ecg_paths: list[str] = field(default_factory=list)
    static_plot_paths: list[str] = field(default_factory=list)
    rendered_assets: list[FourDPaperAsset] = field(default_factory=list)
    render_command: list[str] = field(default_factory=list)
    render_returncode: int | None = None
    render_stdout: str | None = None
    render_stderr: str | None = None
    warnings: list[str] = field(default_factory=list)

    def to_json(self) -> dict:
        return asdict(self)


def _safe_slug(value: str) -> str:
    out = []
    for ch in value:
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        elif ch in ("/", "\\", ".", " "):
            out.append("-")
    slug = "".join(out).strip("-")
    return slug or "case"


def _copy_case_data(source_case: Path, dest_case: Path) -> None:
    if dest_case.exists() or dest_case.is_symlink():
        if dest_case.is_symlink() or dest_case.is_file():
            dest_case.unlink()
        else:
            shutil.rmtree(dest_case)
    ignore = shutil.ignore_patterns(
        "workflow_logs",
        "__pycache__",
        ".pytest_cache",
        ".cache",
    )
    shutil.copytree(source_case, dest_case, ignore=ignore)


def _stage_case(case_root: Path, paper_root: Path, *, copy: bool) -> tuple[Path, list[str]]:
    data_dir = paper_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    dest_case = data_dir / _safe_slug(case_root.name)
    warnings: list[str] = []

    if dest_case.exists() or dest_case.is_symlink():
        if dest_case.is_symlink() or dest_case.is_file():
            dest_case.unlink()
        else:
            shutil.rmtree(dest_case)

    if copy:
        _copy_case_data(case_root, dest_case)
        return dest_case, warnings

    try:
        dest_case.symlink_to(case_root.resolve(), target_is_directory=True)
    except OSError as exc:
        warnings.append(f"Symlink failed ({exc}); copied case data instead.")
        _copy_case_data(case_root, dest_case)
    return dest_case, warnings


def _find_foam_file(case_root: Path) -> Path:
    foam_files = sorted(case_root.glob("*.foam"))
    if foam_files:
        return foam_files[0]
    foam_path = case_root / f"{case_root.name}.foam"
    foam_path.touch()
    return foam_path


def _find_purkinje_series(case_root: Path) -> Path | None:
    candidates = (
        case_root / "postProcessing" / "purkinjeNetworkVTK" / "purkinjeNetwork.vtk.series",
        case_root / "postProcessing" / "purkinjeNetworkVTK" / "combined_purkinje.vtk.series",
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    hits = sorted((case_root / "postProcessing").glob("**/*.vtk.series"))
    return hits[0] if hits else None


def _find_ecg_paths(case_root: Path) -> list[Path]:
    names = ("pseudoECG.dat", "torsoECG.dat", "eikonalECG.dat")
    paths: list[Path] = []
    for name in names:
        paths.extend(sorted((case_root / "postProcessing").glob(f"**/{name}")))
    return paths


def _find_static_plots(case_root: Path) -> list[Path]:
    pp = case_root / "postProcessing"
    if not pp.is_dir():
        return []
    return sorted(p for p in pp.rglob("*.png") if p.is_file())


def _rel_to(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.resolve().relative_to(root.resolve()).as_posix()


def _shortcode_lines(
    *,
    case_slug: str,
    foam_src: str,
    fields: list[str],
    purkinje_src: str | None,
    ecg_paths: list[str],
    static_plot_paths: list[str],
) -> str:
    primary = fields[0] if fields else DEFAULT_FIELDS[0]
    secondary = fields[1] if len(fields) > 1 else primary
    field_list = ",".join(dict.fromkeys(fields or list(DEFAULT_FIELDS)))
    fig = _safe_slug(case_slug)
    lines = [
        f"## {case_slug} Results",
        "",
        "{{< 4d-image "
        f'id="{fig}-vm" src="{foam_src}" field="{primary}" '
        f'fields="{field_list}" time="last" cmap="coolwarm" '
        f'caption="{case_slug}: {primary} at the selected timestep." >}}}}',
        "",
        "{{< 4d-timeseries "
        f'id="{fig}-timeseries" src="{foam_src}" field="{primary}" '
        f'fields="{field_list}" steps="4" height="420px" '
        f'caption="{case_slug}: {primary} evolution across representative timesteps." >}}}}',
        "",
        "{{< 4d-panel "
        f'id="{fig}-panel" layout="2x1" height="520px" camera="sync" '
        f'src1="{foam_src}" id1="{fig}-{primary}" field1="{primary}" '
        f'src2="{foam_src}" id2="{fig}-{secondary}" field2="{secondary}" '
        f'caption="{case_slug}: {primary} and {secondary} with synchronized cameras." >}}}}',
        "",
    ]
    if purkinje_src:
        lines.extend([
            "{{< 4d-multi-image "
            f'id="{fig}-heart-purkinje" '
            f'src1="{foam_src}" field1="{primary}" cmap1="coolwarm" '
            'decimate1="auto" opacity1="0.6" colorbar1="true" '
            f'src2="{purkinje_src}" field2="Vm_V" cmap2="plasma" '
            'decimate2="none" colorbar2="false" time="last" stride="8" '
            f'caption="{case_slug}: myocardium and Purkinje network." >}}}}',
            "",
        ])
    if ecg_paths or static_plot_paths:
        lines.extend(["### Static Results", ""])
    for plot_path in static_plot_paths:
        lines.append(f"![{Path(plot_path).stem}]({plot_path})")
        lines.append("")
    for ecg_path in ecg_paths:
        lines.append(f"- ECG data: `{ecg_path}`")
    if ecg_paths:
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _ensure_include(main_qmd: Path, fragment_rel: str) -> None:
    include_line = f"{{{{< include {fragment_rel} >}}}}"
    if main_qmd.exists():
        text = main_qmd.read_text()
        if include_line in text:
            return
        main_qmd.write_text(text.rstrip() + "\n\n" + include_line + "\n")
        return
    main_qmd.write_text(
        "---\n"
        'title: "cardiacFoam Results"\n'
        "format: html\n"
        "---\n\n"
        f"{include_line}\n"
    )


def _collect_rendered_assets(paper_root: Path, case_slug: str,
                             dashboard_case_dir: Path) -> list[FourDPaperAsset]:
    copied: list[FourDPaperAsset] = []
    out_dir = dashboard_case_dir / "postProcessing" / "4dpaper"
    patterns = ("*.png", "*.html")
    sources = [paper_root / "state" / "figures", paper_root / "_output"]
    for source_root in sources:
        if not source_root.is_dir():
            continue
        for pattern in patterns:
            for src in sorted(source_root.rglob(pattern)):
                if not src.is_file():
                    continue
                name = src.name
                if case_slug not in src.as_posix() and case_slug not in name:
                    continue
                out_dir.mkdir(parents=True, exist_ok=True)
                dest = out_dir / name
                shutil.copy2(src, dest)
                copied.append(FourDPaperAsset(
                    path=dest.relative_to(dashboard_case_dir).as_posix(),
                    kind=src.suffix.lstrip(".") or "asset",
                    source=str(src),
                ))
    return copied


def _run_render(paper_root: Path, *, compose_dir: Path | None,
                render_doc: str) -> tuple[list[str], int, str, str]:
    cwd = compose_dir or paper_root
    command = [
        "docker", "compose", "exec", "-T", "4dpapers",
        "quarto", "render", render_doc, "--to", "html",
    ]
    proc = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    return command, proc.returncode, proc.stdout, proc.stderr


def export_4dpaper(
    *,
    case_root: Path,
    paper_root: Path,
    fields: list[str] | None = None,
    purkinje_series: Path | None = None,
    ecg_paths: list[Path] | None = None,
    static_plot_paths: list[Path] | None = None,
    render: bool = False,
    copy: bool = False,
    compose_dir: Path | None = None,
    render_doc: str = DEFAULT_RENDER_DOC,
    update_main: bool = True,
) -> FourDPaperBridgeResult:
    case_root = Path(case_root).resolve()
    paper_root = Path(paper_root).resolve()
    if not case_root.is_dir():
        raise FileNotFoundError(f"case root not found: {case_root}")
    if not paper_root.is_dir():
        raise FileNotFoundError(f"4Dpapers workspace not found: {paper_root}")

    fields = [f.strip() for f in (fields or list(DEFAULT_FIELDS)) if f.strip()]
    if not fields:
        fields = list(DEFAULT_FIELDS)

    source_foam = _find_foam_file(case_root)
    data_case, warnings = _stage_case(case_root, paper_root, copy=copy)
    staged_foam = data_case / source_foam.name
    foam_src = _rel_to(staged_foam, paper_root)

    purkinje = purkinje_series or _find_purkinje_series(case_root)
    staged_purkinje_src = None
    if purkinje is not None and purkinje.is_file():
        try:
            staged_purkinje_src = _rel_to(data_case / purkinje.relative_to(case_root), paper_root)
        except ValueError:
            staged_purkinje_src = _rel_to(purkinje, paper_root)

    ecg = ecg_paths if ecg_paths is not None else _find_ecg_paths(case_root)
    plots = static_plot_paths if static_plot_paths is not None else _find_static_plots(case_root)
    staged_ecg = []
    for path in ecg:
        try:
            staged_ecg.append(_rel_to(data_case / path.relative_to(case_root), paper_root))
        except ValueError:
            staged_ecg.append(str(path))
    staged_plots = []
    for path in plots:
        try:
            staged_plots.append(_rel_to(data_case / path.relative_to(case_root), paper_root))
        except ValueError:
            staged_plots.append(str(path))

    case_slug = _safe_slug(case_root.name)
    sections = paper_root / "sections"
    sections.mkdir(parents=True, exist_ok=True)
    qmd_path = sections / f"{case_slug}_results.qmd"
    qmd_text = _shortcode_lines(
        case_slug=case_root.name,
        foam_src=foam_src,
        fields=fields,
        purkinje_src=staged_purkinje_src,
        ecg_paths=staged_ecg,
        static_plot_paths=staged_plots,
    )
    qmd_path.write_text(qmd_text)

    if update_main:
        _ensure_include(paper_root / render_doc, qmd_path.relative_to(paper_root).as_posix())

    render_command: list[str] = []
    returncode = None
    stdout = None
    stderr = None
    status = "staged"
    rendered_assets: list[FourDPaperAsset] = []
    if render:
        render_command, returncode, stdout, stderr = _run_render(
            paper_root,
            compose_dir=compose_dir,
            render_doc=render_doc,
        )
        status = "rendered" if returncode == 0 else "render_failed"
        if returncode == 0:
            rendered_assets = _collect_rendered_assets(paper_root, case_slug, case_root)

    manifest_path = case_root / "postProcessing" / "4dpaper" / MANIFEST_NAME
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    result = FourDPaperBridgeResult(
        status=status,
        case_root=str(case_root),
        paper_root=str(paper_root),
        data_path=_rel_to(data_case, paper_root),
        qmd_fragment=str(qmd_path),
        manifest_path=str(manifest_path),
        source_foam=str(source_foam),
        fields=fields,
        purkinje_series=str(purkinje) if purkinje else None,
        ecg_paths=[str(p) for p in ecg],
        static_plot_paths=[str(p) for p in plots],
        rendered_assets=rendered_assets,
        render_command=render_command,
        render_returncode=returncode,
        render_stdout=stdout,
        render_stderr=stderr,
        warnings=warnings,
    )
    manifest_path.write_text(json.dumps(result.to_json(), indent=2))
    return result
