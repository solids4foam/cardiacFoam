from __future__ import annotations

import json
import shutil
from pathlib import Path


CSV_PATTERNS = (
    "*_points_DT*_DX*.csv",
    "*_line_DT*_DX*.csv",
)


def _case_ids_from_manifest(output_dir: Path) -> list[str]:
    manifest_path = output_dir / "run_manifest.json"
    if not manifest_path.exists():
        return []

    payload = json.loads(manifest_path.read_text())
    results = payload.get("results", [])

    case_ids: list[str] = []
    for item in results:
        if not isinstance(item, dict):
            continue
        case_id = item.get("case_id")
        status = item.get("status")
        if isinstance(case_id, str) and case_id and status == "ok":
            case_ids.append(case_id)

    return case_ids


def _clear_old_csv_cache(output_dir: Path) -> None:
    for pattern in CSV_PATTERNS:
        for csv_file in output_dir.glob(pattern):
            csv_file.unlink()


def _copy_case_csvs_to_output(case_cache_dir: Path, output_dir: Path) -> int:
    copied = 0
    for pattern in CSV_PATTERNS:
        for csv_file in case_cache_dir.glob(pattern):
            shutil.copy2(csv_file, output_dir / csv_file.name)
            copied += 1
    return copied


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    cache_root: str | None = None,
    cache_output_subdir: str = "cachedPostProcessing",
    **_: object,
):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if cache_root is not None:
        cache_root_path = Path(cache_root)
    elif setup_root is not None:
        cache_root_path = Path(setup_root) / "cachedCasePostProcessing"
    else:
        raise ValueError("Either cache_root or setup_root must be provided")

    if not cache_root_path.exists():
        print(f"No cached postProcessing folder found at: {cache_root_path}")
        return []

    requested_case_ids = _case_ids_from_manifest(output_path)
    if requested_case_ids:
        case_ids = requested_case_ids
    else:
        case_ids = sorted(path.name for path in cache_root_path.iterdir() if path.is_dir())

    cache_export_dir = output_path / cache_output_subdir
    if cache_export_dir.exists():
        shutil.rmtree(cache_export_dir)
    cache_export_dir.mkdir(parents=True, exist_ok=True)

    _clear_old_csv_cache(output_path)

    copied_cases = 0
    copied_csv = 0
    for case_id in case_ids:
        case_cache_dir = cache_root_path / case_id
        if not case_cache_dir.exists():
            print(f"Skipping missing cache for case: {case_id}")
            continue

        destination = cache_export_dir / case_id
        shutil.copytree(case_cache_dir, destination)
        copied_cases += 1
        copied_csv += _copy_case_csvs_to_output(case_cache_dir, output_path)

    print(
        "Cached postProcessing restore complete: "
        f"{copied_cases} case folder(s), {copied_csv} CSV file(s)."
    )

    return [
        {
            "path": str(cache_export_dir),
            "label": "Per-case cached postProcessing",
            "kind": "data",
            "format": "dir",
        }
    ]
