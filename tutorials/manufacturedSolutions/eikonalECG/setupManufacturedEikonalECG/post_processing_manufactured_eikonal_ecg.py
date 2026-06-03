"""Post-processing helpers for manufactured eikonal ECG runs."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import re
import sys


CASE_PATTERN = re.compile(r"(?P<dimension>\dD)_(?P<cells>\d+)_cells")


def _case_metadata(path: Path) -> dict[str, str]:
    match = CASE_PATTERN.search(path.name)
    if match:
        return {
            "Dimension": match.group("dimension"),
            "N": match.group("cells"),
        }

    for parent in path.parents:
        match = CASE_PATTERN.search(parent.name)
        if match:
            return {
                "Dimension": match.group("dimension"),
                "N": match.group("cells"),
            }
    return {"Dimension": "unknown", "N": "unknown"}


def _as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _safe_rate(coarse_error: float, fine_error: float, coarse_n: int, fine_n: int) -> float:
    if (
        coarse_error <= 0.0
        or fine_error <= 0.0
        or coarse_n <= 0
        or fine_n <= coarse_n
        or not math.isfinite(coarse_error)
        or not math.isfinite(fine_error)
    ):
        return math.nan

    return math.log(coarse_error/fine_error)/math.log(float(fine_n)/float(coarse_n))


def _format_float(value: float) -> str:
    if not math.isfinite(value):
        return ""
    return f"{value:.8e}"


def _parse_activation_summary(path: Path) -> dict[str, str]:
    row: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.split()
        if tokens[0] == "activationTime" and len(tokens) >= 4:
            row.update(
                {
                    "activation_L1": tokens[1],
                    "activation_L2": tokens[2],
                    "activation_Linf": tokens[3],
                }
            )
    return row


def _parse_ecg_summary(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    samples = ""
    tau0 = ""
    grad_tau = ""
    planar_fit_linf = ""

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if tokens[0] == "samples" and len(tokens) >= 2:
            samples = tokens[1]
        elif tokens[0] == "tau0" and len(tokens) >= 2:
            tau0 = tokens[1]
        elif tokens[0] == "gradTau" and len(tokens) >= 2:
            grad_tau = " ".join(tokens[1:])
        elif tokens[0] == "planarFitLinf" and len(tokens) >= 2:
            planar_fit_linf = tokens[1]
        elif tokens[0] == "Electrode":
            header = tokens
        elif header and len(tokens) == len(header):
            row = dict(zip(header, tokens))
            row["samples"] = samples
            row["tau0"] = tau0
            row["gradTau"] = grad_tau
            row["planarFitLinf"] = planar_fit_linf
            rows.append(row)

    return rows


def _activation_convergence_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rates: list[dict[str, str]] = []
    grouped: dict[str, list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        grouped.setdefault(row["Dimension"], []).append(row)

    for dimension, group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["N"]))
        for coarse, fine in zip(ordered, ordered[1:]):
            coarse_n = int(coarse["N"])
            fine_n = int(fine["N"])
            rate_row = {
                "Dimension": dimension,
                "N_coarse": str(coarse_n),
                "N_fine": str(fine_n),
            }

            for field in ("activation_L1", "activation_L2", "activation_Linf"):
                rate_row[f"rate_{field}"] = _format_float(
                    _safe_rate
                    (
                        _as_float(coarse.get(field, "")),
                        _as_float(fine.get(field, "")),
                        coarse_n,
                        fine_n,
                    )
                )

            rates.append(rate_row)

    return rates


def _ecg_convergence_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rates: list[dict[str, str]] = []
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        key = (row["Dimension"], row["Electrode"])
        grouped.setdefault(key, []).append(row)

    for (dimension, electrode), group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["N"]))
        for coarse, fine in zip(ordered, ordered[1:]):
            coarse_n = int(coarse["N"])
            fine_n = int(fine["N"])
            rates.append(
                {
                    "Dimension": dimension,
                    "Electrode": electrode,
                    "N_coarse": str(coarse_n),
                    "N_fine": str(fine_n),
                    "rate_L1_err_ref": _format_float(
                        _safe_rate
                        (
                            _as_float(coarse.get("L1_err_ref", "")),
                            _as_float(fine.get("L1_err_ref", "")),
                            coarse_n,
                            fine_n,
                        )
                    ),
                    "rate_L2_err_ref": _format_float(
                        _safe_rate
                        (
                            _as_float(coarse.get("L2_err_ref", "")),
                            _as_float(fine.get("L2_err_ref", "")),
                            coarse_n,
                            fine_n,
                        )
                    ),
                    "rate_Linf_err_ref": _format_float(
                        _safe_rate
                        (
                            _as_float(coarse.get("Linf_err_ref", "")),
                            _as_float(fine.get("Linf_err_ref", "")),
                            coarse_n,
                            fine_n,
                        )
                    ),
                }
            )

    return rates


_DELTA_Q_PAT = re.compile(r"Linf_delta_q(\d+)_ref")


def _q_checks_from_rows(rows: list[dict[str, str]]) -> list[int]:
    """Extract unique check quadrature orders from row column names."""
    orders: set[int] = set()
    for row in rows:
        for key in row:
            m = _DELTA_Q_PAT.fullmatch(key)
            if m:
                orders.add(int(m.group(1)))
    return sorted(orders)


def _ecg_quadrature_aggregate_rows(
    rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[int]]:
    """Aggregate |ECG(q_check) - ECG(q_ref)| delta columns over electrodes."""
    q_checks = _q_checks_from_rows(rows)
    if not q_checks:
        return [], []

    aggregated: list[dict[str, str]] = []
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        grouped.setdefault((row["Dimension"], row["N"]), []).append(row)

    for (dimension, n_value), group in sorted(
        grouped.items(),
        key=lambda item: (item[0][0], int(item[0][1])),
    ):
        summary = {"Dimension": dimension, "N": n_value}
        for q in q_checks:
            field = f"Linf_delta_q{q}_ref"
            values = [_as_float(row.get(field, "")) for row in group]
            values = [v for v in values if math.isfinite(v)]
            if values:
                summary[f"max_{field}"] = _format_float(max(values))
                summary[f"mean_{field}"] = _format_float(sum(values) / len(values))
            else:
                summary[f"max_{field}"] = ""
                summary[f"mean_{field}"] = ""
        aggregated.append(summary)

    return aggregated, q_checks


def _ecg_aggregate_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    aggregated: list[dict[str, str]] = []
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}

    for row in rows:
        if row.get("N", "unknown") == "unknown":
            continue
        grouped.setdefault((row["Dimension"], row["N"]), []).append(row)

    for (dimension, n_value), group in sorted(
        grouped.items(),
        key=lambda item: (item[0][0], int(item[0][1])),
    ):
        summary = {"Dimension": dimension, "N": n_value}

        for field in ("L1_err_ref", "L2_err_ref", "Linf_err_ref", "planarFitLinf"):
            values = [_as_float(row.get(field, "")) for row in group]
            values = [value for value in values if math.isfinite(value)]
            if values:
                summary[f"max_{field}"] = _format_float(max(values))
                summary[f"mean_{field}"] = _format_float(sum(values)/len(values))
            else:
                summary[f"max_{field}"] = ""
                summary[f"mean_{field}"] = ""

        aggregated.append(summary)

    return aggregated


def _write_report(
    path: Path,
    activation_rows: list[dict[str, str]],
    activation_rate_rows: list[dict[str, str]],
    ecg_aggregate_rows: list[dict[str, str]],
    ecg_rate_rows: list[dict[str, str]],
) -> None:
    lines: list[str] = [
        "# Manufactured eikonal ECG post-processing",
        "",
        "## Activation-time error",
        "",
        "| Dimension | N | L1 | L2 | Linf |",
        "|---|---:|---:|---:|---:|",
    ]

    for row in sorted(
        activation_rows,
        key=lambda item: (item.get("Dimension", ""), int(item.get("N", "0"))),
    ):
        if row.get("N", "unknown") == "unknown":
            continue
        lines.append(
            "| {Dimension} | {N} | {activation_L1} | {activation_L2} | {activation_Linf} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## Activation-time convergence rates",
            "",
            "| Dimension | N coarse | N fine | rate L1 | rate L2 | rate Linf |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )

    for row in activation_rate_rows:
        lines.append(
            "| {Dimension} | {N_coarse} | {N_fine} | {rate_activation_L1} | {rate_activation_L2} | {rate_activation_Linf} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## ECG total error summary",
            "",
            "ECG errors are aggregated over electrodes. They compare numerical ECG computed from numerical psi against the manufactured reference ECG from exact psi.",
            "",
            "| Dimension | N | max Linf ECG | mean Linf ECG | max planar-fit Linf |",
            "|---|---:|---:|---:|---:|",
        ]
    )

    for row in ecg_aggregate_rows:
        lines.append(
            "| {Dimension} | {N} | {max_Linf_err_ref} | {mean_Linf_err_ref} | {max_planarFitLinf} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## ECG convergence rates",
            "",
            "| Dimension | Electrode | N coarse | N fine | rate L1 | rate L2 | rate Linf |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
    )

    for row in ecg_rate_rows:
        lines.append(
            "| {Dimension} | {Electrode} | {N_coarse} | {N_fine} | {rate_L1_err_ref} | {rate_L2_err_ref} | {rate_Linf_err_ref} |".format(**row)
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def postprocess(output_dir: str | Path) -> None:
    root = Path(output_dir)
    activation_rows: list[dict[str, str]] = []
    ecg_rows: list[dict[str, str]] = []

    for summary in root.rglob("*manufacturedEikonalActivationTime.dat"):
        row = _case_metadata(summary)
        row.update(_parse_activation_summary(summary))
        row["source"] = str(summary)
        activation_rows.append(row)

    for summary in root.rglob("*manufacturedEikonalECGSummary.dat"):
        metadata = _case_metadata(summary)
        for row in _parse_ecg_summary(summary):
            row.update(metadata)
            row["source"] = str(summary)
            ecg_rows.append(row)

    activation_rows.sort(
        key=lambda row: (row.get("Dimension", ""), int(row.get("N", "0")))
    )
    ecg_rows.sort(
        key=lambda row:
        (
            row.get("Dimension", ""),
            int(row.get("N", "0")),
            row.get("Electrode", ""),
        )
    )

    activation_rate_rows = _activation_convergence_rows(activation_rows)
    ecg_rate_rows = _ecg_convergence_rows(ecg_rows)
    ecg_aggregate_rows = _ecg_aggregate_rows(ecg_rows)
    ecg_quadrature_rows, q_checks = _ecg_quadrature_aggregate_rows(ecg_rows)

    _write_csv(
        root / "manufacturedEikonalActivationSummary.csv",
        [
            "Dimension",
            "N",
            "activation_L1",
            "activation_L2",
            "activation_Linf",
            "source",
        ],
        activation_rows,
    )

    _write_csv(
        root / "manufacturedEikonalActivationConvergenceRates.csv",
        [
            "Dimension",
            "N_coarse",
            "N_fine",
            "rate_activation_L1",
            "rate_activation_L2",
            "rate_activation_Linf",
        ],
        activation_rate_rows,
    )

    _write_csv(
        root / "manufacturedEikonalECGSummary.csv",
        [
            "Dimension",
            "N",
            "Electrode",
            "samples",
            "tau0",
            "gradTau",
            "planarFitLinf",
            "L1_err_ref",
            "L2_err_ref",
            "Linf_err_ref",
            "source",
        ],
        ecg_rows,
    )

    _write_csv(
        root / "manufacturedEikonalECGAggregateSummary.csv",
        [
            "Dimension",
            "N",
            "max_L1_err_ref",
            "mean_L1_err_ref",
            "max_L2_err_ref",
            "mean_L2_err_ref",
            "max_Linf_err_ref",
            "mean_Linf_err_ref",
            "max_planarFitLinf",
            "mean_planarFitLinf",
        ],
        ecg_aggregate_rows,
    )

    _write_csv(
        root / "manufacturedEikonalECGConvergenceRates.csv",
        [
            "Dimension",
            "Electrode",
            "N_coarse",
            "N_fine",
            "rate_L1_err_ref",
            "rate_L2_err_ref",
            "rate_Linf_err_ref",
        ],
        ecg_rate_rows,
    )

    if ecg_quadrature_rows:
        q_delta_fields = ["Dimension", "N"]
        for q in q_checks:
            q_delta_fields += [f"max_Linf_delta_q{q}_ref", f"mean_Linf_delta_q{q}_ref"]
        _write_csv(
            root / "manufacturedEikonalECGQuadratureSummary.csv",
            q_delta_fields,
            ecg_quadrature_rows,
        )

    _write_report(
        root / "manufacturedEikonalECGReport.md",
        activation_rows,
        activation_rate_rows,
        ecg_aggregate_rows,
        ecg_rate_rows,
    )


def run_postprocessing(
    *,
    output_dir: str,
    setup_root: str | None = None,
    **_: object,
) -> list[dict]:
    """driverFOAM post-processing entry point (called for 'all' runs)."""
    root = Path(output_dir)
    postprocess(root)

    artifacts: list[dict] = [
        {
            "path": str(root / "manufacturedEikonalActivationSummary.csv"),
            "label": "Eikonal activation-time error table",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalActivationConvergenceRates.csv"),
            "label": "Eikonal activation-time convergence rates",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalECGAggregateSummary.csv"),
            "label": "Eikonal ECG aggregate error table",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalECGConvergenceRates.csv"),
            "label": "Eikonal ECG convergence rates",
            "kind": "table",
            "format": "csv",
        },
        {
            "path": str(root / "manufacturedEikonalECGReport.md"),
            "label": "Eikonal ECG verification report",
            "kind": "report",
            "format": "markdown",
        },
        {
            "path": str(root / "manufacturedEikonalECGQuadratureSummary.csv"),
            "label": "Eikonal ECG quadrature delta table",
            "kind": "table",
            "format": "csv",
        },
    ]

    try:
        import importlib.util
        _here = Path(__file__).parent
        spec = importlib.util.spec_from_file_location(
            "plot_convergence", _here / "plot_convergence.py"
        )
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.plot(root)
        mod.plot_quadrature(root)
        for stem, label in (
            ("convergence_plot", "Convergence plot"),
            ("quadrature_plot",  "Quadrature convergence plot"),
        ):
            for suffix in ("pdf", "png"):
                p = root / f"{stem}.{suffix}"
                if p.exists():
                    artifacts.append(
                        {
                            "path": str(p),
                            "label": f"{label} ({suffix})",
                            "kind": "plot",
                            "format": suffix,
                        }
                    )
    except Exception as exc:  # noqa: BLE001
        print(f"Warning: plots skipped ({exc})")

    return artifacts


def main(argv: list[str]) -> int:
    output_dir = Path(argv[1]) if len(argv) > 1 else Path("postProcessing")
    postprocess(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
