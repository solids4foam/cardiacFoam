"""Read each conforming case's native output into canonical rows (no rates yet).

Each function returns a list of dicts carrying the non-rate canonical keys
(case, variant, dim, N, h, field, L1, L2, Linf); rates are applied later by
schema.fill_rates. Native files are only read, never rewritten.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path


def from_tet_scheme_study(path, case="tet"):
    rows = []
    with Path(path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            base = dict(case=case, variant=rec["scheme"], dim="3D",
                        N=rec["N"], h=rec["dx"])
            rows.append({**base, "field": "Vm", "L1": "",
                         "L2": rec["mono_L2"], "Linf": rec["mono_Linf"]})
            rows.append({**base, "field": "Phi_e", "L1": "",
                         "L2": rec["ecg_L2"], "Linf": rec["ecg_Linf"]})
    return rows


def from_eikonal_tet_scheme_study(path, case="eikonal_tet"):
    """eikonalTetMMS scheme_study.csv: activationTime + ecg columns per scheme/N.
    Mirrors from_tet_scheme_study; dx -> h, L1 unavailable (blank)."""
    rows = []
    with Path(path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            base = dict(case=case, variant=rec["scheme"], dim="3D",
                        N=rec["N"], h=rec["dx"])
            rows.append({**base, "field": "activationTime", "L1": "",
                         "L2": rec["activationTime_L2"], "Linf": rec["activationTime_Linf"]})
            rows.append({**base, "field": "Phi_e", "L1": "",
                         "L2": rec["ecg_L2"], "Linf": rec["ecg_Linf"]})
    return rows


def from_coupling_summary(path, regime, case="coupling"):
    rows = []
    with Path(path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            common = dict(case=case, variant=regime, N=rec["N"], h=rec["h"],
                          field="Vm")
            rows.append({**common, "dim": "3D",
                         "L1": rec.get("L1_3D_Vm", ""), "L2": rec.get("L2_3D_Vm", ""),
                         "Linf": rec.get("Linf_3D_Vm", "")})
            rows.append({**common, "dim": "1D",
                         "L1": rec.get("L1_1D_Vm", ""), "L2": rec.get("L2_1D_Vm", ""),
                         "Linf": rec.get("Linf_1D_Vm", "")})
    return rows


def _parse_activation_dat(path):
    """Return (L1, L2, Linf) strings from a line 'activationTime L1 L2 Linf'."""
    for line in Path(path).read_text(errors="ignore").splitlines():
        parts = line.split()
        if parts and parts[0].lower().startswith("activationtime") and len(parts) >= 4:
            return parts[1], parts[2], parts[3]
    raise ValueError(f"no activationTime line in {path}")


def from_eikonal_activation(summary_path, extra_2d_dats=None, case="eikonal"):
    rows = []
    with Path(summary_path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            n = int(rec["N"])
            rows.append(dict(case=case, variant="", dim=rec["Dimension"],
                             N=rec["N"], h=f"{1.0 / n:g}", field="psi",
                             L1=rec["activation_L1"], L2=rec["activation_L2"],
                             Linf=rec["activation_Linf"]))
    for n, dat in (extra_2d_dats or []):
        l1, l2, li = _parse_activation_dat(dat)
        rows.append(dict(case=case, variant="", dim="2D", N=str(n),
                         h=f"{1.0 / int(n):g}", field="psi",
                         L1=l1, L2=l2, Linf=li))
    return rows


def from_eikonal_ecg(aggregate_path, case="eikonal"):
    rows = []
    with Path(aggregate_path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            n = int(rec["N"])
            common = dict(case=case, variant="", dim=rec["Dimension"],
                          N=rec["N"], h=f"{1.0 / n:g}")
            rows.append({**common, "field": "Phi_e_max",
                         "L1": rec.get("max_L1_err_ref", ""),
                         "L2": rec.get("max_L2_err_ref", ""),
                         "Linf": rec.get("max_Linf_err_ref", "")})
            rows.append({**common, "field": "Phi_e_mean",
                         "L1": rec.get("mean_L1_err_ref", ""),
                         "L2": rec.get("mean_L2_err_ref", ""),
                         "Linf": rec.get("mean_Linf_err_ref", "")})
    return rows


_SPATIAL_FILENAME_PATTERN = re.compile(
    r"^(?P<dimension>\dD)_(?P<cells>\d+)_cells_(?P<solver>explicit|implicit)\.dat$"
)
_FIELD_LINE_PATTERN = re.compile(
    r"^(?P<field>\w+)\s+(?P<L1>[-+0-9.eE]+)\s+(?P<L2>[-+0-9.eE]+)\s+(?P<Linf>[-+0-9.eE]+)\s*$"
)
_DX_PATTERN = re.compile(r"Grid spacing \(dx\)\s*=\s*(\S+)")

MONODOMAIN_SPATIAL_FIELDS = ("Vm", "u1", "u2")
BIDOMAIN_FIELDS = ("Vm", "phiE_gauge", "phiI_gauge", "u1", "u2")


def _parse_field_triples(content, allowed_fields):
    """Return {field: (L1,L2,Linf) strings} for lines matching '<field> <L1> <L2> <Linf>'."""
    out = {}
    for line in content.splitlines():
        m = _FIELD_LINE_PATTERN.match(line.strip())
        if not m or m.group("field") not in allowed_fields:
            continue
        out[m.group("field")] = (m.group("L1"), m.group("L2"), m.group("Linf"))
    return out


def _from_spatial_archive(dir_path, allowed_fields, case):
    rows = []
    for path in sorted(Path(dir_path).glob("*.dat")):
        m = _SPATIAL_FILENAME_PATTERN.match(path.name)
        if not m:
            continue
        content = path.read_text(errors="ignore")
        dx_match = _DX_PATTERN.search(content)
        if not dx_match:
            continue
        base = dict(case=case, variant="", dim=m.group("dimension"),
                    N=m.group("cells"), h=dx_match.group(1))
        for field, (l1, l2, linf) in _parse_field_triples(content, allowed_fields).items():
            rows.append({**base, "field": field, "L1": l1, "L2": l2, "Linf": linf})
    return rows


def from_monodomain_spatial_archive(dir_path, case="monodomain-spatial"):
    return _from_spatial_archive(dir_path, MONODOMAIN_SPATIAL_FIELDS, case)


def from_bidomain_archive(dir_path, case="bidomain"):
    return _from_spatial_archive(dir_path, BIDOMAIN_FIELDS, case)


_ECG_SPATIAL_FILENAME_PATTERN = re.compile(
    r"^ECG_(?P<dimension>\dD)_(?P<cells>\d+)_cells_(?P<solver>explicit|implicit)"
    r"(?:_DT[^_]+)?_manufacturedPseudoECGSummary\.dat$"
)


def _parse_ecg_electrode_table(content):
    """Return {'L1_err_ref': [...], 'L2_err_ref': [...], 'Linf_err_ref': [...]} across electrode rows."""
    header = None
    electrode_rows = []
    for line in content.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "Electrode":
            header = parts
            continue
        if header is not None and len(parts) == len(header) and parts[0][:1] == "E":
            electrode_rows.append(parts)
    if header is None:
        return {}
    idx = {name: i for i, name in enumerate(header)}
    out = {}
    for key in ("L1_err_ref", "L2_err_ref", "Linf_err_ref"):
        if key in idx:
            out[key] = [float(row[idx[key]]) for row in electrode_rows]
    return out


# driverFoam's own post-processing (post_processing_manufactured.py) discards 1D/2D
# archived ECG cases as "unsupported": the numerical pseudoECG is accumulated as a 3D
# cell-volume sum, while the 1D/2D manufactured references are lower-dimensional
# integrals, so their errors don't converge under refinement (confirmed against a real
# sweep run 2026-07-17: 1D/2D Phi_e error is flat across N=10..80, not decreasing).
_ECG_SPATIAL_SUPPORTED_DIMENSIONS = ("3D",)


_BATH_STRUCTURED_FIELDS = ("Vm", "phiE", "phiI")


def from_bath_structured(errors_path, case="bath"):
    """bathBidomain/postProcessing/bath_bidomain_errors.csv -> canonical rows.
    Columns: N (or h/dx) + L2_<field> (+ optional L1_/Linf_). h falls back to 1/N.
    NOTE: real header not on disk at plan time; confirm when bathBidomain is run."""
    rows = []
    with Path(errors_path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            n = rec.get("N") or rec.get("cells") or ""
            h = rec.get("h") or rec.get("dx") or (f"{1.0 / int(n):g}" if n else "")
            for field in _BATH_STRUCTURED_FIELDS:
                l2 = rec.get(f"L2_{field}", "")
                if l2 == "":
                    continue
                rows.append(dict(case=case, variant="structured", dim="3D",
                                 N=str(n), h=str(h), field=field,
                                 L1=rec.get(f"L1_{field}", ""), L2=l2,
                                 Linf=rec.get(f"Linf_{field}", "")))
    return rows


_BATH_N_DIR = re.compile(r"N(\d+)", re.IGNORECASE)
# canonical physics identities (confirm names against FINAL_SOLUTION.md):
_BATH_TET_FIELDS = {
    "phiE": "heartPhiE",                          # extracellular potential
    "fluxJump": "x0FluxJump",                     # interface-current continuity
    "intracellularLeak": "x0IntracellularLeak",   # intracellular insulation
}


def from_bath_interface_metrics(study_dir, case="bath_tet"):
    """Glob <study_dir>/N*/bathBidomainInterfaceMetrics.csv; N from the parent dir.
    variant = '<method>/<assembly>'; h = 1/N nominal (tet); L1 unavailable."""
    rows = []
    for csv_path in sorted(Path(study_dir).glob("N*/bathBidomainInterfaceMetrics.csv")):
        m = _BATH_N_DIR.search(csv_path.parent.name)
        if not m:
            continue
        n = int(m.group(1))
        with csv_path.open(newline="") as fh:
            rec = next(csv.DictReader(fh), None)
        if rec is None:
            continue
        variant = f"{rec.get('method', '')}/{rec.get('assembly', '')}"
        base = dict(case=case, variant=variant, dim="3D",
                    N=str(n), h=f"{1.0 / n:g}")
        for field, col in _BATH_TET_FIELDS.items():
            l2 = rec.get(f"{col}_L2", "")
            if l2 == "":
                continue
            rows.append({**base, "field": field, "L1": "",
                         "L2": l2, "Linf": rec.get(f"{col}_Linf", "")})
    return rows


def from_niederer_points(root_dir, case="niederer"):
    """Map Niederer probe activation times into the canonical schema, carrying the
    raw activation time in the L2 slot so check_against_reference compares it by
    relative tolerance. variant=<config dir>, field=<probe label>, h=<DX from dir>.
    Real cached header: Label,Points:0..2,activationTime (seconds)."""
    rows = []
    for csv_path in sorted(Path(root_dir).glob("*/*_points_*.csv")):
        config = csv_path.parent.name
        m = re.search(r"DX([0-9.]+)", config)
        h = m.group(1) if m else ""
        with csv_path.open(newline="") as fh:
            for rec in csv.DictReader(fh):
                probe = rec.get("Label") or rec.get("probe") or rec.get("point") or ""
                t = rec.get("activationTime")
                if t is None:
                    t = rec.get("activationTime_ms") or ""
                if not probe or t == "":
                    continue
                rows.append(dict(case=case, variant=config, dim="3D",
                                 N=str(probe), h=str(h), field=str(probe),
                                 L1="", L2=str(t), Linf=""))
    return rows


def from_bath_parallel_equivalence(comparison_path):
    """Read comparison.csv (metric,serial,parallel,absoluteDifference,tolerance,pass).
    Return (all_pass, failures)."""
    failures = []
    with Path(comparison_path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            if str(rec.get("pass", "")).strip().lower() != "true":
                failures.append(
                    f"{rec.get('metric')}: |delta|={rec.get('absoluteDifference')} "
                    f"tol={rec.get('tolerance')}")
    return (not failures), failures


def from_pseudo_ecg_spatial_archive(dir_path, case="pseudo-ecg-spatial"):
    rows = []
    for path in sorted(Path(dir_path).glob("ECG_*_manufacturedPseudoECGSummary.dat")):
        m = _ECG_SPATIAL_FILENAME_PATTERN.match(path.name)
        if not m or m.group("dimension") not in _ECG_SPATIAL_SUPPORTED_DIMENSIONS:
            continue
        cols = _parse_ecg_electrode_table(path.read_text(errors="ignore"))
        if not cols:
            continue
        n = int(m.group("cells"))
        base = dict(case=case, variant="", dim=m.group("dimension"),
                    N=str(n), h=f"{1.0 / n:g}")
        max_row = {**base, "field": "Phi_e_max"}
        mean_row = {**base, "field": "Phi_e_mean"}
        for out_key, col_key in (("L1", "L1_err_ref"), ("L2", "L2_err_ref"), ("Linf", "Linf_err_ref")):
            values = cols.get(col_key, [])
            max_row[out_key] = f"{max(values):g}" if values else ""
            mean_row[out_key] = f"{sum(values) / len(values):g}" if values else ""
        rows.append(max_row)
        rows.append(mean_row)
    return rows
