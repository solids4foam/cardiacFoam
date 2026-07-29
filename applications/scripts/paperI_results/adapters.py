"""Read each conforming case's native output into canonical rows (no rates yet).

Each function returns a list of dicts carrying the non-rate canonical keys
(case, variant, dim, N, h, field, L1, L2, Linf); rates are applied later by
schema.fill_rates. Native files are only read, never rewritten.

The hex-sweep readers (bidomain/monodomain/bath/eikonal) share one pattern:
read raw per-case output directly from a sweep-run's sweepCases/<case_id>/
archive (openfoam_driver's generic snapshot/diff collector -- see
core/runtime/output_collection.py), with N/dimension sourced from that same
sweep's sweep_manifest.json (resolved_axis_values), never from filename or
file content. Not every verifier's raw output name is case-parameter-
qualified -- eikonal's activation-time and ECG summaries use a fixed name
regardless of N/dimension (confirmed directly in
src/verificationModels/eikonalVerification/manufacturedEikonalVerifier.C and
src/verificationModels/ecgVerification/eikonalECGManufacturedVerifier.C) --
so the per-case subfolder plus the manifest are the only source that works
uniformly across all four.
"""
from __future__ import annotations

import csv
import json
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


def from_bidomain_tet_scheme_study(path, case="bidomain_tet"):
    """bidomain's setup/mesh/tet/ scheme_study.csv: Vm + phiE_gauge columns
    per scheme/N. Mirrors from_tet_scheme_study; dx -> h, L1 unavailable."""
    rows = []
    with Path(path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            base = dict(case=case, variant=rec["scheme"], dim="3D",
                        N=rec["N"], h=rec["dx"])
            rows.append({**base, "field": "Vm", "L1": "",
                         "L2": rec["vm_L2"], "Linf": rec["vm_Linf"]})
            rows.append({**base, "field": "Phi_e", "L1": "",
                         "L2": rec["phiE_L2"], "Linf": rec["phiE_Linf"]})
    return rows


def from_eikonal_tet_scheme_study(path, case="eikonal_tet"):
    """eikonalECG's tet overlay (setup/mesh/tet, formerly eikonalTetMMS)
    scheme_study.csv: activationTime + ecg columns per scheme/N.
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


_FIELD_LINE_PATTERN = re.compile(
    r"^(?P<field>\w+)\s+(?P<L1>[-+0-9.eE]+)\s+(?P<L2>[-+0-9.eE]+)\s+(?P<Linf>[-+0-9.eE]+)\s*$"
)

MONODOMAIN_SPATIAL_FIELDS = ("Vm", "u1", "u2")
BIDOMAIN_FIELDS = ("Vm", "phiE_gauge", "phiI_gauge", "u1", "u2")
BATH_HEX_FIELDS = ("Vm", "phiE", "phiI")
EIKONAL_ACTIVATION_FIELDS = ("activationTime",)


def _parse_field_triples(content, allowed_fields):
    """Return {field: (L1,L2,Linf) strings} for lines matching '<field> <L1> <L2> <Linf>'."""
    out = {}
    for line in content.splitlines():
        m = _FIELD_LINE_PATTERN.match(line.strip())
        if not m or m.group("field") not in allowed_fields:
            continue
        out[m.group("field")] = (m.group("L1"), m.group("L2"), m.group("Linf"))
    return out


def _load_case_axis_values(manifest_path):
    """case_id -> resolved_axis_values, straight from a sweep-run's own
    sweep_manifest.json (core/runtime/sweep_manifest.py's CaseManifestEntry)."""
    manifest = json.loads(Path(manifest_path).read_text())
    return {entry["case_id"]: entry["resolved_axis_values"] for entry in manifest["cases"]}


def _scalar(value):
    """Unwrap a possibly-singleton-list axis value (zip-mode sweep convention,
    e.g. number_cells=[10]) to a plain scalar."""
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def _iter_sweep_case_dirs(sweep_cases_dir, manifest_path, dim_axis, n_axis):
    """Yield (case_dir, dim, N) for every case in the manifest whose
    sweepCases/<case_id>/ subfolder actually exists on disk."""
    for case_id, values in _load_case_axis_values(manifest_path).items():
        case_dir = Path(sweep_cases_dir) / case_id
        if not case_dir.is_dir():
            continue
        dim = _scalar(values.get(dim_axis))
        n = _scalar(values.get(n_axis))
        if dim is None or n is None:
            continue
        yield case_dir, str(dim), int(n)


def from_sweep_cases_field_triples(
    sweep_cases_dir, manifest_path, *, filename_glob, allowed_fields, case,
    dim_axis="dimensions", n_axis="number_cells", variant="",
):
    """Generic reader for the '<field> L1 L2 Linf' .dat format shared by the
    manufactured verifiers (bidomain/monodomain/bath/eikonal-activation).
    N/dimension come from that case's own sweep_manifest.json
    resolved_axis_values -- never from filename or file content, since not
    every verifier's raw output name is case-parameter-qualified (eikonal's
    activation/ECG summaries use a fixed name; only the per-case
    sweepCases/<case_id>/ subfolder disambiguates them)."""
    rows = []
    for case_dir, dim, n in _iter_sweep_case_dirs(sweep_cases_dir, manifest_path, dim_axis, n_axis):
        base = dict(case=case, variant=variant, dim=dim, N=str(n), h=f"{1.0 / n:g}")
        for path in sorted(case_dir.glob(filename_glob)):
            content = path.read_text(errors="ignore")
            for field, (l1, l2, linf) in _parse_field_triples(content, allowed_fields).items():
                rows.append({**base, "field": field, "L1": l1, "L2": l2, "Linf": linf})
    return rows


def from_eikonal_activation(sweep_cases_dir, manifest_path, extra_2d_dats=None, case="eikonal"):
    rows = from_sweep_cases_field_triples(
        sweep_cases_dir, manifest_path,
        filename_glob="manufacturedEikonalActivationTime.dat",
        allowed_fields=EIKONAL_ACTIVATION_FIELDS, case=case,
    )
    for row in rows:
        row["field"] = "psi"  # canonical schema name; raw file says "activationTime"
    for n, dat in (extra_2d_dats or []):
        l1, l2, li = _parse_activation_dat(dat)
        rows.append(dict(case=case, variant="", dim="2D", N=str(n),
                         h=f"{1.0 / int(n):g}", field="psi",
                         L1=l1, L2=l2, Linf=li))
    return rows


def from_eikonal_ecg(sweep_cases_dir, manifest_path, case="eikonal"):
    return from_sweep_cases_electrode_table(
        sweep_cases_dir, manifest_path,
        filename_glob="manufacturedEikonalECGSummary.dat", case=case,
    )


def from_monodomain_spatial_archive(sweep_cases_dir, manifest_path, case="monodomain-spatial"):
    return from_sweep_cases_field_triples(
        sweep_cases_dir, manifest_path,
        filename_glob="*_cells_*.dat", allowed_fields=MONODOMAIN_SPATIAL_FIELDS, case=case,
    )


def from_bidomain_archive(sweep_cases_dir, manifest_path, case="bidomain"):
    return from_sweep_cases_field_triples(
        sweep_cases_dir, manifest_path,
        filename_glob="*_cells_*.dat", allowed_fields=BIDOMAIN_FIELDS, case=case,
    )


def from_bath_hex_archive(sweep_cases_dir, manifest_path, case="bath"):
    return from_sweep_cases_field_triples(
        sweep_cases_dir, manifest_path,
        filename_glob="bathBidomain_*_cells_*.dat", allowed_fields=BATH_HEX_FIELDS, case=case,
        variant="structured",
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


def from_sweep_cases_electrode_table(
    sweep_cases_dir, manifest_path, *, filename_glob, case,
    dim_axis="dimensions", n_axis="number_cells", allowed_dims=None,
):
    """Generic reader for the electrode-table ECG summary format (pseudo-ECG,
    eikonal ECG) shared with from_sweep_cases_field_triples: N/dimension come
    from that case's own sweep_manifest.json resolved_axis_values, never from
    filename or file content. allowed_dims restricts which dimensions get
    reported at all (pseudo-ECG only supports 3D -- see its caller)."""
    rows = []
    for case_dir, dim, n in _iter_sweep_case_dirs(sweep_cases_dir, manifest_path, dim_axis, n_axis):
        if allowed_dims is not None and dim not in allowed_dims:
            continue
        base = dict(case=case, variant="", dim=dim, N=str(n), h=f"{1.0 / n:g}")
        for path in sorted(case_dir.glob(filename_glob)):
            cols = _parse_ecg_electrode_table(path.read_text(errors="ignore"))
            if not cols:
                continue
            max_row = {**base, "field": "Phi_e_max"}
            mean_row = {**base, "field": "Phi_e_mean"}
            for out_key, col_key in (("L1", "L1_err_ref"), ("L2", "L2_err_ref"), ("Linf", "Linf_err_ref")):
                values = cols.get(col_key, [])
                max_row[out_key] = f"{max(values):g}" if values else ""
                mean_row[out_key] = f"{sum(values) / len(values):g}" if values else ""
            rows.append(max_row)
            rows.append(mean_row)
    return rows


# driverFoam's own post-processing (post_processing_manufactured.py) discards 1D/2D
# archived ECG cases as "unsupported": the numerical pseudoECG is accumulated as a 3D
# cell-volume sum, while the 1D/2D manufactured references are lower-dimensional
# integrals, so their errors don't converge under refinement (confirmed against a real
# sweep run 2026-07-17: 1D/2D Phi_e error is flat across N=10..80, not decreasing).
_ECG_SPATIAL_SUPPORTED_DIMENSIONS = ("3D",)


_BATH_N_DIR = re.compile(r"N(\d+)", re.IGNORECASE)
# The eight diagnostics reported in tbl-bath-bidomain-tet. Each maps a canonical
# field name to its column prefix in bathBidomainInterfaceMetrics.csv; the _L2 /
# _Linf error columns are read from that prefix.
_BATH_TET_FIELDS = {
    "heartPhiE":            "heartPhiE",            # myocardium extracellular potential
    "bathPhiE":            "bathPhiE",              # bath potential
    "x0FluxJump":          "x0FluxJump",            # x=0 interface-current continuity
    "x1FluxJump":          "x1FluxJump",            # x=1 interface-current continuity
    "x0IntracellularLeak": "x0IntracellularLeak",   # x=0 intracellular insulation
    "x1IntracellularLeak": "x1IntracellularLeak",   # x=1 intracellular insulation
    "x0AssembledFlux":     "x0AssembledFlux",       # x=0 assembled-current constitutive error
    "x1AssembledFlux":     "x1AssembledFlux",       # x=1 assembled-current constitutive error
}


def from_bath_interface_metrics(study_dir, case="bath_tet"):
    """Glob <study_dir>/N*/bathBidomainInterfaceMetrics.csv; N from the parent dir.
    variant = '<method>/<assembly>'; h = 1/N nominal (tet)."""
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
            rows.append({**base, "field": field, "L1": rec.get(f"{col}_L1", ""),
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


def from_pseudo_ecg_spatial_archive(sweep_cases_dir, manifest_path, case="pseudo-ecg-spatial"):
    return from_sweep_cases_electrode_table(
        sweep_cases_dir, manifest_path,
        filename_glob="manufacturedPseudoECGSummary.dat", case=case,
        allowed_dims=_ECG_SPATIAL_SUPPORTED_DIMENSIONS,
    )
