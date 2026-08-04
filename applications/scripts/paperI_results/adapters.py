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


def _parse_field_triples(content, allowed_fields):
    """Return {field: (L1,L2,Linf) strings} for lines matching '<field> <L1> <L2> <Linf>'."""
    out = {}
    for line in content.splitlines():
        m = _FIELD_LINE_PATTERN.match(line.strip())
        if not m or m.group("field") not in allowed_fields:
            continue
        out[m.group("field")] = (m.group("L1"), m.group("L2"), m.group("Linf"))
    return out


def _load_manifest_entries(manifest_path):
    """Load complete case entries from one sweep manifest."""
    manifest = json.loads(Path(manifest_path).read_text())
    return manifest["cases"]


def _sweep_sources(sweep_cases_dir, manifest_path):
    """Return paired archive/manifest sources for one or several run batches."""
    dirs = ([Path(sweep_cases_dir)] if isinstance(sweep_cases_dir, (str, Path))
            else [Path(p) for p in sweep_cases_dir])
    manifests = ([Path(manifest_path)] if isinstance(manifest_path, (str, Path))
                 else [Path(p) for p in manifest_path])
    if len(dirs) == 1 and len(manifests) > 1:
        dirs *= len(manifests)
    if len(dirs) != len(manifests):
        raise ValueError(
            f"sweep archive/manifest count mismatch: {len(dirs)} != {len(manifests)}"
        )
    return list(zip(dirs, manifests))


def _scalar(value):
    """Unwrap a possibly-singleton-list axis value (zip-mode sweep convention,
    e.g. number_cells=[10]) to a plain scalar."""
    if isinstance(value, list) and len(value) == 1:
        return value[0]
    return value


def _iter_sweep_case_dirs(
    sweep_cases_dir, manifest_path, dim_axis, n_axis, *, fixed_dim=None,
    require_completed=False, require_case_dirs=False,
):
    """Yield (case_dir, dim, N, resolved_axis_values) for every case in the
    manifest whose sweepCases/<case_id>/ subfolder actually exists on disk.

    fixed_dim overrides dim_axis lookup entirely -- some studies (every tet
    convergence study so far) fix dimensions=["3D"] in the sweep.json's base
    rather than sweeping it as an independent axis, so it is simply absent
    from resolved_axis_values (which only ever carries independent+dependent
    axis values, never base). Pass fixed_dim explicitly rather than silently
    falling back within this function, so a genuinely missing dim elsewhere
    still means 'skip', not 'assume 3D'."""
    seen = set()
    for archive_dir, source_manifest in _sweep_sources(sweep_cases_dir, manifest_path):
        for entry in _load_manifest_entries(source_manifest):
            case_id = entry["case_id"]
            if case_id in seen:
                raise ValueError(f"duplicate case_id across sweep manifests: {case_id}")
            seen.add(case_id)
            if require_completed and entry.get("status") != "completed":
                raise ValueError(
                    f"paper sweep case is not completed: {case_id} "
                    f"({entry.get('status', 'missing status')})"
                )
            values = entry["resolved_axis_values"]
            case_dir = archive_dir / case_id
            if not case_dir.is_dir():
                if require_case_dirs:
                    raise FileNotFoundError(f"missing archived sweep case: {case_dir}")
                continue
            dim = fixed_dim if fixed_dim is not None else _scalar(values.get(dim_axis))
            n = _scalar(values.get(n_axis))
            if dim is None or n is None:
                continue
            yield case_dir, str(dim), int(n), values


_GRID_SPACING_PATTERN = re.compile(r"Grid spacing \(dx\)\s*=\s*([0-9.eE+-]+)")
_NUMBER_OF_CELLS_PATTERN = re.compile(r"Number of cells\s*=\s*(\d+)")


def _measured_h_from_grid_spacing(content):
    """bidomain/monodomain tet verifiers print their own measured
    'Grid spacing (dx) = <value>' line (src/verificationModels/
    verificationUtils.H's structuredManufacturedDx), computed from the
    ACTUAL tet cell count -- not 1/N. A tet mesh's real cell count at a
    given gmsh characteristic length is not exactly N^3, so nominal 1/N
    (correct for hex block meshes by construction) is simply wrong here."""
    m = _GRID_SPACING_PATTERN.search(content)
    return f"{float(m.group(1)):g}" if m else None


def _measured_h_from_cell_count(content):
    """eikonal's activation-time verifier prints only 'Number of cells',
    not a derived Grid spacing -- replicate verificationUtils.H's
    structuredCellsPerDirection/structuredManufacturedDx formula exactly
    (max(1, int(totalCells**(1/3) + 0.5)), then 1/that) rather than reading
    a value the file never writes. Confirmed byte-for-byte against the
    original run_eikonal_tet.sh's own awk: dx=1.0/int((n)^(1/3)+0.5)."""
    m = _NUMBER_OF_CELLS_PATTERN.search(content)
    if not m:
        return None
    n_per_direction = max(1, int(float(m.group(1)) ** (1.0 / 3.0) + 0.5))
    return f"{1.0 / n_per_direction:g}"


_GRAD_SCHEME_VARIANT_LABELS = {"gauss_linear": "GaussLinear", "least_squares": "leastSquares"}


def _grad_scheme_variant(values):
    """Map a tet case's own grad_scheme axis value to the display variant
    label used throughout the committed tet references (a data label
    carried over from the original bash studies' own CSV column, distinct
    from the literal OpenFOAM dict token -- see _GRAD_SCHEME_TOKENS in
    manufactured_fda.py)."""
    return _GRAD_SCHEME_VARIANT_LABELS[_scalar(values["grad_scheme"])]


def _generic_eikonal_variant(values):
    """Select the axis-aligned, non-advection generic eikonal experiment."""
    if values.get("conductivity_label", "axis") != "axis":
        return "__exclude__"
    if str(values.get("eikonal_advection_diffusion_approach", "false")).lower() != "false":
        return "__exclude__"
    return _grad_scheme_variant(values)


_MONODOMAIN_TENSOR_LABELS = {
    "manufacturedFDAMonodomainVerifier": "diagonal",
    "manufacturedAnisotropicMonodomainVerifier": "rotated",
}


def frontal_monodomain_variant(values):
    """Canonical ``tensor/gradient`` label for the Frontal MMS matrix."""
    model = _scalar(values["verification_model_type"])
    try:
        tensor = _MONODOMAIN_TENSOR_LABELS[model]
    except KeyError as exc:
        raise ValueError(f"unsupported Frontal monodomain verifier: {model}") from exc
    return f"{tensor}/{_grad_scheme_variant(values)}"


def from_sweep_cases_field_triples(
    sweep_cases_dir, manifest_path, *, filename_glob, allowed_fields, case,
    dim_axis="dimensions", n_axis="number_cells", variant="",
    fixed_dim=None, measured_h=None, variant_of=None, h_by_n=None,
    require_completed=False, require_case_dirs=False,
):
    """Generic reader for the '<field> L1 L2 Linf' .dat format shared by the
    manufactured verifiers (bidomain/monodomain/bath/eikonal-activation).
    N/dimension come from that case's own sweep_manifest.json
    resolved_axis_values -- never from filename or file content, since not
    every verifier's raw output name is case-parameter-qualified (eikonal's
    activation/ECG summaries use a fixed name; only the per-case
    sweepCases/<case_id>/ subfolder disambiguates them).

    measured_h, when given, parses the real per-case h out of each file's
    own content (tet studies) instead of using the nominal 1/N (correct
    for hex block meshes only). variant_of, when given, derives the
    variant label per-case from resolved_axis_values (tet's grad_scheme)
    instead of the single fixed `variant` string (hex's convention).

    L1 is always kept, for every caller -- the raw verifier always computes
    it. An earlier version of this reader discarded it for tet studies to
    match the original bash pipeline's incomplete aggregation (which never
    extracted an L1 column), but the statistic set a reader reports must be
    uniform across hex/tet regardless of what the old pipeline happened to
    capture; only the underlying field variables are allowed to differ."""
    rows = []
    for case_dir, dim, n, values in _iter_sweep_case_dirs(
        sweep_cases_dir, manifest_path, dim_axis, n_axis, fixed_dim=fixed_dim,
        require_completed=require_completed, require_case_dirs=require_case_dirs,
    ):
        case_variant = variant_of(values) if variant_of is not None else variant
        nominal_h = f"{1.0 / n:g}"
        for path in sorted(case_dir.glob(filename_glob)):
            content = path.read_text(errors="ignore")
            h = None if h_by_n is None else h_by_n.get(n)
            h = h or (measured_h(content) if measured_h is not None else None) or nominal_h
            base = dict(case=case, variant=case_variant, dim=dim, N=str(n), h=h)
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


def _parse_ecg_electrode_table_by_name(content):
    """Return {electrode: {'L1_err_ref': v, 'L2_err_ref': v, 'Linf_err_ref': v}}.

    The aggregate columns collapse the electrodes to max/mean/min, which hides
    the fact that the reported maximum is set by whichever electrode sits
    closest to the source region and therefore carries the sharpest lead-field
    kernel.  Keeping the names lets that be checked rather than assumed.
    """
    header = None
    out = {}
    for line in content.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "Electrode":
            header = parts
            continue
        if header is not None and len(parts) == len(header) and parts[0][:1] == "E":
            idx = {name: i for i, name in enumerate(header)}
            out[parts[0]] = {
                key: float(parts[idx[key]])
                for key in ("L1_err_ref", "L2_err_ref", "Linf_err_ref")
                if key in idx
            }
    return out


def from_sweep_cases_electrode_table(
    sweep_cases_dir, manifest_path, *, filename_glob, case,
    dim_axis="dimensions", n_axis="number_cells", allowed_dims=None,
    fixed_dim=None, variant_of=None, h_source_glob=None, h_parser=None,
    field_name="Phi_e", per_electrode=False,
):
    """Generic reader for the electrode-table ECG summary format (pseudo-ECG,
    eikonal ECG, tet ECG summaries) shared with from_sweep_cases_field_triples:
    N/dimension come from that case's own sweep_manifest.json
    resolved_axis_values, never from filename or file content. allowed_dims
    restricts which dimensions get reported at all (pseudo-ECG only supports
    3D -- see its caller).

    Reports max/mean/min across electrodes uniformly for every caller --
    hex and tet studies alike -- even though the underlying field variable
    differs by solver; the statistic set must not (see this reader's own
    history: tet originally only tracked a running max, matching its
    original bash script's awk, but that was an incompleteness inherited
    from the old pipeline, not a real convention worth preserving).

    variant_of/h_source_glob/h_parser mirror from_sweep_cases_field_triples'
    tet support: variant_of derives the variant label per-case (tet's
    grad_scheme) instead of the fixed empty string hex uses; h_source_glob/
    h_parser read the real measured h from a companion file in the same
    case_dir when the ECG summary itself doesn't report one (neither tet
    ECG summary file prints its own cell count or grid spacing)."""
    rows = []
    for case_dir, dim, n, values in _iter_sweep_case_dirs(
        sweep_cases_dir, manifest_path, dim_axis, n_axis, fixed_dim=fixed_dim
    ):
        if allowed_dims is not None and dim not in allowed_dims:
            continue
        variant = variant_of(values) if variant_of is not None else ""
        h = None
        if h_source_glob is not None:
            for h_path in sorted(case_dir.glob(h_source_glob)):
                h = h_parser(h_path.read_text(errors="ignore"))
                if h is not None:
                    break
        h = h or f"{1.0 / n:g}"
        base = dict(case=case, variant=variant, dim=dim, N=str(n), h=h)
        for path in sorted(case_dir.glob(filename_glob)):
            cols = _parse_ecg_electrode_table(path.read_text(errors="ignore"))
            if not cols:
                continue
            max_row = {**base, "field": f"{field_name}_max"}
            mean_row = {**base, "field": f"{field_name}_mean"}
            min_row = {**base, "field": f"{field_name}_min"}
            for out_key, col_key in (("L1", "L1_err_ref"), ("L2", "L2_err_ref"), ("Linf", "Linf_err_ref")):
                col_values = cols.get(col_key, [])
                max_row[out_key] = f"{max(col_values):g}" if col_values else ""
                mean_row[out_key] = f"{sum(col_values) / len(col_values):g}" if col_values else ""
                min_row[out_key] = f"{min(col_values):g}" if col_values else ""
            rows.append(max_row)
            rows.append(mean_row)
            rows.append(min_row)

            if per_electrode:
                by_name = _parse_ecg_electrode_table_by_name(
                    path.read_text(errors="ignore")
                )
                for electrode in sorted(by_name):
                    vals = by_name[electrode]
                    row = {**base, "field": f"{field_name}_{electrode}"}
                    for out_key, col_key in (
                        ("L1", "L1_err_ref"),
                        ("L2", "L2_err_ref"),
                        ("Linf", "Linf_err_ref"),
                    ):
                        row[out_key] = (
                            f"{vals[col_key]:g}" if col_key in vals else ""
                        )
                    rows.append(row)
    return rows


def from_eikonal_tet_ecg(sweep_cases_dir, manifest_path, case="eikonal_tet"):
    rows = from_sweep_cases_electrode_table(
        sweep_cases_dir, manifest_path,
        filename_glob="manufacturedEikonalECGSummary.dat", case=case,
        fixed_dim="3D", variant_of=_generic_eikonal_variant,
        h_source_glob="manufacturedEikonalActivationTime.dat",
        h_parser=_measured_h_from_cell_count,
    )
    # The tetrahedral convergence observable is the worst electrode, matching
    # the original study definition. Keep one stable Phi_e row per case; the
    # max/mean/min expansion is useful for exploratory electrode analysis but
    # is not the convergence experiment represented by this adapter.
    return [{**row, "field": "Phi_e", "L1": ""} for row in rows
            if row["field"] == "Phi_e_max" and row["variant"] != "__exclude__"]


def from_eikonal_tet_activation(sweep_cases_dir, manifest_path, case="eikonal_tet"):
    # Unlike hex's from_eikonal_activation, the tet reference keeps the raw
    # field name "activationTime" -- no "psi" rename here (confirmed against
    # the committed reference CSV).
    rows = from_sweep_cases_field_triples(
        sweep_cases_dir, manifest_path,
        filename_glob="manufacturedEikonalActivationTime.dat",
        allowed_fields=EIKONAL_ACTIVATION_FIELDS, case=case,
        fixed_dim="3D", measured_h=_measured_h_from_cell_count,
        variant_of=_generic_eikonal_variant,
    )
    return [{**row, "L1": ""} for row in rows
            if row["variant"] != "__exclude__"]


_ACTIVATION_SPLIT_PATTERN = re.compile(
    r"^activationTimeSplit\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*$"
)


def _parse_activation_split(content):
    """Return (L2_bulk, L2_boundary, L2_total) strings from the verifier's
    'activationTimeSplit <bulk> <boundary> <total>' line. Only written when
    eikonalSolverCoeffs.verificationModel.writeErrorField is enabled -- see
    manufacturedEikonalVerifier.C's computeBoundaryBulkNorms call, added
    alongside the standalone gradientReconstructionOrder utility's own
    bulk/boundary split so the solved-field decomposition and the
    reconstruction-only decomposition can be compared level by level on the
    same mesh ladder. Returns None if the case's .dat file predates that
    verifier change or writeErrorField was off."""
    for line in content.splitlines():
        m = _ACTIVATION_SPLIT_PATTERN.match(line.strip())
        if m:
            return m.group(1), m.group(2), m.group(3)
    return None


def from_eikonal_bulk_boundary(sweep_cases_dir, manifest_path, case="eikonal_tet_split"):
    """Bulk/boundary L2 decomposition of the SOLVED activation-time error
    (as opposed to eikonal_gradient_tet.csv, which decomposes the GRADIENT
    OPERATOR's own reconstruction error against an exact analytic field, no
    solve involved). Requires the driving sweep to enable
    eikonalSolverCoeffs.verificationModel.writeErrorField -- see
    setup/studies/errorLocalisation/sweep_tet_error_localisation.json.

    Column set (case,variant,dim,N,h,L2_bulk,L2_boundary,L2_total,
    boundary_energy_fraction) does not match schema.CANONICAL_FIELDS --
    this is a diagnostic without a committed numerical reference, not a
    keyset-gated convergence table -- so callers write it directly (see
    aggregate_bulk_boundary.py) rather than through schema.write_canonical."""
    rows = []
    for case_dir, dim, n, values in _iter_sweep_case_dirs(
        sweep_cases_dir, manifest_path, "dimensions", "number_cells", fixed_dim="3D",
    ):
        variant = _grad_scheme_variant(values)
        for path in sorted(case_dir.glob("manufacturedEikonalActivationTime.dat")):
            content = path.read_text(errors="ignore")
            split = _parse_activation_split(content)
            if split is None:
                continue
            l2_bulk, l2_boundary, l2_total = (float(v) for v in split)
            h = _measured_h_from_cell_count(content) or f"{1.0 / n:g}"
            fraction = (l2_boundary / l2_total) ** 2 if l2_total else 0.0
            rows.append({
                "case": case, "variant": variant, "dim": dim, "N": str(n), "h": h,
                "L2_bulk": f"{l2_bulk:g}", "L2_boundary": f"{l2_boundary:g}",
                "L2_total": f"{l2_total:g}",
                "boundary_energy_fraction": f"{fraction:g}",
            })
    return rows


def from_monodomain_tet_ecg(
    sweep_cases_dir, manifest_path, case="tet", per_electrode=False
):
    # per_electrode is opt-in, not the default. The reported tetrahedral
    # pseudo-ECG diagnostic is a maximum over the five electrodes, and knowing
    # which electrode sets it is useful -- but keyset_gate.py requires the
    # fresh and committed reference key sets to be exactly equal, so emitting
    # Phi_e_E1..E5 into the canonical CSV would fail reproduction against the
    # existing reference. Callers that want the breakdown request it and write
    # it to a separate artifact.
    return from_sweep_cases_electrode_table(
        sweep_cases_dir, manifest_path,
        filename_glob="manufacturedPseudoECGSummary.dat", case=case,
        fixed_dim="3D", variant_of=_grad_scheme_variant,
        h_source_glob="*_cells_*.dat",
        h_parser=_measured_h_from_grid_spacing,
        per_electrode=per_electrode,
    )


def from_monodomain_tet_vm(
    sweep_cases_dir, manifest_path, case="tet", *, variant_of=_grad_scheme_variant,
    h_by_n=None, require_completed=False, require_case_dirs=False,
):
    # case="tet" (not "mono_tet") -- confirmed against the committed
    # reference CSV's own "case" column, a naming artifact carried over
    # from the original bash-era scheme_study.csv this replaces.
    return from_sweep_cases_field_triples(
        sweep_cases_dir, manifest_path,
        filename_glob="*_cells_*.dat", allowed_fields=("Vm",), case=case,
        fixed_dim="3D", measured_h=_measured_h_from_grid_spacing, variant_of=variant_of,
        h_by_n=h_by_n, require_completed=require_completed,
        require_case_dirs=require_case_dirs,
    )


def from_bidomain_tet_archive(sweep_cases_dir, manifest_path, case="bidomain_tet"):
    # Only Vm + phiE_gauge (renamed Phi_e) are reported for tet -- the raw
    # verifier also emits phiI_gauge/u1/u2, but the original bash-era
    # scheme_study.csv this replaces only ever extracted these two columns
    # (confirmed via its own awk), and the committed reference has no rows
    # for the others.
    rows = from_sweep_cases_field_triples(
        sweep_cases_dir, manifest_path,
        filename_glob="*_cells_*.dat", allowed_fields=("Vm", "phiE_gauge"), case=case,
        fixed_dim="3D", measured_h=_measured_h_from_grid_spacing, variant_of=_grad_scheme_variant,
    )
    for row in rows:
        if row["field"] == "phiE_gauge":
            row["field"] = "Phi_e"
    return rows


def from_bath_interface_metrics(sweep_cases_dir, manifest_path, case="bath_tet"):
    """bath_tet's own CSV-row format (method,assembly,fieldSource,time,
    interfaceFaces,<field>_L1,<field>_L2,<field>_Linf,...) is written by
    applications/utilities/bathBidomainInterfaceMetrics (a post-hoc pass over
    the reconstructed mesh, run as its own workflow_dag step -- see
    manufactured_fda_bath_bidomain.py's interfaceMetrics step), structurally
    different from the '<field> L1 L2 Linf' text format every other tet
    verifier uses, so it needs its own reader rather than
    from_sweep_cases_field_triples.

    Unlike bidomain/monodomain/eikonal tet, this utility does not report a
    measured Grid spacing -- h is nominal 1/N here (confirmed against the
    committed reference: h=0.1 at N=10, not 0.0588235 like the others)."""
    rows = []
    for case_dir, dim, n, _values in _iter_sweep_case_dirs(
        sweep_cases_dir, manifest_path, "dimensions", "number_cells", fixed_dim="3D"
    ):
        base = dict(case=case, dim=dim, N=str(n), h=f"{1.0 / n:g}")
        for path in sorted(case_dir.glob("bathBidomainInterfaceMetrics.csv")):
            with path.open(newline="") as fh:
                rec = next(csv.DictReader(fh), None)
            if rec is None:
                continue
            variant = f"{rec.get('method', '')}/{rec.get('assembly', '')}"
            for field, col in _BATH_TET_FIELDS.items():
                l2 = rec.get(f"{col}_L2", "")
                if l2 == "":
                    continue
                rows.append({**base, "variant": variant, "field": field,
                             "L1": rec.get(f"{col}_L1", ""),
                             "L2": l2, "Linf": rec.get(f"{col}_Linf", "")})
    return rows


def from_bath_interface_metric_files(files_by_n, case="bath_tet_reported"):
    """Read the reported predictor-corrector interface metrics directly.

    The same-mesh coupling study predates the sweep archive convention but
    writes the identical utility CSV.  Keeping this small adapter makes the
    paper's accepted N=10,20,40 ladder reproducible without inventing a second
    numerical result or copying values from the manuscript table.
    """
    rows = []
    for n, path in sorted(files_by_n.items()):
        with Path(path).open(newline="") as fh:
            rec = next(csv.DictReader(fh), None)
        if rec is None:
            raise ValueError(f"empty bath interface metrics file: {path}")
        base = dict(
            case=case, variant="predictor/distanceWeightedHarmonic/matchedSubmesh",
            dim="3D", N=str(n), h=f"{1.0 / int(n):g}",
        )
        for field, col in _BATH_TET_FIELDS.items():
            l2 = rec.get(f"{col}_L2", "")
            if l2 == "":
                continue
            rows.append({**base, "field": field,
                         "L1": rec.get(f"{col}_L1", ""),
                         "L2": l2, "Linf": rec.get(f"{col}_Linf", "")})
    return rows


def from_tet_scheme_study(path, case, field_columns):
    """Map the native hand-run tetrahedral scheme table to canonical rows.

    ``field_columns`` maps the canonical field name to its native
    ``(<L2 column>, <Linf column>)`` pair.  These tables are what the legacy
    tet runners actually regenerate, so reading them directly avoids a false
    dependency on an obsolete driverFoam sweep archive.
    """
    rows = []
    with Path(path).open(newline="") as fh:
        for rec in csv.DictReader(fh):
            for field, (l2_col, linf_col) in field_columns.items():
                l2 = rec.get(l2_col, "")
                if l2 == "":
                    continue
                rows.append({
                    "case": case,
                    "variant": rec["scheme"],
                    "dim": "3D",
                    "N": rec["N"],
                    "h": rec["dx"],
                    "field": field,
                    "L1": "",
                    "L2": l2,
                    "Linf": rec.get(linf_col, ""),
                })
    return rows


# driverFoam's own post-processing (post_processing_manufactured.py) discards 1D/2D
# archived ECG cases as "unsupported": the numerical pseudoECG is accumulated as a 3D
# cell-volume sum, while the 1D/2D manufactured references are lower-dimensional
# integrals, so their errors don't converge under refinement (confirmed against a real
# sweep run 2026-07-17: 1D/2D Phi_e error is flat across N=10..80, not decreasing).
_ECG_SPATIAL_SUPPORTED_DIMENSIONS = ("3D",)


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
