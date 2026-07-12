"""Read each conforming case's native output into canonical rows (no rates yet).

Each function returns a list of dicts carrying the non-rate canonical keys
(case, variant, dim, N, h, field, L1, L2, Linf); rates are applied later by
schema.fill_rates. Native files are only read, never rewritten.
"""
from __future__ import annotations

import csv
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
