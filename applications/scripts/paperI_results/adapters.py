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
