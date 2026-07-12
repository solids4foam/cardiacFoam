"""CLI: derive canonical <case>_convergence.csv from a case's native output.

Run as a script (not `-m`):
    python applications/scripts/paperI_results/aggregate.py <tet|coupling|eikonal>
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import adapters  # noqa: E402
import schema    # noqa: E402

_TUT = "tutorials/manufacturedSolutions"


def _tet(root: Path):
    src = root / _TUT / "monodomainTetMMS/setup/results/scheme_study.csv"
    return schema.fill_rates(adapters.from_tet_scheme_study(src))


def _coupling(root: Path):
    base = root / _TUT / "monodomain1D3D/outputs"
    regimes = {
        "decoupled": "coupled1D3DConvergence_rpvj1e6/coupled_convergence_summary.csv",
        "active": "coupled1D3DConvergence/coupled_convergence_summary.csv",
        "bidirectional": "coupled1D3DConvergence_bidirectional/coupled_convergence_summary.csv",
    }
    rows = []
    for regime, rel in regimes.items():
        p = base / rel
        if p.exists():
            rows += adapters.from_coupling_summary(p, regime=regime)
    return schema.fill_rates(rows)


def _eikonal(root: Path):
    pp = root / _TUT / "eikonalECG/postProcessing"
    excl = root / _TUT / "eikonalECG/excluded_from_pureEikonal_report"
    fname = "2D_{n}_cells_eikonal_manufacturedEikonalActivationTime.dat"
    extra_2d = [(n, excl / fname.format(n=n)) for n in (10, 20, 40, 80)
                if (excl / fname.format(n=n)).exists()]
    rows = adapters.from_eikonal_activation(
        pp / "manufacturedEikonalActivationSummary.csv", extra_2d_dats=extra_2d)
    ecg = pp / "manufacturedEikonalECGAggregateSummary.csv"
    if ecg.exists():
        rows += adapters.from_eikonal_ecg(ecg)
    return schema.fill_rates(rows)


CASES = {"tet": _tet, "coupling": _coupling, "eikonal": _eikonal}

_OUT = {
    "tet": "monodomainTetMMS/setup/results/tet_convergence.csv",
    "coupling": "monodomain1D3D/setup/results/coupling_convergence.csv",
    "eikonal": "eikonalECG/setup/results/eikonal_convergence.csv",
}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=sorted(CASES))
    ap.add_argument("--repo-root", type=Path,
                    default=Path(__file__).resolve().parents[3])
    args = ap.parse_args(argv)
    rows = CASES[args.case](args.repo_root)
    dest = args.repo_root / _TUT / _OUT[args.case]
    schema.write_canonical(dest, rows)
    print(f"wrote {len(rows)} rows -> {dest}")


if __name__ == "__main__":
    main()
