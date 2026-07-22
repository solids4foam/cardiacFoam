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
    src = root / _TUT / "monodomainPseudoECG/setup/results/scheme_study.csv"
    return schema.fill_rates(adapters.from_tet_scheme_study(src))


def _eikonal_tet(root: Path):
    src = root / _TUT / "eikonalECG/setup/results/scheme_study.csv"
    return schema.fill_rates(adapters.from_eikonal_tet_scheme_study(src))


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


def _mono_spatial(root: Path):
    src = root / _TUT / "monodomainPseudoECG/driverPostProcessingArchive_postProcessing"
    return schema.fill_rates(adapters.from_monodomain_spatial_archive(src))


def _pseudo_ecg_spatial(root: Path):
    src = root / _TUT / "monodomainPseudoECG/driverPostProcessingArchive_postProcessing"
    return schema.fill_rates(adapters.from_pseudo_ecg_spatial_archive(src))


def _bidomain(root: Path):
    src = root / _TUT / "bidomain/driverPostProcessingArchive_postProcessing"
    return schema.fill_rates(adapters.from_bidomain_archive(src))


def _bidomain_tet(root: Path):
    src = root / _TUT / "bidomain/setup/results/scheme_study.csv"
    return schema.fill_rates(adapters.from_bidomain_tet_scheme_study(src))


def _bath(root: Path):
    src = root / _TUT / "bathBidomain/postProcessing/bath_bidomain_errors.csv"
    return schema.fill_rates(adapters.from_bath_structured(src))


def _bath_tet(root: Path):
    # run_parallel_interface_sweep.sh (matchedSubmesh / distanceWeightedHarmonic)
    # writes N10/20/40/80 here, one bathBidomainInterfaceMetrics.csv per N.
    src = (root / _TUT / "bathBidomain/setup/mesh/tet"
                 "/interfaceStudy/matchedSubmesh/distanceWeightedHarmonic")
    return schema.fill_rates(adapters.from_bath_interface_metrics(src))


def _niederer(root: Path):
    src = (root / "tutorials/NiedererEtAl2011/NiedererEtAl2011verification"
                  "/setup/cachedCasePostProcessing")
    return schema.fill_rates(adapters.from_niederer_points(src))


CASES = {
    "tet": _tet, "eikonal_tet": _eikonal_tet, "coupling": _coupling,
    "eikonal": _eikonal, "mono_spatial": _mono_spatial,
    "pseudo_ecg_spatial": _pseudo_ecg_spatial, "bidomain": _bidomain,
    "bidomain_tet": _bidomain_tet,
    "bath": _bath, "bath_tet": _bath_tet, "niederer": _niederer,
}

_OUT = {
    "tet": "monodomainPseudoECG/setup/results/tet_convergence.csv",
    "eikonal_tet": "eikonalECG/setup/results/eikonal_tet_convergence.csv",
    "coupling": "monodomain1D3D/setup/results/coupling_convergence.csv",
    "eikonal": "eikonalECG/setup/results/eikonal_convergence.csv",
    "mono_spatial": "monodomainPseudoECG/setup/results/mono_spatial_convergence.csv",
    "pseudo_ecg_spatial": "monodomainPseudoECG/setup/results/pseudo_ecg_spatial_convergence.csv",
    "bidomain": "bidomain/setup/results/bidomain_convergence.csv",
    "bidomain_tet": "bidomain/setup/results/bidomain_tet_convergence.csv",
    "bath": "bathBidomain/setup/results/bath_convergence.csv",
    "bath_tet": "bathBidomain/setup/mesh/tet/results/bath_tet_convergence.csv",
    "niederer": ("../NiedererEtAl2011/NiedererEtAl2011verification"
                 "/setup/results/niederer_activation.csv"),
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
