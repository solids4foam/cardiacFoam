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


def _mono_tet(root: Path):
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "monodomainPseudoECG", "tetConvergence")
    rows = adapters.from_monodomain_tet_vm(sweep_cases, manifest)
    rows += adapters.from_monodomain_tet_ecg(sweep_cases, manifest)
    return schema.fill_rates(rows)


def _eikonal_tet(root: Path):
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "eikonalECG", "tetConvergence")
    rows = adapters.from_eikonal_tet_activation(sweep_cases, manifest)
    rows += adapters.from_eikonal_tet_ecg(sweep_cases, manifest)
    return schema.fill_rates(rows)


def _coupling1D3D_hex(root: Path):
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


def _sweep_cases_and_manifest(root: Path, tutorial: str, study: str):
    """Every sweep (hex or tet, any tutorial) archives to its own
    setup/studies/<study>/sweepCases/ (openfoam_driver's generic snapshot/
    diff collector) with the driving sweep-run's own --output-dir
    (sweep_manifest.json) at setup/studies/<study>/sweepRun/ -- one shared
    convention regardless of tutorial or mesh family, since a tutorial
    folder can hold more than just hex/tet (see run_study_sweep_common.sh,
    which derives this same path from the sweep.json's own location)."""
    study_dir = root / _TUT / tutorial / "setup" / "studies" / study
    return study_dir / "sweepCases", study_dir / "sweepRun" / "sweep_manifest.json"


def _eikonal_hex(root: Path):
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "eikonalECG", "spatialConvergence")
    excl = root / _TUT / "eikonalECG/excluded_from_pureEikonal_report"
    fname = "2D_{n}_cells_eikonal_manufacturedEikonalActivationTime.dat"
    extra_2d = [(n, excl / fname.format(n=n)) for n in (10, 20, 40, 80)
                if (excl / fname.format(n=n)).exists()]
    rows = adapters.from_eikonal_activation(sweep_cases, manifest, extra_2d_dats=extra_2d)
    rows += adapters.from_eikonal_ecg(sweep_cases, manifest)
    return schema.fill_rates(rows)


def _mono_hex(root: Path):
    # One monodomain sweep produces both the Vm/auxiliary fields and the
    # pseudo-ECG samples, so both are aggregated into a single CSV, matching
    # how the eikonal rows carry their activation field and ECG functional
    # together. The two adapters read the same sweepCases archive and emit
    # disjoint fields: Vm/u1/u2 in 1D-3D, and Phi_e_max/Phi_e_mean in 3D.
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "monodomainPseudoECG", "spatialConvergence")
    rows = adapters.from_monodomain_spatial_archive(sweep_cases, manifest)
    rows += adapters.from_pseudo_ecg_spatial_archive(sweep_cases, manifest)
    return schema.fill_rates(rows)


def _bidomain_hex(root: Path):
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "bidomain", "spatialConvergence")
    return schema.fill_rates(adapters.from_bidomain_archive(sweep_cases, manifest))


def _bidomain_tet(root: Path):
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "bidomain", "tetConvergence")
    return schema.fill_rates(adapters.from_bidomain_tet_archive(sweep_cases, manifest))


def _bath_hex(root: Path):
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "bathBidomain", "spatialConvergence")
    return schema.fill_rates(adapters.from_bath_hex_archive(sweep_cases, manifest))


def _bath_tet(root: Path):
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "bathBidomain", "tetConvergence")
    return schema.fill_rates(adapters.from_bath_interface_metrics(sweep_cases, manifest))


def _niederer_hex(root: Path):
    src = (root / "tutorials/NiedererEtAl2011/NiedererEtAl2011verification"
                  "/setup/cachedCasePostProcessing")
    return schema.fill_rates(adapters.from_niederer_points(src))


CASES = {
    "mono_tet": _mono_tet, "eikonal_tet": _eikonal_tet, "coupling1D3D_hex": _coupling1D3D_hex,
    "eikonal_hex": _eikonal_hex, "mono_hex": _mono_hex,
    "bidomain_hex": _bidomain_hex,
    "bidomain_tet": _bidomain_tet,
    "bath_hex": _bath_hex, "bath_tet": _bath_tet, "niederer_hex": _niederer_hex,
}

_OUT = {
    "mono_tet": "monodomainPseudoECG/setup/results/mono_tet_convergence.csv",
    "eikonal_tet": "eikonalECG/setup/results/eikonal_tet_convergence.csv",
    "coupling1D3D_hex": "monodomain1D3D/setup/results/coupling1D3D_hex_convergence.csv",
    "eikonal_hex": "eikonalECG/setup/results/eikonal_hex_convergence.csv",
    "mono_hex": "monodomainPseudoECG/setup/results/mono_hex_convergence.csv",
    "bidomain_hex": "bidomain/setup/results/bidomain_hex_convergence.csv",
    "bidomain_tet": "bidomain/setup/results/bidomain_tet_convergence.csv",
    "bath_hex": "bathBidomain/setup/results/bath_hex_convergence.csv",
    "bath_tet": "bathBidomain/setup/mesh/tet/results/bath_tet_convergence.csv",
    "niederer_hex": ("../NiedererEtAl2011/NiedererEtAl2011verification"
                 "/setup/results/niederer_hex_activation.csv"),
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
