"""CLI: derive canonical <case>_convergence.csv from a case's native output.

Run as a script (not `-m`):
    python applications/scripts/paperI_results/aggregate.py <tet|coupling|eikonal>
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import adapters  # noqa: E402
import schema    # noqa: E402

_TUT = "tutorials/manufacturedSolutions"


def _mono_tet(root: Path):
    native = root / _TUT / "monodomainPseudoECG/setup/results/scheme_study.csv"
    if native.exists():
        rows = adapters.from_tet_scheme_study(
            native, "tet", {
                "Vm": ("mono_L2", "mono_Linf"),
                "Phi_e": ("ecg_L2", "ecg_Linf"),
            },
        )
        return schema.fill_rates(rows)
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "monodomainPseudoECG", "tetConvergence")
    rows = adapters.from_monodomain_tet_vm(sweep_cases, manifest)
    rows += adapters.from_monodomain_tet_ecg(sweep_cases, manifest)
    return schema.fill_rates(rows)


def _mesh_h_by_nominal_n(path: Path):
    """Read actual tet cell counts and return h_eff=nCells^(-1/3)."""
    records = json.loads(path.read_text())["levels"]
    return {
        int(rec["nominal_N"]): float(rec["n_cells"]) ** (-1.0 / 3.0)
        for rec in records
    }


def _mono_tet_frontal(root: Path):
    study = root / _TUT / "monodomainPseudoECG/setup/studies/tetConvergence"
    archive = study / "sweepCasesOptimised"
    manifest = study / "sweepRunOptimised/sweep_manifest.json"
    h_by_n = _mesh_h_by_nominal_n(study / "mesh_metadata_optimised.json")
    rows = adapters.from_monodomain_tet_vm(
        archive, manifest, case="mono_tet_frontal",
        variant_of=adapters.frontal_monodomain_variant,
        h_by_n=h_by_n, require_completed=True, require_case_dirs=True,
    )
    return schema.fill_rates(rows)


def _eikonal_tet(root: Path):
    native = root / _TUT / "eikonalECG/setup/results/scheme_study.csv"
    if native.exists():
        rows = adapters.from_tet_scheme_study(
            native, "eikonal_tet", {
                "activationTime": ("activationTime_L2", "activationTime_Linf"),
                "Phi_e": ("ecg_L2", "ecg_Linf"),
            },
        )
        schemes = {row["variant"] for row in rows}
        levels = {(row["variant"], row["N"]) for row in rows}
        expected = {(scheme, str(n)) for scheme in ("GaussLinear", "leastSquares")
                    for n in (10, 20, 40, 80)}
        if schemes == {"GaussLinear", "leastSquares"} and levels == expected:
            return schema.fill_rates(rows)
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


def _hex_sweep_cases_and_manifest(root: Path, tutorial: str):
    """Use the current hexConvergence name, retaining test/legacy fallback."""
    current = _sweep_cases_and_manifest(root, tutorial, "hexConvergence")
    if current[1].exists():
        return current
    return _sweep_cases_and_manifest(root, tutorial, "spatialConvergence")


def _eikonal_hex(root: Path):
    sweep_cases, manifest = _hex_sweep_cases_and_manifest(root, "eikonalECG")
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
    sweep_cases, manifest = _hex_sweep_cases_and_manifest(root, "monodomainPseudoECG")
    rows = adapters.from_monodomain_spatial_archive(sweep_cases, manifest)
    rows += adapters.from_pseudo_ecg_spatial_archive(sweep_cases, manifest)
    return schema.fill_rates(rows)


def _bidomain_hex(root: Path):
    sweep_cases, manifest = _hex_sweep_cases_and_manifest(root, "bidomain")
    return schema.fill_rates(adapters.from_bidomain_archive(sweep_cases, manifest))


def _bidomain_tet(root: Path):
    native = root / _TUT / "bidomain/setup/results/scheme_study.csv"
    if native.exists():
        rows = adapters.from_tet_scheme_study(
            native, "bidomain_tet", {
                "Vm": ("vm_L2", "vm_Linf"),
                "Phi_e": ("phiE_L2", "phiE_Linf"),
            },
        )
        return schema.fill_rates(rows)
    sweep_cases, manifest = _sweep_cases_and_manifest(root, "bidomain", "tetConvergence")
    return schema.fill_rates(adapters.from_bidomain_tet_archive(sweep_cases, manifest))


def _bath_hex(root: Path):
    sweep_cases, manifest = _hex_sweep_cases_and_manifest(root, "bathBidomain")
    return schema.fill_rates(adapters.from_bath_hex_archive(sweep_cases, manifest))


def _bath_tet(root: Path):
    base = root / _TUT / "bathBidomain/setup/mesh/tet/studies/coupling/results"
    files = {n: base / f"N{n}_predictor/metrics.csv" for n in (10, 20, 40)}
    return schema.fill_rates(adapters.from_bath_interface_metric_files(files))


def _niederer_hex(root: Path):
    src = (root / "tutorials/NiedererEtAl2011/NiedererEtAl2011verification"
                  "/setup/cachedCasePostProcessing")
    return schema.fill_rates(adapters.from_niederer_points(src))


def _as_experiment(rows, experiment_id):
    """Assign the stable experiment identifier independently of legacy files."""
    return [{**row, "case": experiment_id} for row in rows]


def _normalized(builder, experiment_id):
    return lambda root: _as_experiment(builder(root), experiment_id)


CASES = {
    "mono_tet": _mono_tet, "mono_tet_frontal": _mono_tet_frontal,
    "eikonal_tet": _eikonal_tet, "coupling1D3D_hex": _coupling1D3D_hex,
    "eikonal_hex": _eikonal_hex, "mono_hex": _mono_hex,
    "bidomain_hex": _bidomain_hex,
    "bidomain_tet": _bidomain_tet,
    "bath_hex": _bath_hex, "bath_tet": _bath_tet, "niederer_hex": _niederer_hex,
    "monodomain_cartesian": _normalized(_mono_hex, "monodomain_cartesian"),
    "monodomain_tet_generic": _normalized(_mono_tet, "monodomain_tet_generic"),
    "monodomain_tet_frontal": _normalized(_mono_tet_frontal, "monodomain_tet_frontal"),
    "eikonal_cartesian": _normalized(_eikonal_hex, "eikonal_cartesian"),
    "eikonal_tet_generic": _normalized(_eikonal_tet, "eikonal_tet_generic"),
    "bidomain_cartesian": _normalized(_bidomain_hex, "bidomain_cartesian"),
    "bidomain_tet_generic": _normalized(_bidomain_tet, "bidomain_tet_generic"),
    "bath_bidomain_cartesian": _normalized(_bath_hex, "bath_bidomain_cartesian"),
    "bath_bidomain_tet_conformal": _normalized(_bath_tet, "bath_bidomain_tet_conformal"),
    "purkinje_monodomain_coupled": _normalized(_coupling1D3D_hex, "purkinje_monodomain_coupled"),
    "niederer_cartesian": _normalized(_niederer_hex, "niederer_cartesian"),
}

_OUT = {
    "mono_tet": "monodomainPseudoECG/setup/results/mono_tet_convergence.csv",
    "mono_tet_frontal": "monodomainPseudoECG/setup/results/monodomain_tet_optimised_convergence.csv",
    "eikonal_tet": "eikonalECG/setup/results/eikonal_tet_convergence.csv",
    "coupling1D3D_hex": "monodomain1D3D/setup/results/coupling1D3D_hex_convergence.csv",
    "eikonal_hex": "eikonalECG/setup/results/eikonal_hex_convergence.csv",
    "mono_hex": "monodomainPseudoECG/setup/results/mono_hex_convergence.csv",
    "bidomain_hex": "bidomain/setup/results/bidomain_hex_convergence.csv",
    "bidomain_tet": "bidomain/setup/results/bidomain_tet_convergence.csv",
    "bath_hex": "bathBidomain/setup/results/bath_hex_convergence.csv",
    "bath_tet": "bathBidomain/setup/mesh/tet/results/bath_tet_reported_convergence.csv",
    "niederer_hex": ("../NiedererEtAl2011/NiedererEtAl2011verification"
                 "/setup/results/niederer_hex_activation.csv"),
    "monodomain_cartesian": "monodomainPseudoECG/setup/results/monodomain_cartesian.csv",
    "monodomain_tet_generic": "monodomainPseudoECG/setup/results/monodomain_tet_generic.csv",
    "monodomain_tet_frontal": "monodomainPseudoECG/setup/results/monodomain_tet_frontal.csv",
    "eikonal_cartesian": "eikonalECG/setup/results/eikonal_cartesian.csv",
    "eikonal_tet_generic": "eikonalECG/setup/results/eikonal_tet_generic.csv",
    "bidomain_cartesian": "bidomain/setup/results/bidomain_cartesian.csv",
    "bidomain_tet_generic": "bidomain/setup/results/bidomain_tet_generic.csv",
    "bath_bidomain_cartesian": "bathBidomain/setup/results/bath_bidomain_cartesian.csv",
    "bath_bidomain_tet_conformal": "bathBidomain/setup/results/bath_bidomain_tet_conformal.csv",
    "purkinje_monodomain_coupled": "monodomain1D3D/setup/results/purkinje_monodomain_coupled.csv",
    "niederer_cartesian": ("../NiedererEtAl2011/NiedererEtAl2011verification"
                 "/setup/results/niederer_cartesian.csv"),
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
