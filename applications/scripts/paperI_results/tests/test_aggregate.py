import json

import aggregate


def _write_manifest(path, cases):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        "schema_version": "1.0", "sweep_spec_hash": "x",
        "created_at": "x", "updated_at": "x",
        "cases": [
            {
                "case_id": case_id, "resolved_axis_values": values,
                "override_hash": "x", "run_document_path": "x",
                "workflow_state_path": "x", "status": "completed",
                "outcome": "fresh", "started_at": "x", "updated_at": "x",
            }
            for case_id, values in cases.items()
        ],
    }))


def test_build_tet_rows_have_rates(tmp_path):
    sweep_cases = (
        tmp_path / "tutorials/manufacturedSolutions/monodomainPseudoECG"
        "/setup/studies/tetConvergence/sweepCases"
    )
    manifest = (
        tmp_path / "tutorials/manufacturedSolutions/monodomainPseudoECG"
        "/setup/studies/tetConvergence/sweepRun/sweep_manifest.json"
    )
    _write_manifest(manifest, {
        "least_squares_20": {"grad_scheme": "least_squares", "number_cells": [20]},
        "least_squares_40": {"grad_scheme": "least_squares", "number_cells": [40]},
    })
    (sweep_cases / "least_squares_20").mkdir(parents=True)
    (sweep_cases / "least_squares_20" / "3D_20_cells_implicit.dat").write_text(
        "Vm     1e-3   0.00298745   0.0213467\nGrid spacing (dx)     = 0.030303\n"
    )
    (sweep_cases / "least_squares_40").mkdir(parents=True)
    (sweep_cases / "least_squares_40" / "3D_40_cells_implicit.dat").write_text(
        "Vm     1e-4   0.000711825   0.00598599\nGrid spacing (dx)     = 0.0151515\n"
    )
    rows = aggregate.CASES["mono_tet"](tmp_path)
    fine_vm = [r for r in rows if r["field"] == "Vm" and r["N"] == "40"][0]
    assert fine_vm["rate_L2"] != ""      # rate computed between 20 and 40


def test_build_bidomain_tet_rows_have_rates(tmp_path):
    sweep_cases = (
        tmp_path / "tutorials/manufacturedSolutions/bidomain"
        "/setup/studies/tetConvergence/sweepCases"
    )
    manifest = (
        tmp_path / "tutorials/manufacturedSolutions/bidomain"
        "/setup/studies/tetConvergence/sweepRun/sweep_manifest.json"
    )
    _write_manifest(manifest, {
        "least_squares_20": {"grad_scheme": "least_squares", "number_cells": [20]},
        "least_squares_40": {"grad_scheme": "least_squares", "number_cells": [40]},
    })
    (sweep_cases / "least_squares_20").mkdir(parents=True)
    (sweep_cases / "least_squares_20" / "3D_20_cells_implicit.dat").write_text(
        "Vm     1e-3   0.00184609   0.00203236\n"
        "phiE_gauge     1e-3   0.00298745   0.0213467\n"
        "Grid spacing (dx)     = 0.030303\n"
    )
    (sweep_cases / "least_squares_40").mkdir(parents=True)
    (sweep_cases / "least_squares_40" / "3D_40_cells_implicit.dat").write_text(
        "Vm     1e-4   0.000381311   0.000396357\n"
        "phiE_gauge     1e-4   0.000711825   0.00598599\n"
        "Grid spacing (dx)     = 0.0151515\n"
    )
    rows = aggregate.CASES["bidomain_tet"](tmp_path)
    fine_vm = [r for r in rows if r["field"] == "Vm" and r["N"] == "40"][0]
    assert fine_vm["rate_L2"] != ""      # rate computed between 20 and 40


def test_build_mono_spatial_rows_have_rates(tmp_path):
    sweep_cases = (
        tmp_path / "tutorials/manufacturedSolutions/monodomainPseudoECG"
        "/setup/studies/spatialConvergence/sweepCases"
    )
    manifest = (
        tmp_path / "tutorials/manufacturedSolutions/monodomainPseudoECG"
        "/setup/studies/spatialConvergence/sweepRun/sweep_manifest.json"
    )
    _write_manifest(manifest, {
        "10_3D": {"dimensions": ["3D"], "number_cells": [10]},
        "20_3D": {"dimensions": ["3D"], "number_cells": [20]},
    })
    (sweep_cases / "10_3D").mkdir(parents=True)
    (sweep_cases / "10_3D" / "3D_10_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\nVm     4e-3   4e-3   8e-3\n"
    )
    (sweep_cases / "20_3D").mkdir(parents=True)
    (sweep_cases / "20_3D" / "3D_20_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\nVm     1e-3   1e-3   2e-3\n"
    )
    rows = aggregate.CASES["mono_hex"](tmp_path)
    fine = [r for r in rows if r["N"] == "20"][0]
    assert fine["rate_L2"] == "2.00"


def test_build_bidomain_rows_have_rates(tmp_path):
    sweep_cases = (
        tmp_path / "tutorials/manufacturedSolutions/bidomain"
        "/setup/studies/spatialConvergence/sweepCases"
    )
    manifest = (
        tmp_path / "tutorials/manufacturedSolutions/bidomain"
        "/setup/studies/spatialConvergence/sweepRun/sweep_manifest.json"
    )
    _write_manifest(manifest, {
        "10_3D": {"dimensions": ["3D"], "number_cells": [10]},
        "20_3D": {"dimensions": ["3D"], "number_cells": [20]},
    })
    (sweep_cases / "10_3D").mkdir(parents=True)
    (sweep_cases / "10_3D" / "3D_10_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\nVm     4e-3   4e-3   8e-3\n"
    )
    (sweep_cases / "20_3D").mkdir(parents=True)
    (sweep_cases / "20_3D" / "3D_20_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\nVm     1e-3   1e-3   2e-3\n"
    )
    rows = aggregate.CASES["bidomain_hex"](tmp_path)
    fine = [r for r in rows if r["N"] == "20"][0]
    assert fine["rate_L2"] == "2.00"
