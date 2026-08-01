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


def test_build_frontal_monodomain_uses_one_sweep_and_cell_count_h(tmp_path):
    study = (
        tmp_path / "tutorials/manufacturedSolutions/monodomainPseudoECG"
        "/setup/studies/tetConvergence"
    )
    archive = study / "sweepCasesOptimised"
    cases = [
        ("least_squares_40_manufacturedFDAMonodomainVerifier", 40, "0.00051079"),
        ("least_squares_80_manufacturedFDAMonodomainVerifier", 80, "0.000141195"),
    ]
    _write_manifest(study / "sweepRunOptimised/sweep_manifest.json", {
        case_id: {
            "grad_scheme": "least_squares", "number_cells": [n],
            "verification_model_type": "manufacturedFDAMonodomainVerifier",
        }
        for case_id, n, _ in cases
    })
    for case_id, n, l2 in cases:
        case_dir = archive / case_id
        case_dir.mkdir(parents=True)
        (case_dir / f"3D_{n}_cells_implicit.dat").write_text(
            f"Vm 1e-4 {l2} 1e-3\nGrid spacing (dx) = {1/n}\n"
        )
    (study / "mesh_metadata_optimised.json").write_text(json.dumps({
        "h_definition": "n_cells^(-1/3)",
        "levels": [
            {"nominal_N": 40, "n_cells": 349109},
            {"nominal_N": 80, "n_cells": 2662487},
        ],
    }))

    rows = aggregate.CASES["mono_tet_frontal"](tmp_path)
    assert len(rows) == 2
    fine = [r for r in rows if r["N"] == "80"][0]
    assert fine["variant"] == "diagonal/leastSquares"
    assert abs(float(fine["h"]) - 2662487 ** (-1 / 3)) < 1e-15
    assert fine["rate_L2"] == "1.90"


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


def test_build_bath_tet_uses_reported_predictor_files(tmp_path):
    results = (
        tmp_path / "tutorials/manufacturedSolutions/bathBidomain"
        "/setup/mesh/tet/studies/coupling/results"
    )
    for n, heart_l2 in ((10, "3e-3"), (20, "1e-3"), (40, "3e-4")):
        path = results / f"N{n}_predictor" / "metrics.csv"
        path.parent.mkdir(parents=True)
        path.write_text(
            "method,assembly,heartPhiE_L1,heartPhiE_L2,heartPhiE_Linf\n"
            f"distanceWeightedHarmonic,matchedSubmesh,1e-3,{heart_l2},4e-3\n"
        )

    rows = aggregate.CASES["bath_tet"](tmp_path)
    heart = [r for r in rows if r["field"] == "heartPhiE"]
    assert [r["N"] for r in heart] == ["10", "20", "40"]
    assert heart[-1]["rate_L2"] != ""


def test_eikonal_tet_prefers_complete_native_scheme_matrix(tmp_path):
    result = (
        tmp_path / "tutorials/manufacturedSolutions/eikonalECG"
        "/setup/results/scheme_study.csv"
    )
    result.parent.mkdir(parents=True)
    lines = [
        "scheme,N,dx,activationTime_L2,activationTime_Linf,ecg_L2,ecg_Linf,outerIterations"
    ]
    for scheme in ("GaussLinear", "leastSquares"):
        for n in (10, 20, 40, 80):
            lines.append(f"{scheme},{n},{1/n:g},{1/n:g},{2/n:g},{3/n:g},{4/n:g},10")
    result.write_text("\n".join(lines) + "\n")

    rows = aggregate.CASES["eikonal_tet_generic"](tmp_path)
    assert len(rows) == 16
    assert {row["case"] for row in rows} == {"eikonal_tet_generic"}
    assert {row["variant"] for row in rows} == {"GaussLinear", "leastSquares"}
