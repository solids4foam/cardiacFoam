import json

import adapters

# ---- Task 3: structured bath (sweepCases archive) ---------------------------


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


def test_bath_hex_archive_three_fields(tmp_path):
    # manufacturedFDABathBidomainVerifier.C writes
    # bathBidomain_<dim>_<N>_cells_<algo>.dat directly (confirmed in src/) --
    # N/dimension still come from the sweep's own manifest, matching every
    # other hex reader's pattern, not from this filename's own embedded values.
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {
        "10_3D": {"dimensions": ["3D"], "number_cells": [10]},
        "20_3D": {"dimensions": ["3D"], "number_cells": [20]},
    })
    for case_id, l2 in (("10_3D", "4e-3"), ("20_3D", "1e-3")):
        case_dir = sweep_cases / case_id
        case_dir.mkdir(parents=True)
        (case_dir / f"bathBidomain_3D_{case_id.split('_')[0]}_cells_implicit.dat").write_text(
            "# Bath-bidomain manufactured solution error summary\n"
            "# dimension 3D\n"
            f"Vm {l2} {l2} {l2}\n"
            f"phiE {l2} {l2} {l2}\n"
            f"phiI {l2} {l2} {l2}\n"
            "u1 9e-9 9e-9 9e-9\n"
            "u2 9e-9 9e-9 9e-9\n"
            "u3 9e-9 9e-9 9e-9\n"
        )
    rows = adapters.from_bath_hex_archive(sweep_cases, manifest)
    assert {r["field"] for r in rows} == {"Vm", "phiE", "phiI"}
    vm20 = next(r for r in rows if r["field"] == "Vm" and r["N"] == "20")
    assert vm20["variant"] == "structured"
    assert vm20["dim"] == "3D" and vm20["L2"] == "1e-3" and vm20["h"] == f"{1.0 / 20:g}"


def test_bath_hex_archive_preserves_dimension_when_present(tmp_path):
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {
        "10_1D": {"dimensions": ["1D"], "number_cells": [10]},
        "10_3D": {"dimensions": ["3D"], "number_cells": [10]},
    })
    for case_id, dim in (("10_1D", "1D"), ("10_3D", "3D")):
        case_dir = sweep_cases / case_id
        case_dir.mkdir(parents=True)
        (case_dir / f"bathBidomain_{dim}_10_cells_implicit.dat").write_text(
            f"# dimension {dim}\nVm 4e-3 4e-3 4e-3\nphiE 5e-3 5e-3 5e-3\nphiI 6e-3 6e-3 6e-3\n"
        )
    rows = adapters.from_bath_hex_archive(sweep_cases, manifest)
    assert {r["dim"] for r in rows} == {"1D", "3D"}

# ---- Task 4: tet bath interface metrics -------------------------------------

def _metrics_row(**over):
    cols = {"method": "distanceWeightedHarmonic", "assembly": "matchedSubmesh",
            "heartPhiE_L2": "5e-3", "heartPhiE_Linf": "7e-3",
            "bathPhiE_L2": "6e-3", "bathPhiE_Linf": "8e-3",
            "x0FluxJump_L2": "1e-3", "x0FluxJump_Linf": "2e-3",
            "x1FluxJump_L2": "1.1e-3", "x1FluxJump_Linf": "2.1e-3",
            "x0IntracellularLeak_L2": "1e-9", "x0IntracellularLeak_Linf": "3e-9",
            "x1IntracellularLeak_L2": "1.1e-9", "x1IntracellularLeak_Linf": "3.1e-9",
            "x0AssembledFlux_L2": "4e-4", "x0AssembledFlux_Linf": "5e-4",
            "x1AssembledFlux_L2": "4.1e-4", "x1AssembledFlux_Linf": "5.1e-4"}
    cols.update(over)
    return cols

def _write_metrics(path, row):
    import csv as _csv
    with open(path, "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=list(row.keys()))
        w.writeheader(); w.writerow(row)

def test_bath_tet_interface_identities(tmp_path):
    # bath_tet is fixed at dimensions=["3D"] in the sweep.json's base (never
    # swept as an independent axis), so it is absent from resolved_axis_values
    # entirely -- fixed_dim="3D" inside from_bath_interface_metrics handles
    # this, same as every other tet reader.
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {
        "10_0.00892857": {"number_cells": [10], "dt_values": [0.00892857]},
        "20_0.00224215": {"number_cells": [20], "dt_values": [0.00224215]},
    })
    for n, l2 in (("10", "4e-3"), ("20", "1e-3")):
        case_id = {"10": "10_0.00892857", "20": "20_0.00224215"}[n]
        d = sweep_cases / case_id
        d.mkdir(parents=True)
        _write_metrics(d / "bathBidomainInterfaceMetrics.csv",
                       _metrics_row(heartPhiE_L2=l2))
    rows = adapters.from_bath_interface_metrics(sweep_cases, manifest)
    assert {r["field"] for r in rows} == {
        "heartPhiE", "bathPhiE", "x0FluxJump", "x1FluxJump",
        "x0IntracellularLeak", "x1IntracellularLeak",
        "x0AssembledFlux", "x1AssembledFlux"}
    phie20 = next(r for r in rows if r["field"] == "heartPhiE" and r["N"] == "20")
    assert phie20["variant"] == "distanceWeightedHarmonic/matchedSubmesh"
    assert phie20["dim"] == "3D" and phie20["h"] == "0.05" and phie20["L2"] == "1e-3"


def test_bath_reported_metric_files_preserve_paper_ladder(tmp_path):
    files = {}
    for n, heart_l2 in ((10, "3e-3"), (20, "1e-3"), (40, "3e-4")):
        path = tmp_path / f"N{n}_predictor.csv"
        _write_metrics(path, _metrics_row(heartPhiE_L2=heart_l2))
        files[n] = path

    rows = adapters.from_bath_interface_metric_files(files)
    heart = [r for r in rows if r["field"] == "heartPhiE"]
    assert [r["N"] for r in heart] == ["10", "20", "40"]
    assert [r["L2"] for r in heart] == ["3e-3", "1e-3", "3e-4"]
    assert all(
        r["variant"] == "predictor/distanceWeightedHarmonic/matchedSubmesh"
        for r in rows
    )

# ---- Task 5: tet bath parallel equivalence ----------------------------------

def test_parallel_equivalence_all_pass(tmp_path):
    p = tmp_path / "comparison.csv"
    p.write_text("metric,serial,parallel,absoluteDifference,tolerance,pass\n"
                 "interfaceFaces,488.0,488.0,0.0,4.9e-06,True\n"
                 "heartPhiE_L1,0.0056,0.0056,0.0,1.6e-10,True\n")
    ok, failures = adapters.from_bath_parallel_equivalence(p)
    assert ok and failures == []

def test_parallel_equivalence_flags_failure(tmp_path):
    p = tmp_path / "comparison.csv"
    p.write_text("metric,serial,parallel,absoluteDifference,tolerance,pass\n"
                 "heartPhiE_L1,0.0056,0.0060,4e-4,1.6e-10,False\n")
    ok, failures = adapters.from_bath_parallel_equivalence(p)
    assert not ok and "heartPhiE_L1" in failures[0]
