import json

import adapters


def _write_manifest(path, cases):
    """cases: {case_id: resolved_axis_values}. Mirrors sweep_manifest.json's
    real shape (core/runtime/sweep_manifest.py) closely enough for the
    adapters, which only read case_id + resolved_axis_values."""
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


def test_bidomain_tet_archive_renames_phie_gauge_and_reads_measured_h(tmp_path):
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    sweep_cases = tmp_path / "sweepCases"
    _write_manifest(manifest, {
        "gauss_linear_10": {"grad_scheme": "gauss_linear", "number_cells": [10]},
        "least_squares_80": {"grad_scheme": "least_squares", "number_cells": [80]},
    })
    (sweep_cases / "gauss_linear_10").mkdir(parents=True)
    (sweep_cases / "gauss_linear_10" / "3D_10_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     1.2e-05   2.3e-05   3.4e-05\n"
        "phiE_gauge     4.5e-05   0.010246   0.055506\n"
        "\nGrid spacing (dx)     = 0.0588235\n"
    )
    (sweep_cases / "least_squares_80").mkdir(parents=True)
    (sweep_cases / "least_squares_80" / "3D_80_cells_implicit.dat").write_text(
        "Vm     1e-06   2e-06   3e-06\n"
        "phiE_gauge     4e-06   0.00269257   0.0191601\n"
        "Grid spacing (dx)     = 0.00757576\n"
    )
    rows = adapters.from_bidomain_tet_archive(sweep_cases, manifest)
    phi = [r for r in rows if r["field"] == "Phi_e" and r["N"] == "10"][0]
    assert phi["variant"] == "GaussLinear" and phi["dim"] == "3D"
    assert phi["h"] == "0.0588235" and phi["L2"] == "0.010246" and phi["Linf"] == "0.055506"
    assert not any(r["field"] == "phiE_gauge" for r in rows)
    phi80 = [r for r in rows if r["field"] == "Phi_e" and r["N"] == "80"][0]
    assert phi80["variant"] == "leastSquares" and phi80["h"] == "0.00757576"


def test_monodomain_tet_vm_and_ecg_share_measured_h_from_grid_spacing(tmp_path):
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    sweep_cases = tmp_path / "sweepCases"
    _write_manifest(manifest, {
        "gauss_linear_10": {"grad_scheme": "gauss_linear", "number_cells": [10]},
    })
    case_dir = sweep_cases / "gauss_linear_10"
    case_dir.mkdir(parents=True)
    (case_dir / "3D_10_cells_implicit.dat").write_text(
        "Vm     1e-05   0.0224588   0.152636\n"
        "Grid spacing (dx)     = 0.0588235\n"
    )
    (case_dir / "manufacturedPseudoECGSummary.dat").write_text(
        "header line 1\nheader line 2\nheader line 3\nheader line 4\nheader line 5\n"
        "Electrode L1_err_ref L2_err_ref Linf_err_ref\n"
        "E1 0.001 0.0038452 0.00426517\n"
        "E2 0.0005 0.002 0.003\n"
    )
    vm_rows = adapters.from_monodomain_tet_vm(sweep_cases, manifest)
    assert vm_rows == [{"case": "tet", "variant": "GaussLinear", "dim": "3D",
                        "N": "10", "h": "0.0588235", "field": "Vm",
                        "L1": "1e-05", "L2": "0.0224588", "Linf": "0.152636"}]
    ecg_rows = adapters.from_monodomain_tet_ecg(sweep_cases, manifest)
    by_field = {r["field"]: r for r in ecg_rows}
    assert set(by_field) == {"Phi_e_max", "Phi_e_mean", "Phi_e_min"}
    assert by_field["Phi_e_max"]["L2"] == "0.0038452" and by_field["Phi_e_max"]["Linf"] == "0.00426517"
    assert by_field["Phi_e_min"]["L2"] == "0.002" and by_field["Phi_e_min"]["Linf"] == "0.003"
    assert by_field["Phi_e_mean"]["L2"] == f"{(0.0038452 + 0.002) / 2:g}"
    assert all(r["h"] == "0.0588235" and r["variant"] == "GaussLinear" for r in ecg_rows)


def test_eikonal_tet_activation_keeps_raw_field_name_and_derives_h_from_cell_count(tmp_path):
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    sweep_cases = tmp_path / "sweepCases"
    _write_manifest(manifest, {
        "least_squares_10": {"grad_scheme": "least_squares", "number_cells": [10]},
    })
    case_dir = sweep_cases / "least_squares_10"
    case_dir.mkdir(parents=True)
    (case_dir / "manufacturedEikonalActivationTime.dat").write_text(
        "# Eikonal manufactured activation-time summary\n"
        "Number of cells = 4913\n"
        "# field L1 L2 Linf\n"
        "activationTime 0.02 0.068 0.30\n"
    )
    rows = adapters.from_eikonal_tet_activation(sweep_cases, manifest)
    assert rows == [{"case": "eikonal_tet", "variant": "leastSquares", "dim": "3D",
                     "N": "10", "h": "0.0588235", "field": "activationTime",
                     "L1": "0.02", "L2": "0.068", "Linf": "0.30"}]


def test_eikonal_tet_ecg_reads_h_from_companion_activation_file(tmp_path):
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    sweep_cases = tmp_path / "sweepCases"
    _write_manifest(manifest, {
        "gauss_linear_10": {"grad_scheme": "gauss_linear", "number_cells": [10]},
    })
    case_dir = sweep_cases / "gauss_linear_10"
    case_dir.mkdir(parents=True)
    (case_dir / "manufacturedEikonalActivationTime.dat").write_text(
        "Number of cells = 4913\nactivationTime 0.02 0.068 0.30\n"
    )
    (case_dir / "manufacturedEikonalECGSummary.dat").write_text(
        "header 1\nheader 2\nheader 3\nheader 4\nheader 5\nheader 6\n"
        "Electrode L1_err_ref L2_err_ref Linf_err_ref\n"
        "E1 0.01 0.0446989 0.0632351\n"
        "E2 0.005 0.02 0.03\n"
    )
    rows = adapters.from_eikonal_tet_ecg(sweep_cases, manifest)
    by_field = {r["field"]: r for r in rows}
    assert set(by_field) == {"Phi_e_max", "Phi_e_mean", "Phi_e_min"}
    assert by_field["Phi_e_max"]["L2"] == "0.0446989" and by_field["Phi_e_max"]["Linf"] == "0.0632351"
    assert by_field["Phi_e_max"]["L1"] == "0.01"
    assert by_field["Phi_e_min"]["L2"] == "0.02" and by_field["Phi_e_min"]["Linf"] == "0.03"
    assert all(r["h"] == "0.0588235" for r in rows)


def test_bath_interface_metrics_reads_from_sweep_cases_with_nominal_h(tmp_path):
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    sweep_cases = tmp_path / "sweepCases"
    _write_manifest(manifest, {"10_0.00892857": {"number_cells": [10], "dt_values": [0.00892857]}})
    case_dir = sweep_cases / "10_0.00892857"
    case_dir.mkdir(parents=True)
    header = ("method,assembly,fieldSource,time,interfaceFaces,"
              "heartPhiE_L1,heartPhiE_L2,heartPhiE_Linf,bathPhiE_L1,bathPhiE_L2,bathPhiE_Linf")
    row = "distanceWeightedHarmonic,matchedSubmesh,numerical,0.196429,120,0.00575671,0.00606607,0.00912926,0.00585127,0.00821426,0.0141732"
    (case_dir / "bathBidomainInterfaceMetrics.csv").write_text(header + "\n" + row + "\n")
    rows = adapters.from_bath_interface_metrics(sweep_cases, manifest)
    heart = [r for r in rows if r["field"] == "heartPhiE"][0]
    assert heart["variant"] == "distanceWeightedHarmonic/matchedSubmesh"
    assert heart["dim"] == "3D" and heart["N"] == "10" and heart["h"] == "0.1"
    assert heart["L2"] == "0.00606607" and heart["Linf"] == "0.00912926"
    bath = [r for r in rows if r["field"] == "bathPhiE"][0]
    assert bath["L2"] == "0.00821426"


def test_coupling_summary(tmp_path):
    p = tmp_path / "coupled_convergence_summary.csv"
    p.write_text(
        "N,h,nodes_1D,L1_3D_Vm,L2_3D_Vm,Linf_3D_Vm,L1_1D_Vm,L2_1D_Vm,Linf_1D_Vm\n"
        "10,0.1,11,0.00132,0.00185713,0.00496035,0.0027,0.00300658,0.00423\n"
    )
    rows = adapters.from_coupling_summary(p, regime="decoupled")
    vm3d = [r for r in rows if r["dim"] == "3D"][0]
    assert vm3d["variant"] == "decoupled" and vm3d["L2"] == "0.00185713"
    assert any(r["dim"] == "1D" and r["L2"] == "0.00300658" for r in rows)


def test_eikonal_activation_reads_sweep_cases_archive(tmp_path):
    # manufacturedEikonalVerifier.C writes a FIXED filename regardless of N/
    # dimension (confirmed directly in src/) -- the sweepCases/<case_id>/
    # subfolder is what disambiguates cases, and N/dimension come from the
    # sweep's own manifest, not from the file's name or content.
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {
        "10_1D": {"dimensions": ["1D"], "number_cells": [10]},
        "10_3D": {"dimensions": ["3D"], "number_cells": [10]},
    })
    (sweep_cases / "10_1D").mkdir(parents=True)
    (sweep_cases / "10_1D" / "manufacturedEikonalActivationTime.dat").write_text(
        "activationTime 0.00158772 0.00215763 0.00486533\n"
    )
    (sweep_cases / "10_3D").mkdir(parents=True)
    (sweep_cases / "10_3D" / "manufacturedEikonalActivationTime.dat").write_text(
        "activationTime 0.00483561 0.00752141 0.0310654\n"
    )

    dat2d = tmp_path / "2D_10_cells_eikonal_manufacturedEikonalActivationTime.dat"
    dat2d.write_text("activationTime 0.00288871 0.00415336 0.0121584\n")

    rows = adapters.from_eikonal_activation(sweep_cases, manifest, extra_2d_dats=[(10, dat2d)])
    r1d = [r for r in rows if r["dim"] == "1D"][0]
    assert r1d["L2"] == "0.00215763" and r1d["field"] == "psi" and r1d["h"] == "0.1"
    r3d = [r for r in rows if r["dim"] == "3D"][0]
    assert r3d["L2"] == "0.00752141"
    r2d = [r for r in rows if r["dim"] == "2D"][0]
    assert r2d["L2"] == "0.00415336" and r2d["h"] == "0.1"


def test_eikonal_ecg_max_and_mean_reads_sweep_cases_archive(tmp_path):
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {"80_3D": {"dimensions": ["3D"], "number_cells": [80]}})
    case_dir = sweep_cases / "80_3D"
    case_dir.mkdir(parents=True)
    (case_dir / "manufacturedEikonalECGSummary.dat").write_text(
        "Manufactured eikonal ECG summary\n"
        "samples 89\n"
        "dimension 3D\n"
        "Electrode  L1_err_ref  L2_err_ref  Linf_err_ref\n"
        "E1 1e-6 1e-6 1.11883e-05\n"
        "E2 1e-6 1e-6 6.65081e-06\n"
    )
    rows = adapters.from_eikonal_ecg(sweep_cases, manifest)
    assert any(r["field"] == "Phi_e_max" and r["Linf"] == "1.11883e-05" for r in rows)
    assert any(r["field"] == "Phi_e_mean" and r["Linf"] == "8.91955e-06" for r in rows)
    assert any(r["field"] == "Phi_e_min" and r["Linf"] == "6.65081e-06" for r in rows)


def test_monodomain_spatial_archive_reads_sweep_cases_archive(tmp_path):
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {"80_3D": {"dimensions": ["3D"], "number_cells": [80]}})
    case_dir = sweep_cases / "80_3D"
    case_dir.mkdir(parents=True)
    (case_dir / "3D_80_cells_implicit.dat").write_text(
        "Manufactured-solution error summary (t = 0.199551):\n"
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     5.54093e-05   7.92099e-05   0.00030695\n"
        "u1     4.7947e-05   6.65211e-05   0.000234554\n"
        "u2     2.08544e-05   2.87229e-05   8.75773e-05\n"
        "-------------------------------------------------\n"
        "\n"
        "Simulation summary:\n"
        "-------------------\n"
        "Number of cells (N)   = 80\n"
        "Solver type           = Implicit\n"
        "Grid spacing (dx)     = 0.0125\n"
        "Time step (dt)        = 0.00224215\n"
    )
    rows = adapters.from_monodomain_spatial_archive(sweep_cases, manifest)
    vm = [r for r in rows if r["field"] == "Vm"][0]
    # dim/N/h come from the manifest's resolved_axis_values, not the file's
    # own "Grid spacing (dx)"/"Number of cells" comments.
    assert vm["dim"] == "3D" and vm["N"] == "80" and vm["h"] == f"{1.0 / 80:g}"
    assert vm["L2"] == "7.92099e-05" and vm["Linf"] == "0.00030695"
    assert {r["field"] for r in rows} == {"Vm", "u1", "u2"}


def test_bidomain_archive_reads_sweep_cases_archive(tmp_path):
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {"20_3D": {"dimensions": ["3D"], "number_cells": [20]}})
    case_dir = sweep_cases / "20_3D"
    case_dir.mkdir(parents=True)
    (case_dir / "3D_20_cells_implicit.dat").write_text(
        "Bidomain manufactured-solution error summary (t = 0.200028):\n"
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm        0.00025015   0.00034074   0.000976693\n"
        "phiE      0.0814577   0.081458   0.0821454\n"
        "phiE_gauge 0.000176892   0.000240922   0.000687749\n"
        "phiI      0.0814536   0.0814537   0.0817234\n"
        "phiI_gauge 7.33107e-05   9.98768e-05   0.000288945\n"
        "u1        3.09526e-05   4.29283e-05   0.000152202\n"
        "u2        1.34773e-05   1.84508e-05   5.38718e-05\n"
        "\n"
        "Number of cells (N)   = 20\n"
        "Solver type           = Implicit\n"
        "Grid spacing (dx)     = 0.05\n"
    )
    rows = adapters.from_bidomain_archive(sweep_cases, manifest)
    assert {r["field"] for r in rows} == {"Vm", "phiE_gauge", "phiI_gauge", "u1", "u2"}
    phie = [r for r in rows if r["field"] == "phiE_gauge"][0]
    assert phie["L2"] == "0.000240922" and phie["h"] == f"{1.0 / 20:g}" and phie["N"] == "20"


def test_pseudo_ecg_spatial_archive_reads_sweep_cases_archive(tmp_path):
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {"80_3D": {"dimensions": ["3D"], "number_cells": [80]}})
    case_dir = sweep_cases / "80_3D"
    case_dir.mkdir(parents=True)
    (case_dir / "manufacturedPseudoECGSummary.dat").write_text(
        "Manufactured pseudo-ECG summary\n"
        "samples 89\n"
        "dimension 3D\n"
        "qChecks 6\n"
        "qReference 96\n"
        "Electrode  L1_err_ref  L2_err_ref  Linf_err_ref  L1_err_q6  L2_err_q6  Linf_err_q6\n"
        "E1 3.10522e-07 3.80622e-07 7.4854e-07 3.10522e-07 3.80622e-07 7.4854e-07\n"
        "E2 1.87535e-07 2.39847e-07 5.00594e-07 1.87535e-07 2.39847e-07 5.00594e-07\n"
    )
    rows = adapters.from_pseudo_ecg_spatial_archive(sweep_cases, manifest)
    maxrow = [r for r in rows if r["field"] == "Phi_e_max"][0]
    meanrow = [r for r in rows if r["field"] == "Phi_e_mean"][0]
    minrow = [r for r in rows if r["field"] == "Phi_e_min"][0]
    assert maxrow["dim"] == "3D" and maxrow["N"] == "80" and maxrow["h"] == f"{1.0 / 80:g}"
    assert float(maxrow["L1"]) == 3.10522e-07
    assert abs(float(meanrow["L1"]) - 2.490285e-07) < 1e-12
    assert float(minrow["L1"]) == 1.87535e-07


def test_pseudo_ecg_spatial_archive_excludes_unsupported_1d_2d(tmp_path):
    sweep_cases = tmp_path / "sweepCases"
    manifest = tmp_path / "sweepRun" / "sweep_manifest.json"
    _write_manifest(manifest, {
        "10_1D": {"dimensions": ["1D"], "number_cells": [10]},
        "10_3D": {"dimensions": ["3D"], "number_cells": [10]},
    })
    header = "Electrode  L1_err_ref  L2_err_ref  Linf_err_ref\n"
    (sweep_cases / "10_1D").mkdir(parents=True)
    (sweep_cases / "10_1D" / "manufacturedPseudoECGSummary.dat").write_text(
        "dimension 1D\n" + header + "E1 6.9e-01 6.9e-01 7.2e-01\n"
    )
    (sweep_cases / "10_3D").mkdir(parents=True)
    (sweep_cases / "10_3D" / "manufacturedPseudoECGSummary.dat").write_text(
        "dimension 3D\n" + header + "E1 3.7e-04 3.7e-04 3.7e-04\n"
    )
    rows = adapters.from_pseudo_ecg_spatial_archive(sweep_cases, manifest)
    assert {r["dim"] for r in rows} == {"3D"}
