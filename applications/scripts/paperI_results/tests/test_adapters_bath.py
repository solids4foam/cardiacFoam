import adapters

# ---- Task 3: structured bath ------------------------------------------------

def test_bath_structured_three_fields(tmp_path):
    p = tmp_path / "bath_bidomain_errors.csv"
    p.write_text(
        "N,h,L2_Vm,L2_phiE,L2_phiI\n"
        "10,0.1,4e-3,5e-3,6e-3\n"
        "20,0.05,1e-3,1.25e-3,1.5e-3\n"
    )
    rows = adapters.from_bath_structured(p)
    assert {r["field"] for r in rows} == {"Vm", "phiE", "phiI"}
    vm20 = next(r for r in rows if r["field"] == "Vm" and r["N"] == "20")
    assert vm20["variant"] == "structured" and vm20["dim"] == "3D"
    assert vm20["L2"] == "1e-3" and vm20["h"] == "0.05"

def test_bath_structured_h_from_N_when_absent(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text("N,L2_Vm\n10,4e-3\n")
    rows = adapters.from_bath_structured(p)
    assert rows[0]["h"] == "0.1"   # 1/10 nominal

def test_bath_structured_preserves_dimension_when_present(tmp_path):
    p = tmp_path / "bath_bidomain_errors.csv"
    p.write_text(
        "Dimension,N,L2_Vm\n"
        "1D,10,4e-3\n"
        "3D,10,5e-3\n"
    )
    rows = adapters.from_bath_structured(p)
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
    for n, l2 in (("10", "4e-3"), ("20", "1e-3")):
        d = tmp_path / f"N{n}"; d.mkdir()
        _write_metrics(d / "bathBidomainInterfaceMetrics.csv",
                       _metrics_row(heartPhiE_L2=l2))
    rows = adapters.from_bath_interface_metrics(tmp_path)
    assert {r["field"] for r in rows} == {
        "heartPhiE", "bathPhiE", "x0FluxJump", "x1FluxJump",
        "x0IntracellularLeak", "x1IntracellularLeak",
        "x0AssembledFlux", "x1AssembledFlux"}
    phie20 = next(r for r in rows if r["field"] == "heartPhiE" and r["N"] == "20")
    assert phie20["variant"] == "distanceWeightedHarmonic/matchedSubmesh"
    assert phie20["dim"] == "3D" and phie20["h"] == "0.05" and phie20["L2"] == "1e-3"

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
