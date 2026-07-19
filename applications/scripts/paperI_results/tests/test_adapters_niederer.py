import adapters

def test_niederer_points_maps_probe_time_into_L2(tmp_path):
    # real cached header: Label,Points:0,Points:1,Points:2,activationTime
    cfg = tmp_path / "explicit_TNNP_epicardialCells_DT0.001_DX0.1"
    cfg.mkdir()
    (cfg / "explicit_TNNP_epicardialCells_points_DT0001_DX01.csv").write_text(
        "Label,Points:0,Points:1,Points:2,activationTime\n"
        "P1,0.0,0.0,0.007,0.00121026\n"
        "P8,0.02,0.02,0.007,0.0458\n")
    rows = adapters.from_niederer_points(tmp_path)
    p8 = next(r for r in rows if r["field"] == "P8")
    assert p8["case"] == "niederer"
    assert p8["variant"] == "explicit_TNNP_epicardialCells_DT0.001_DX0.1"
    assert p8["dim"] == "3D" and p8["h"] == "0.1"
    assert p8["L2"] == "0.0458" and p8["L1"] == "" and p8["Linf"] == ""
