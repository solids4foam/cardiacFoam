import adapters

def test_eikonal_tet_scheme_study_two_fields(tmp_path):
    p = tmp_path / "scheme_study.csv"
    p.write_text(
        "scheme,N,dx,activationTime_L2,activationTime_Linf,ecg_L2,ecg_Linf\n"
        "leastSquares,10,0.0588235,0.068,0.30,0.071,0.10\n"
        "leastSquares,20,0.030303,0.017,0.15,0.017,0.03\n"
    )
    rows = adapters.from_eikonal_tet_scheme_study(p)
    fields = {r["field"] for r in rows}
    assert fields == {"activationTime", "Phi_e"}
    at10 = next(r for r in rows if r["field"] == "activationTime" and r["N"] == "10")
    assert at10 == {"case": "eikonal_tet", "variant": "leastSquares", "dim": "3D",
                    "N": "10", "h": "0.0588235", "field": "activationTime",
                    "L1": "", "L2": "0.068", "Linf": "0.30"}
