import aggregate


def test_build_tet_rows_have_rates(tmp_path):
    case_dir = tmp_path / "tutorials/manufacturedSolutions/monodomainTetMMS/setup/results"
    case_dir.mkdir(parents=True)
    (case_dir / "scheme_study.csv").write_text(
        "scheme,N,dx,mono_L2,mono_Linf,ecg_L2,ecg_Linf\n"
        "leastSquares,20,0.030303,0.00298745,0.0213467,0.000184609,0.000203236\n"
        "leastSquares,40,0.0151515,0.000711825,0.00598599,3.81311e-05,3.96357e-05\n"
    )
    rows = aggregate.CASES["tet"](tmp_path)
    fine_vm = [r for r in rows if r["field"] == "Vm" and r["N"] == "40"][0]
    assert fine_vm["rate_L2"] != ""      # rate computed between 20 and 40
