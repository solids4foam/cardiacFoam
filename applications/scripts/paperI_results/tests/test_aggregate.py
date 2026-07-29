import aggregate


def test_build_tet_rows_have_rates(tmp_path):
    case_dir = tmp_path / "tutorials/manufacturedSolutions/monodomainPseudoECG/setup/results"
    case_dir.mkdir(parents=True)
    (case_dir / "scheme_study.csv").write_text(
        "scheme,N,dx,mono_L2,mono_Linf,ecg_L2,ecg_Linf\n"
        "leastSquares,20,0.030303,0.00298745,0.0213467,0.000184609,0.000203236\n"
        "leastSquares,40,0.0151515,0.000711825,0.00598599,3.81311e-05,3.96357e-05\n"
    )
    rows = aggregate.CASES["mono_tet"](tmp_path)
    fine_vm = [r for r in rows if r["field"] == "Vm" and r["N"] == "40"][0]
    assert fine_vm["rate_L2"] != ""      # rate computed between 20 and 40


def test_build_bidomain_tet_rows_have_rates(tmp_path):
    case_dir = tmp_path / "tutorials/manufacturedSolutions/bidomain/setup/results"
    case_dir.mkdir(parents=True)
    (case_dir / "scheme_study.csv").write_text(
        "scheme,N,dx,vm_L2,vm_Linf,phiE_L2,phiE_Linf\n"
        "leastSquares,20,0.030303,0.00298745,0.0213467,0.00184609,0.00203236\n"
        "leastSquares,40,0.0151515,0.000711825,0.00598599,0.000381311,0.000396357\n"
    )
    rows = aggregate.CASES["bidomain_tet"](tmp_path)
    fine_vm = [r for r in rows if r["field"] == "Vm" and r["N"] == "40"][0]
    assert fine_vm["rate_L2"] != ""      # rate computed between 20 and 40


def test_build_mono_spatial_rows_have_rates(tmp_path):
    archive = (
        tmp_path / "tutorials/manufacturedSolutions/monodomainPseudoECG"
        / "driverPostProcessingArchive_postProcessing"
    )
    archive.mkdir(parents=True)
    (archive / "3D_10_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     4e-3   4e-3   8e-3\n"
        "\nGrid spacing (dx)     = 0.1\n"
    )
    (archive / "3D_20_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     1e-3   1e-3   2e-3\n"
        "\nGrid spacing (dx)     = 0.05\n"
    )
    rows = aggregate.CASES["mono_hex"](tmp_path)
    fine = [r for r in rows if r["N"] == "20"][0]
    assert fine["rate_L2"] == "2.00"


def test_build_bidomain_rows_have_rates(tmp_path):
    archive = tmp_path / "tutorials/manufacturedSolutions/bidomain/driverPostProcessingArchive_postProcessing"
    archive.mkdir(parents=True)
    (archive / "3D_10_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     4e-3   4e-3   8e-3\n"
        "\nGrid spacing (dx)     = 0.1\n"
    )
    (archive / "3D_20_cells_implicit.dat").write_text(
        "Field     L1-error       L2-error       Linf-error\n"
        "Vm     1e-3   1e-3   2e-3\n"
        "\nGrid spacing (dx)     = 0.05\n"
    )
    rows = aggregate.CASES["bidomain_hex"](tmp_path)
    fine = [r for r in rows if r["N"] == "20"][0]
    assert fine["rate_L2"] == "2.00"


