import adapters


def test_tet_scheme_study(tmp_path):
    p = tmp_path / "scheme_study.csv"
    p.write_text(
        "scheme,N,dx,mono_L2,mono_Linf,ecg_L2,ecg_Linf\n"
        "GaussLinear,10,0.0588235,0.0224588,0.152636,0.0038452,0.00426517\n"
        "leastSquares,80,0.00757576,0.000170952,0.00191104,6.95364e-06,8.29839e-06\n"
    )
    rows = adapters.from_tet_scheme_study(p)
    gl_vm = [r for r in rows if r["variant"] == "GaussLinear" and r["field"] == "Vm"][0]
    assert gl_vm["L2"] == "0.0224588" and gl_vm["dim"] == "3D" and gl_vm["h"] == "0.0588235"
    assert any(r["field"] == "Phi_e" and r["variant"] == "leastSquares" for r in rows)


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


def test_eikonal_activation_folds_2d(tmp_path):
    summary = tmp_path / "act.csv"
    summary.write_text(
        "Dimension,N,activation_L1,activation_L2,activation_Linf,source\n"
        "1D,10,0.00158772,0.00215763,0.00486533,x\n"
        "3D,10,0.00483561,0.00752141,0.0310654,x\n"
    )
    dat2d = tmp_path / "2D_10_cells_eikonal_manufacturedEikonalActivationTime.dat"
    dat2d.write_text("activationTime 0.00288871 0.00415336 0.0121584\n")
    rows = adapters.from_eikonal_activation(summary, extra_2d_dats=[(10, dat2d)])
    r2d = [r for r in rows if r["dim"] == "2D"][0]
    assert r2d["L2"] == "0.00415336" and r2d["field"] == "psi" and r2d["h"] == "0.1"
    r1d = [r for r in rows if r["dim"] == "1D"][0]
    assert r1d["h"] == "0.1"


def test_eikonal_ecg_max_and_mean(tmp_path):
    p = tmp_path / "ecg.csv"
    p.write_text(
        "Dimension,N,max_L1_err_ref,mean_L1_err_ref,max_L2_err_ref,mean_L2_err_ref,"
        "max_Linf_err_ref,mean_Linf_err_ref\n"
        "3D,80,1e-6,1e-6,1e-6,1e-6,1.11883e-05,9.28982e-06\n"
    )
    rows = adapters.from_eikonal_ecg(p)
    assert any(r["field"] == "Phi_e_max" and r["Linf"] == "1.11883e-05" for r in rows)
    assert any(r["field"] == "Phi_e_mean" and r["Linf"] == "9.28982e-06" for r in rows)
