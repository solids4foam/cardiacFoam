import csv
import check_against_reference as chk

_HEADER = ["case","variant","dim","N","h","field","L1","L2","Linf","rate_L2","rate_Linf"]

def _write(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_HEADER)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in _HEADER})

def _row(**kw):
    base = dict(case="tet", variant="leastSquares", dim="3D", N="20", h="0.03",
                field="Vm", L1="", L2="1e-3", Linf="2e-3", rate_L2="2.00", rate_Linf="1.80")
    base.update(kw)
    return base

def test_identical_passes(tmp_path):
    f, r = tmp_path / "f.csv", tmp_path / "r.csv"
    _write(f, [_row()]); _write(r, [_row()])
    ok, failures = chk.compare(f, r)
    assert ok and failures == []

def test_within_rel_tolerance_passes(tmp_path):
    f, r = tmp_path / "f.csv", tmp_path / "r.csv"
    _write(f, [_row(L2="1.02e-3")]); _write(r, [_row(L2="1e-3")])  # 2% < 5%
    ok, _ = chk.compare(f, r)
    assert ok

def test_error_out_of_tolerance_fails(tmp_path):
    f, r = tmp_path / "f.csv", tmp_path / "r.csv"
    _write(f, [_row(L2="1.5e-3")]); _write(r, [_row(L2="1e-3")])  # 50% > 5%
    ok, failures = chk.compare(f, r)
    assert not ok and any("L2" in m for m in failures)

def test_rate_out_of_abs_tolerance_fails(tmp_path):
    f, r = tmp_path / "f.csv", tmp_path / "r.csv"
    _write(f, [_row(rate_L2="1.50")]); _write(r, [_row(rate_L2="2.00")])  # 0.5 > 0.15
    ok, failures = chk.compare(f, r)
    assert not ok and any("rate_L2" in m for m in failures)

def test_missing_reference_row_fails(tmp_path):
    f, r = tmp_path / "f.csv", tmp_path / "r.csv"
    _write(f, [_row()]); _write(r, [_row(), _row(field="u1")])
    ok, failures = chk.compare(f, r)
    assert not ok and any("missing row" in m for m in failures)

def test_cli_skip_on_missing_reference(tmp_path, capsys):
    f = tmp_path / "f.csv"; _write(f, [_row()])
    code = chk.main([str(f), str(tmp_path / "nope.csv")])
    assert code == 2
    assert "SKIP" in capsys.readouterr().out
