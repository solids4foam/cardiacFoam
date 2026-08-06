import csv, keyset_gate as kg

def _write(p, rows):
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case","variant","dim","N","field","L2"])
        w.writeheader()
        for r in rows: w.writerow(r)

def test_equal_keysets_pass(tmp_path):
    a, b = tmp_path/"a.csv", tmp_path/"b.csv"
    row = dict(case="c", variant="structured", dim="3D", N="10", field="Vm", L2="1e-3")
    _write(a, [row]); _write(b, [dict(row, L2="9e9")])  # value differs, key same
    assert kg.main([str(a), str(b)]) == 0

def test_extra_fresh_field_row_fails(tmp_path):
    a, b = tmp_path/"a.csv", tmp_path/"b.csv"
    base = dict(case="c", variant="structured", dim="3D", N="10", field="Vm", L2="1e-3")
    extra = dict(base, field="phiE")
    _write(a, [base, extra]); _write(b, [base])
    assert kg.main([str(a), str(b)]) == 1

def test_extra_fresh_N_row_fails(tmp_path):  # N must be in the key
    a, b = tmp_path/"a.csv", tmp_path/"b.csv"
    base = dict(case="c", variant="structured", dim="3D", N="10", field="Vm", L2="1e-3")
    extra = dict(base, N="20")
    _write(a, [base, extra]); _write(b, [base])
    assert kg.main([str(a), str(b)]) == 1
