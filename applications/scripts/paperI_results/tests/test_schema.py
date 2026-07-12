import csv

import schema


def test_observed_order_second_order():
    # error / 4 while h / 2  ->  p = 2
    assert round(schema.observed_order(4e-3, 1e-3, 0.1, 0.05), 3) == 2.0


def test_observed_order_guards_bad_input():
    assert schema.observed_order(0, 1e-3, 0.1, 0.05) is None
    assert schema.observed_order(1e-3, 1e-3, 0.1, 0.1) is None
    assert schema.observed_order("x", 1e-3, 0.1, 0.05) is None


def test_fill_rates_blank_on_coarsest_then_order_two():
    rows = [
        {"case": "c", "variant": "", "dim": "1D", "field": "Vm",
         "N": "10", "h": "0.1", "L1": "", "L2": "4e-3", "Linf": "8e-3"},
        {"case": "c", "variant": "", "dim": "1D", "field": "Vm",
         "N": "20", "h": "0.05", "L1": "", "L2": "1e-3", "Linf": "2e-3"},
    ]
    out = schema.fill_rates(rows)
    assert out[0]["rate_L2"] == ""            # coarsest
    assert out[1]["rate_L2"] == "2.00"
    assert out[1]["rate_Linf"] == "2.00"


def test_write_canonical_roundtrip(tmp_path):
    rows = schema.fill_rates([
        {"case": "c", "variant": "", "dim": "1D", "field": "Vm",
         "N": "10", "h": "0.1", "L1": "", "L2": "4e-3", "Linf": "8e-3"},
    ])
    dest = tmp_path / "sub" / "out.csv"
    schema.write_canonical(dest, rows)
    with dest.open() as fh:
        read = list(csv.DictReader(fh))
    assert list(read[0].keys()) == schema.CANONICAL_FIELDS
    assert read[0]["field"] == "Vm"
