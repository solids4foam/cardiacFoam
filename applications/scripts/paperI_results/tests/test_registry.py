import csv
from pathlib import Path
import aggregate

REG = Path(__file__).resolve().parents[1] / "paperI_cases.tsv"
COLS = ["key", "case_dir", "run_entry", "agg_key", "fresh_rel", "ref_rel"]

def _rows():
    with REG.open(newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if r["key"] and not r["key"].startswith("#")]

def test_header_is_uniform_bash_schema():
    with REG.open(newline="") as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    assert header == COLS            # no run_kind: every case runs via `bash run_entry`

def test_agg_keys_are_registered_in_aggregate():
    known = set(aggregate.CASES)
    assert {r["agg_key"] for r in _rows()} <= known

def test_expected_case_keys_present():
    keys = {r["key"] for r in _rows()}
    assert keys == {"mono_tet", "eikonal_tet", "coupling1D3D_hex", "eikonal_hex", "mono_hex",
                    "bidomain_hex", "bidomain_tet", "bath_hex",
                    "bath_tet", "niederer_hex"}
