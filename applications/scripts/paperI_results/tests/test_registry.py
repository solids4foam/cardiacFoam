import csv
import json
from pathlib import Path
import aggregate

ROOT = Path(__file__).resolve().parents[4]
REG = ROOT / "applications/scripts/driverFoam/verification_experiments.json"

def _rows():
    payload = json.loads(REG.read_text())
    return [{
        "experiment_id": item["experiment_id"],
        "case_dir": item["case_dir"],
        "runner": item["execution"]["runner"],
        "aggregator": item["aggregation"].get("key") or "-",
        "result": item["aggregation"]["result"],
        "reference": item["aggregation"].get("reference") or "-",
    } for item in payload["experiments"]]

def test_contract_schema_and_required_sections():
    payload = json.loads(REG.read_text())
    assert payload["schema_version"] == "1"
    for item in payload["experiments"]:
        assert {"experiment_id", "description", "case_dir", "execution",
                "matrix", "observables", "aggregation"} <= set(item)
        assert item["execution"]["kind"] in {"driver_sweep", "legacy_bash"}

def test_agg_keys_are_registered_in_aggregate():
    known = set(aggregate.CASES)
    assert {r["aggregator"] for r in _rows() if r["aggregator"] != "-"} <= known

def test_expected_case_keys_present():
    keys = {r["experiment_id"] for r in _rows()}
    assert keys == {
        "monodomain_cartesian", "monodomain_temporal",
        "monodomain_tet_generic", "monodomain_tet_frontal",
        "eikonal_cartesian", "eikonal_tet_generic", "eikonal_tet_frontal",
        "eikonal_gradient_tet", "eikonal_bulk_boundary_tet",
        "bidomain_cartesian", "bidomain_temporal", "bidomain_tet_generic",
        "bath_bidomain_cartesian", "bath_bidomain_tet_conformal",
        "purkinje_monodomain_coupled", "niederer_cartesian",
    }


def test_normalized_ids_and_result_names():
    for row in _rows():
        experiment_id = row["experiment_id"]
        assert experiment_id == experiment_id.lower()
        assert row["result"].endswith(f"/{experiment_id}.csv")


def test_registry_paths_are_real():
    for row in _rows():
        case_root = ROOT / row["case_dir"]
        assert (case_root / row["runner"]).is_file()
        assert (case_root / row["result"]).is_file()
        if row["reference"] != "-":
            assert (case_root / row["reference"]).is_file()


def test_normalized_gradient_result_is_complete():
    result = (ROOT / "tutorials/manufacturedSolutions/eikonalECG"
              / "setup/results/eikonal_gradient_tet.csv")
    with result.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 8
    assert {row["mesh_family"] for row in rows} == {"generic", "frontal"}
    assert {row["N"] for row in rows} == {"10", "20", "40", "80"}
    assert all(row["n_cells"] and row["L2_total"] for row in rows)
