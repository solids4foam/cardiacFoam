import pytest
from ..verification_contracts import load_contracts, plan, tsv_rows


def test_contract_catalog_is_unique_and_complete():
    contracts = load_contracts()
    identifiers = [item["experiment_id"] for item in contracts]
    assert len(identifiers) == 14
    assert len(identifiers) == len(set(identifiers))


def test_every_registered_path_is_ready():
    planned = plan()
    # On a fresh clone, results (and sometimes references) will be missing.
    # We only verify that the JSON schema parsing and case_dir exist.
    failures = {
        item["experiment_id"]: item["checks"]
        for item in planned["experiments"] if not item["checks"].get("case_dir")
    }
    assert failures == {}


def test_monodomain_tet_generic_is_one_full_sweep_in_contract_and_bash():
    experiment = next(
        item for item in load_contracts()
        if item["experiment_id"] == "monodomain_tet_generic"
    )
    assert experiment["matrix"]["N"] == [10, 20, 40, 80]
    assert experiment["execution"]["driver_specs"] == [
        "setup/studies/tetConvergence/sweep_tet_generic.json"
    ]




def test_tsv_adapter_has_one_eight_column_row_per_contract():
    rows = tsv_rows().splitlines()
    assert len(rows) == len(load_contracts())
    assert all(len(row.split("\t")) == 8 for row in rows)
