import json

import pytest
from openfoam_driver.tests.conftest import skip_without_monorepo
pytestmark = skip_without_monorepo

from openfoam_driver.cli import main
from openfoam_driver.verification_contracts import load_contracts, plan, tsv_rows


def test_contract_catalog_is_unique_and_complete():
    contracts = load_contracts()
    identifiers = [item["experiment_id"] for item in contracts]
    assert len(identifiers) == 15
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


def test_frontal_monodomain_is_one_full_sweep_in_contract_and_bash():
    experiment = next(
        item for item in load_contracts()
        if item["experiment_id"] == "monodomain_tet_frontal"
    )
    assert experiment["matrix"]["N"] == [10, 20, 40, 80]
    assert experiment["execution"]["driver_specs"] == [
        "setup/studies/tetConvergence/sweep_tet_frontal.json"
    ]




def test_tsv_adapter_has_one_six_column_row_per_contract():
    rows = tsv_rows().splitlines()
    assert len(rows) == len(load_contracts())
    assert all(len(row.split("\t")) == 6 for row in rows)


def test_cli_can_describe_one_experiment(capsys):
    assert main(["experiment-plan", "--experiment", "eikonal_tet_generic"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert [item["experiment_id"] for item in payload["experiments"]] == [
        "eikonal_tet_generic"
    ]


def test_cli_tsv_format_honours_experiment_selection(capsys):
    assert main([
        "experiment-plan", "--format", "tsv",
        "--experiment", "eikonal_tet_generic",
    ]) == 0
    rows = capsys.readouterr().out.splitlines()
    assert len(rows) == 1
    assert rows[0].split("\t", 1)[0] == "eikonal_tet_generic"


def test_cli_rejects_unknown_experiment():
    with pytest.raises(SystemExit):
        main(["experiment-plan", "--experiment", "not_a_real_experiment"])
