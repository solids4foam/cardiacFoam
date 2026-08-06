import sys
from openfoam_driver.specs.validation import validate_run, ValidationError
from openfoam_driver.dict_entries import DictEntry
from openfoam_driver.core.runtime.run_model import RunDocument

entry = DictEntry(
    driver_path="test",
    description="Test",
    phases=frozenset({"physics"}),
    forbidden_when={"myocardiumSolver": "eikonalSolver"}
)

run = RunDocument(
    id="1", name="1", status="draft",
    config={"physics": {"test": "val", "myocardiumSolver": "eikonalSolver"}}
)

errors = validate_run(run, entries=[entry])
print(f"Errors: {errors}")
