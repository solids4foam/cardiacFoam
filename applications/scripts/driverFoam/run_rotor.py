import sys
from pathlib import Path
from openfoam_driver.core.services.runner import run_workflow
from openfoam_driver.core.runtime.models import CaseRunSpec
from openfoam_driver.core.services.registry_service import resolve_case_folder

# Resolve generic case
spec = resolve_case_folder(Path("../../tutorials/electrophysiologyProtocols/rotorInstability"))

# Override
spec.config_overrides["system/controlDict"] = {
    "endTime": 0.04,
    "writeInterval": 0.04
}
spec.config_overrides["constant/electroProperties"] = {
    "ionic": {
        "export": "(Vm FakeJsi)"
    }
}

result = run_workflow(spec, output_dir=Path("/tmp/rotor_test"))
print(result)
