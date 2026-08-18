import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent.parent))

from openfoam_driver.tests.core.test_cardiac_tutorial_characterization import _current_characterization, _FIXTURE

def update_fixture():
    data = _current_characterization()
    expected = json.loads(_FIXTURE.read_text())
    expected["tutorials"] = data
    _FIXTURE.write_text(json.dumps(expected, indent=2) + "\n")
    print(f"Updated {_FIXTURE}")

if __name__ == "__main__":
    update_fixture()
