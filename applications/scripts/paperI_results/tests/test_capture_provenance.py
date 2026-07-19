import json
import os
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "capture_provenance.sh"

def test_emits_valid_json_with_keys(tmp_path):
    case = tmp_path / "case"
    (case / "constant").mkdir(parents=True)
    (case / "system").mkdir()
    (case / "constant" / "electroProperties").write_text("foo\n")
    out = subprocess.run(["bash", str(SCRIPT), str(case)],
                         capture_output=True, text=True, check=True,
                         env={**os.environ, "WM_PROJECT_VERSION": "v2412"})
    doc = json.loads(out.stdout)                      # must be valid JSON
    assert doc["environment"]["openfoam_version"] == "v2412"
    assert doc["inputs"]["electroProperties"].startswith("sha256:")
    assert doc["inputs"]["controlDict"] == "absent"   # not created
    assert set(doc) == {"environment", "inputs", "provenance"}

def test_sha_is_64_hex(tmp_path):
    case = tmp_path / "c"; (case / "system").mkdir(parents=True)
    (case / "system" / "fvSolution").write_text("bar\n")
    out = subprocess.run(["bash", str(SCRIPT), str(case)],
                         capture_output=True, text=True, check=True)
    h = json.loads(out.stdout)["inputs"]["fvSolution"].split(":", 1)[1]
    assert len(h) == 64 and all(c in "0123456789abcdef" for c in h)
