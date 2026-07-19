import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]        # repo root
SCRIPT = ROOT / "reproduce_paperI.sh"

def test_dry_run_lists_all_cases():
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run"],
                         capture_output=True, text=True)
    assert out.returncode == 0
    for key in ("tet", "eikonal_tet", "coupling", "eikonal", "mono_spatial",
                "pseudo_ecg_spatial", "bidomain", "bath", "bath_tet", "niederer"):
        assert f"== {key} ==" in out.stdout

def test_dry_run_selection_filters():
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run", "tet"],
                         capture_output=True, text=True)
    assert "== tet ==" in out.stdout and "== coupling ==" not in out.stdout
