import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]        # repo root
SCRIPT = ROOT / "reproduce_paperI.sh"

def test_dry_run_lists_all_cases():
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run"],
                         capture_output=True, text=True)
    assert out.returncode == 0
    for key in ("mono_tet", "eikonal_tet", "coupling1D3D_hex", "eikonal_hex", "mono_hex",
                "bidomain_hex", "bath_hex", "bath_tet", "niederer_hex"):
        assert f"== {key} ==" in out.stdout

def test_dry_run_selection_filters():
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run", "mono_tet"],
                         capture_output=True, text=True)
    assert "== mono_tet ==" in out.stdout and "== coupling1D3D_hex ==" not in out.stdout
