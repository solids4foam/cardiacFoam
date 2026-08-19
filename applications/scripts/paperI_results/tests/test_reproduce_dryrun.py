import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]        # repo root
SCRIPT = ROOT / "reproduce_verification.sh"

def test_dry_run_lists_all_cases():
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run"],
                         capture_output=True, text=True)
    assert out.returncode == 0
    for key in (
        "monodomain_cartesian", "monodomain_tet_generic",
        "eikonal_cartesian", "eikonal_tet_generic",
        "eikonal_gradient_tet", "bidomain_cartesian", "bidomain_tet_generic",
        "bath_bidomain_cartesian", "bath_bidomain_tet_conformal",
        "purkinje_monodomain_coupled", "niederer_cartesian",
    ):
        assert f"== {key} ==" in out.stdout

def test_dry_run_selection_filters():
    out = subprocess.run(["bash", str(SCRIPT), "--dry-run", "monodomain_tet_generic"],
                         capture_output=True, text=True)
    assert "== monodomain_tet_generic ==" in out.stdout
    assert "== purkinje_monodomain_coupled ==" not in out.stdout
