"""Contract: the threat-model doc exists and AGENT_GUIDE points to it."""
from pathlib import Path


def test_security_doc_exists_and_is_linked_from_agent_guide() -> None:
    driver_root = Path(__file__).resolve().parents[3]  # .../driverFoam
    security = driver_root / "SECURITY.md"
    guide = driver_root / "AGENT_GUIDE.md"
    assert security.is_file(), "SECURITY.md must exist"
    assert "SECURITY.md" in guide.read_text(), "AGENT_GUIDE must link SECURITY.md"
