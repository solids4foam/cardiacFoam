from pathlib import Path

from checkmesh_parse import parse_checkmesh_log

FIXTURE = Path(__file__).parent / "fixtures" / "log.checkMesh.sample"


def test_parses_real_checkmesh_output():
    result = parse_checkmesh_log(FIXTURE.read_text())
    assert result.mesh_ok is True
    assert result.has_negative_volume is False
    assert result.max_non_orthogonality is not None
    assert result.max_non_orthogonality < 1.0  # undistorted mesh from Task 3
