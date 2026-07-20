import checkmesh_parse as cm


def test_parses_nonortho_and_ok():
    log = "Mesh non-orthogonality Max: 42.7 average: 11.3\nMesh OK.\n"
    r = cm.parse_checkmesh_log(log)
    assert r.max_non_orthogonality == 42.7
    assert r.average_non_orthogonality == 11.3
    assert r.mesh_ok is True
    assert r.has_negative_volume is False


def test_detects_negative_volume_and_not_ok():
    log = "Mesh non-orthogonality Max: 5 average: 1\n***Cells with negative volumes found.\n"
    r = cm.parse_checkmesh_log(log)
    assert r.has_negative_volume is True
    assert r.mesh_ok is False
