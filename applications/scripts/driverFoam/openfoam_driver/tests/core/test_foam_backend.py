import pytest

from openfoam_driver.core.runtime import foam_backend


HEADER = "FoamFile { version 2.0; class dictionary; object d; }\n"


def _dict(tmp_path, body):
    path = tmp_path / "d"
    path.write_text(HEADER + body)
    return path


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("1e-6", 1e-6),
        ("0.0", 0.0),
        ("250", 250),
        ("yes", True),
        ("no", False),
        ("PCG", "PCG"),
        ("constant/purkinjeGraph", "constant/purkinjeGraph"),
    ],
)
def test_coerce_value_maps_strings_to_foamlib_types(raw, expected):
    assert foam_backend.coerce_value(raw) == expected


def test_coerce_value_passes_non_strings_through():
    assert foam_backend.coerce_value(1e-6) == 1e-6
    assert foam_backend.coerce_value(True) is True


def test_update_entry_writes_four_space_separator(tmp_path):
    path = _dict(tmp_path, "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n")
    foam_backend.update_entry(path, "tolerance", "1e-12", scope=["solvers", "Vm"])
    assert "tolerance    1e-12;" in path.read_text()


def test_update_entry_does_not_insert_blank_line(tmp_path):
    path = _dict(tmp_path, "// documented\nendTime 0.3;\n")
    before = path.read_text().count("\n\n")
    foam_backend.update_entry(path, "endTime", "0.4")
    assert path.read_text().count("\n\n") == before


def test_update_entry_fails_closed_on_unknown_key(tmp_path):
    path = _dict(tmp_path, "endTime 0.3;\n")
    with pytest.raises(KeyError):
        foam_backend.update_entry(path, "endTimee", "999")
    assert "endTimee" not in path.read_text()


def test_update_entry_creates_key_when_add_if_missing(tmp_path):
    path = _dict(tmp_path, "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n")
    foam_backend.update_entry(
        path, "purgeWrite", "0", scope=["solvers", "Vm"], add_if_missing=True
    )
    assert "purgeWrite    0;" in path.read_text()


def test_update_entry_add_if_missing_requires_scope(tmp_path):
    """Mirrors tier 1's mutators.py:434 guard so the two tiers agree.

    Tier 1's ``update_foam_entry`` raises this *before* ever falling
    through to this adapter, so a mismatched contract here would be
    untestable through the real call path and misleading in isolation.
    """
    path = _dict(tmp_path, "endTime 0.3;\n")
    with pytest.raises(ValueError, match="add_if_missing requires a scope"):
        foam_backend.update_entry(path, "purgeWrite", "0", add_if_missing=True)


def test_update_entry_handles_brace_inside_quoted_string(tmp_path):
    """The case tier 1 cannot parse: a brace inside a quoted value."""
    path = _dict(
        tmp_path,
        'note  "a value with { an unbalanced brace";\n'
        "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n",
    )
    foam_backend.update_entry(path, "tolerance", "1e-12", scope=["solvers", "Vm"])
    text = path.read_text()
    assert "1e-12" in text
    assert 'note  "a value with { an unbalanced brace";' in text


def test_update_entry_rejects_injected_second_entry(tmp_path):
    path = _dict(tmp_path, "deltaT 1e-06;\n")
    with pytest.raises(ValueError):
        foam_backend.update_entry(path, "deltaT", "1e-6;  rogue  1")
    assert "rogue" not in path.read_text()


def test_update_entry_rejects_directive_value(tmp_path):
    path = _dict(tmp_path, "deltaT 1e-06;\n")
    with pytest.raises(ValueError):
        foam_backend.update_entry(path, "deltaT", '#calc "2*3"')


def test_remove_dict_deletes_block(tmp_path):
    path = _dict(tmp_path, "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n")
    foam_backend.remove_dict(path, "Vm", scope=["solvers"])
    assert "Vm" not in path.read_text()


def test_remove_dict_missing_ok_false_raises(tmp_path):
    path = _dict(tmp_path, "solvers\n{\n}\n")
    with pytest.raises(KeyError):
        foam_backend.remove_dict(path, "Vm", scope=["solvers"], missing_ok=False)


def test_remove_dict_missing_ok_true_is_silent(tmp_path):
    path = _dict(tmp_path, "solvers\n{\n}\n")
    foam_backend.remove_dict(path, "Vm", scope=["solvers"], missing_ok=True)


def test_remove_dict_leaves_no_whitespace_only_line(tmp_path):
    """foamlib's ``del`` leaves the emptied block's line as spaces, not gone.

    Measured directly: deleting ``solvers/Vm`` from
    ``solvers\\n{\\n    Vm {...}\\n}\\n`` leaves ``solvers\\n{\\n    \\n}\\n`` --
    the brace lines survive but the interior line is whitespace, not absent.
    ``"Vm" not in text`` alone does not catch this.
    """
    path = _dict(tmp_path, "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n")
    foam_backend.remove_dict(path, "Vm", scope=["solvers"])
    text = path.read_text()
    assert "Vm" not in text
    assert not any(line.strip() == "" and line != "\n" for line in text.splitlines(keepends=True))


def test_remove_dict_maps_decode_error_to_value_error(tmp_path):
    """``del`` re-parses the whole file, so unrelated garbage elsewhere can

    surface as a ``FoamFileDecodeError`` even when the delete target itself
    is well-formed. Measured directly: deleting a perfectly valid
    ``solvers/Vm`` from a file that *also* contains unrelated malformed text
    raises ``foamlib.FoamFileDecodeError`` (a ``ValueError`` subclass) from
    the ``del`` call, not a ``KeyError`` -- the bare ``except KeyError:`` an
    earlier draft of this function had would let it escape unmapped.
    """
    path = _dict(
        tmp_path,
        "solvers\n{\n    Vm { tolerance 1e-11; }\n}\n"
        "garbage {{{ not valid @@@ ;;; \n",
    )
    with pytest.raises(ValueError):
        foam_backend.remove_dict(path, "Vm", scope=["solvers"])


def test_missing_file_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        foam_backend.update_entry(tmp_path / "nope", "k", "1")
