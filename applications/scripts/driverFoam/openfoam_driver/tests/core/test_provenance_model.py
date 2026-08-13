"""The snapshot model: structured, consumption-aware, honest about strength."""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.runtime.provenance import (
    CONTENT_HASH_MAX_BYTES,
    SCHEMA_VERSION,
    ProvenanceSnapshot,
    compare,
    component_for_path,
    snapshot_from_components,
)


def _component(path: Path, root: Path, **kw):
    return component_for_path(path, kind="case_file", relative_to=root, **kw)


def _snapshot(tmp_path: Path, text: str) -> ProvenanceSnapshot:
    target = tmp_path / "controlDict"
    target.write_text(text)
    return snapshot_from_components(
        (_component(target, tmp_path),),
        workflow_digest="sha256:wf",
        plugin_identity={"id": "org.example", "version": "1"},
    )


def test_a_file_defaults_to_required_input(tmp_path: Path) -> None:
    """Unknown means required -- a spurious refusal is recoverable, a silent
    stale replay is the incident this phase exists to prevent."""
    target = tmp_path / "controlDict"
    target.write_text("deltaT 0.001;\n")
    component = _component(target, tmp_path)
    assert component.role == "required_input"
    assert component.method == "sha256"
    assert component.strength == "content"


def test_a_generated_file_the_solver_reads_is_still_a_required_input(
    tmp_path: Path,
) -> None:
    """constant/polyMesh is written by blockMesh and read by the solver;
    origin must not decide severity."""
    target = tmp_path / "owner"
    target.write_text("mesh")
    component = _component(target, tmp_path, role="required_input", origin="generated")
    assert component.role == "required_input"
    assert component.origin == "generated"


def test_a_snapshot_is_stable_across_repeated_computation(tmp_path: Path) -> None:
    assert _snapshot(tmp_path, "a").aggregate_digest == _snapshot(tmp_path, "a").aggregate_digest
    assert compare(_snapshot(tmp_path, "a"), _snapshot(tmp_path, "a")) == ()


def test_a_diff_names_the_path_and_the_change(tmp_path: Path) -> None:
    before = _snapshot(tmp_path, "deltaT 0.001;\n")
    after = _snapshot(tmp_path, "deltaT 0.002;\n")
    diffs = compare(before, after)
    assert len(diffs) == 1
    assert diffs[0].path == "controlDict"
    assert diffs[0].change == "modified"
    assert diffs[0].role == "required_input"


def test_added_and_removed_components_are_distinguished(tmp_path: Path) -> None:
    a = tmp_path / "a"; a.write_text("1")
    b = tmp_path / "b"; b.write_text("2")
    one = snapshot_from_components(
        (_component(a, tmp_path),), workflow_digest="w", plugin_identity={},
    )
    two = snapshot_from_components(
        (_component(a, tmp_path), _component(b, tmp_path)),
        workflow_digest="w", plugin_identity={},
    )
    assert [d.change for d in compare(one, two)] == ["added"]
    assert [d.change for d in compare(two, one)] == ["removed"]


def test_a_changed_workflow_produces_an_actionable_diff(tmp_path: Path) -> None:
    """Without a synthetic entry the aggregate changes while compare()
    returns nothing an agent can act on."""
    target = tmp_path / "f"; target.write_text("x")
    component = _component(target, tmp_path)
    a = snapshot_from_components(
        (component,), workflow_digest="sha256:a", plugin_identity={"id": "p"},
    )
    b = snapshot_from_components(
        (component,), workflow_digest="sha256:b", plugin_identity={"id": "p"},
    )
    diffs = compare(a, b)
    assert [(d.kind, d.path, d.change) for d in diffs] == [
        ("workflow", "<workflow_dag>", "modified")
    ]


def test_a_changed_plugin_identity_produces_an_actionable_diff(tmp_path: Path) -> None:
    target = tmp_path / "f"; target.write_text("x")
    component = _component(target, tmp_path)
    a = snapshot_from_components(
        (component,), workflow_digest="w", plugin_identity={"id": "p", "version": "1"},
    )
    b = snapshot_from_components(
        (component,), workflow_digest="w", plugin_identity={"id": "p", "version": "2"},
    )
    diffs = compare(a, b)
    assert [(d.kind, d.path, d.change) for d in diffs] == [
        ("plugin", "<plugin_identity>", "modified")
    ]


def test_a_file_over_the_threshold_degrades_and_says_so(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(
        "openfoam_driver.core.runtime.provenance.CONTENT_HASH_MAX_BYTES", 4
    )
    big = tmp_path / "mesh"; big.write_bytes(b"0123456789")
    component = _component(big, tmp_path)
    assert component.method == "metadata"
    assert component.strength == "metadata"
    assert component.digest is None
    assert component.size == 10 and component.mtime_ns is not None


def test_a_content_hash_ignores_a_timestamp_only_touch(tmp_path: Path) -> None:
    """A checkout or rsync changes mtime without changing content. A
    content-based snapshot must not report that as tampering."""
    import os

    target = tmp_path / "controlDict"
    target.write_text("deltaT 0.001;\n")
    before = _component(target, tmp_path)
    os.utime(target, (1_800_000_000, 1_800_000_000))
    after = _component(target, tmp_path)
    assert before.digest == after.digest


def test_an_unavailable_component_makes_the_snapshot_partial(tmp_path: Path) -> None:
    missing = _component(tmp_path / "gone", tmp_path)
    assert missing.method == "unavailable"
    assert missing.strength == "unavailable"
    assert snapshot_from_components(
        (missing,), workflow_digest="w", plugin_identity={},
    ).is_complete is False
    present = tmp_path / "here"; present.write_text("x")
    assert snapshot_from_components(
        (_component(present, tmp_path),), workflow_digest="w", plugin_identity={},
    ).is_complete is True


def test_an_external_symlink_records_the_target_it_fingerprinted(tmp_path: Path) -> None:
    outside = tmp_path / "outside"; outside.mkdir()
    external = outside / "shared_mesh"; external.write_text("shared")
    case = tmp_path / "case"; case.mkdir()
    link = case / "mesh"; link.symlink_to(external)
    component = component_for_path(link, kind="case_file", relative_to=case)
    assert component.kind == "external_link"
    assert component.path == "mesh"
    assert component.link_target == str(external)
    assert component.digest is not None


def test_the_schema_version_encodes_the_hashing_policy() -> None:
    assert SCHEMA_VERSION.startswith("2.")
    assert "sha256" in SCHEMA_VERSION
    assert CONTENT_HASH_MAX_BYTES == 256 * 1024 * 1024


def test_a_snapshot_round_trips_through_json(tmp_path: Path) -> None:
    original = _snapshot(tmp_path, "deltaT 0.001;\n")
    assert ProvenanceSnapshot.from_json(original.to_json()) == original


def test_a_nondeterministic_payload_type_is_refused_not_stringified():
    """The digest is only meaningful if it reproduces byte-for-byte.

    A permissive ``default=str`` would silently accept an object whose repr
    embeds a memory address, reintroducing the non-determinism that sorting and
    key-canonicalisation exist to remove -- and it would surface as a spurious
    stale_inputs refusal rather than as an error. Fail at the source instead.
    """
    import pytest

    class Opaque:
        pass

    with pytest.raises(TypeError, match="deterministic"):
        snapshot_from_components(
            (),
            workflow_digest="sha256:wf",
            plugin_identity={"id": Opaque()},
        )


def test_the_three_component_outcomes_share_one_construction_site(tmp_path, monkeypatch):
    """All three branches must carry every field, so a later addition cannot
    reach two of three."""
    monkeypatch.setattr(
        "openfoam_driver.core.runtime.provenance.CONTENT_HASH_MAX_BYTES", 4
    )
    small = tmp_path / "small"; small.write_text("x")
    big = tmp_path / "big"; big.write_bytes(b"0123456789")
    missing = tmp_path / "gone"

    built = [
        component_for_path(p, kind="case_file", relative_to=tmp_path)
        for p in (small, big, missing)
    ]
    assert [c.method for c in built] == ["sha256", "metadata", "unavailable"]
    assert [c.strength for c in built] == ["content", "metadata", "unavailable"]
    # Every branch populates the identity fields, not just the varying three.
    for component in built:
        assert component.kind == "case_file"
        assert component.role == "required_input"
        assert component.path
