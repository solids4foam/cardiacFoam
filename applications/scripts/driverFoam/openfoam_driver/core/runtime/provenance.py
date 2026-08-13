#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     provenance
#
# Description
#     The provenance snapshot model: what a run actually consumed, so a
#     resume can be refused when it no longer matches. Pure -- filesystem
#     stat/read and hashing only, no state, no CLI, no plugin imports, no
#     subprocess. Nothing in this module is wired into checkpointing or the
#     resume policy yet; that lands in later phase-2 tasks.
#
#     ``role`` is about consumption, not authorship. ``required_input``
#     (something the workflow reads) versus ``generated_output``. A file can
#     be written by an earlier step and still be a mandatory input --
#     constant/polyMesh is written by blockMesh and then read by the solver.
#     ``origin`` may be recorded for diagnostics but never drives severity.
#     An unclassified file defaults to ``required_input``: a spurious
#     refusal from over-classifying is recoverable, a silent stale replay
#     from under-classifying is the incident this phase exists to prevent.
#
#     A snapshot never overstates its own strength. A file at or under
#     CONTENT_HASH_MAX_BYTES is fingerprinted by content (sha256); a file
#     over the threshold degrades to metadata (size + mtime) rather than
#     paying for a hash of a file that large; a file that cannot be stat'd
#     or read becomes unavailable. Any unavailable component makes the
#     snapshot's ``is_complete`` False -- a check that could not run must
#     never read as a check that passed.
#
#     The 256 MiB threshold is a measured choice, not an arbitrary one:
#     SHA-256 runs at roughly 1.6 GB/s on the hardware this was profiled on,
#     so hashing a file at the threshold costs on the order of 150-200 ms --
#     acceptable inline latency for a provenance check, while multi-GB mesh
#     or field files fall back to the cheap metadata path instead of
#     stalling a resume.
#
#     Identity for both ``compare()`` and the aggregate digest deliberately
#     excludes raw mtime: a checkout or rsync changes mtime without changing
#     content, and reporting that as tampering would make the model
#     untrustworthy. mtime is still recorded on each component for
#     diagnostics, it just never decides whether something changed.
#
#     Diffs name paths, not just kinds -- ``("case_file",)`` tells an agent
#     nothing actionable. Two synthetic diff entries cover the non-file
#     scalars (``workflow_digest`` and ``plugin_identity``) so the aggregate
#     digest can never move while ``compare()`` reports nothing.
#
#     SCHEMA_VERSION encodes the hashing policy, not just the field layout:
#     changing CONTENT_HASH_MAX_BYTES changes what a digest means, so a
#     schema bump is required alongside a threshold change. What a consumer
#     does when two snapshots carry different schema versions is a resume
#     policy question for a later task, not this module's concern.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

SCHEMA_VERSION = "2.0-sha256-256mib"
"""Encodes the hashing policy (algorithm + degrade threshold), not just the
field layout of the model. Bump this whenever CONTENT_HASH_MAX_BYTES or the
hashing algorithm changes -- a snapshot computed under a different policy is
unrecorded, not a differently-valued observation of the same thing."""

CONTENT_HASH_MAX_BYTES = 256 * 1024 * 1024
"""Files at or under this size are fingerprinted by sha256 content hash.
Larger files degrade to a metadata (size + mtime) fingerprint instead of
paying for a hash of a file that large. See the module docstring for the
throughput measurement behind the 256 MiB figure."""


@dataclass(frozen=True, kw_only=True)
class ProvenanceComponent:
    """One fingerprinted filesystem entry (or degraded/unavailable stand-in)
    that a workflow run consumes or produces.

    ``role`` drives severity and defaults to ``"required_input"``. ``origin``
    is diagnostic only. ``strength`` is honest about how much the fingerprint
    actually proves: ``"content"`` (sha256 of the bytes), ``"metadata"``
    (size + mtime only, file too large to hash), or ``"unavailable"`` (could
    not be stat'd or read at all).
    """

    kind: str
    path: str
    role: str = "required_input"
    origin: str | None = None
    method: str = "unavailable"
    strength: str = "unavailable"
    digest: str | None = None
    size: int | None = None
    mtime_ns: int | None = None
    link_target: str | None = None


@dataclass(frozen=True, kw_only=True)
class ProvenanceDiff:
    """One actionable change between two snapshots: what kind of thing,
    which path (or synthetic marker for the non-file scalars), and whether
    it was added, removed, or modified."""

    kind: str
    path: str
    change: str
    role: str | None = None


@dataclass(frozen=True, kw_only=True)
class ProvenanceSnapshot:
    """A complete fingerprint of what a workflow run consumed: its file
    components plus the two non-file scalars (workflow DAG digest and
    plugin identity) that also determine what a resume would actually be
    resuming."""

    schema_version: str
    components: tuple[ProvenanceComponent, ...]
    workflow_digest: str
    plugin_identity: dict[str, Any]
    aggregate_digest: str
    is_complete: bool

    def to_json(self) -> str:
        """Serialize to a JSON string. Inverse of :meth:`from_json`."""
        payload = {
            "schema_version": self.schema_version,
            "components": [asdict(component) for component in self.components],
            "workflow_digest": self.workflow_digest,
            "plugin_identity": self.plugin_identity,
            "aggregate_digest": self.aggregate_digest,
            "is_complete": self.is_complete,
        }
        return json.dumps(payload, sort_keys=True)

    @classmethod
    def from_json(cls, data: str) -> "ProvenanceSnapshot":
        """Reconstruct a :class:`ProvenanceSnapshot` from :meth:`to_json` output."""
        payload = json.loads(data)
        components = tuple(
            ProvenanceComponent(**component) for component in payload["components"]
        )
        return cls(
            schema_version=payload["schema_version"],
            components=components,
            workflow_digest=payload["workflow_digest"],
            plugin_identity=payload["plugin_identity"],
            aggregate_digest=payload["aggregate_digest"],
            is_complete=payload["is_complete"],
        )


def _relative_path_str(path: Path, relative_to: Path) -> str:
    return path.relative_to(relative_to).as_posix()


def component_for_path(
    path: Path,
    *,
    kind: str,
    relative_to: Path,
    role: str = "required_input",
    origin: str | None = None,
) -> ProvenanceComponent:
    """Fingerprint a single filesystem entry, stat'ing (and reading, if small
    enough) it exactly once. Never raises on a missing or unreadable path --
    that degrades to an ``"unavailable"`` component instead.

    A symlink whose resolved target falls outside ``relative_to`` is
    reclassified as ``kind="external_link"`` and fingerprinted through the
    target it points at, with ``link_target`` recording that target for
    diagnostics. A symlink resolving inside the case tree is treated like an
    ordinary file at its own path -- reading it naturally follows the link.
    """
    rel_path = _relative_path_str(path, relative_to)
    resolved_kind = kind
    link_target: str | None = None

    try:
        is_symlink = path.is_symlink()
    except OSError:
        is_symlink = False

    if is_symlink:
        try:
            link_target = str(path.readlink())
        except OSError:
            link_target = None

        try:
            root_resolved = relative_to.resolve()
            target_resolved = path.resolve()
            is_external = not target_resolved.is_relative_to(root_resolved)
        except OSError:
            is_external = True

        if is_external:
            resolved_kind = "external_link"

    try:
        stat_result = path.stat()
    except OSError:
        return ProvenanceComponent(
            kind=resolved_kind,
            path=rel_path,
            role=role,
            origin=origin,
            method="unavailable",
            strength="unavailable",
            digest=None,
            size=None,
            mtime_ns=None,
            link_target=link_target,
        )

    size = stat_result.st_size
    mtime_ns = stat_result.st_mtime_ns

    def _build(method: str, strength: str, digest: str | None) -> ProvenanceComponent:
        """One construction site for all three outcomes.

        The three branches differ only in method/strength/digest. Building them
        separately invites a later field being added to two of three -- and a
        component that silently loses a field is exactly the class of defect
        this module exists to detect.
        """
        return ProvenanceComponent(
            kind=resolved_kind,
            path=rel_path,
            role=role,
            origin=origin,
            method=method,
            strength=strength,
            digest=digest,
            size=size,
            mtime_ns=mtime_ns,
            link_target=link_target,
        )

    if size > CONTENT_HASH_MAX_BYTES:
        return _build("metadata", "metadata", None)

    try:
        content = path.read_bytes()
    except OSError:
        return _build("unavailable", "unavailable", None)

    return _build("sha256", "content", "sha256:" + hashlib.sha256(content).hexdigest())


def _component_digest_payload(component: ProvenanceComponent) -> dict[str, Any]:
    """Fields that participate in identity -- everything except ``mtime_ns``,
    which is diagnostic only (see module docstring: a checkout or rsync
    changes mtime without changing content)."""
    return {
        "kind": component.kind,
        "path": component.path,
        "role": component.role,
        "origin": component.origin,
        "method": component.method,
        "strength": component.strength,
        "digest": component.digest,
        "size": component.size,
        "link_target": component.link_target,
    }


def _component_identity(component: ProvenanceComponent) -> tuple[Any, ...]:
    payload = _component_digest_payload(component)
    return tuple(payload[key] for key in sorted(payload))


def snapshot_from_components(
    components: Sequence[ProvenanceComponent],
    *,
    workflow_digest: str,
    plugin_identity: dict[str, Any],
) -> ProvenanceSnapshot:
    """Assemble a :class:`ProvenanceSnapshot` from already-fingerprinted
    components plus the two non-file scalars that also determine whether a
    resume is resuming the same thing.

    Preserves the order ``components`` was given in -- ``compare()`` walks
    that order so a diff reads in the same order as the file it came from.
    The aggregate digest, in contrast, hashes a canonical (key-sorted,
    component-sorted) representation so it is stable regardless of
    construction order.
    """
    components = tuple(components)
    is_complete = all(component.strength != "unavailable" for component in components)
    aggregate_digest = _compute_aggregate_digest(components, workflow_digest, plugin_identity)
    return ProvenanceSnapshot(
        schema_version=SCHEMA_VERSION,
        components=components,
        workflow_digest=workflow_digest,
        plugin_identity=dict(plugin_identity),
        aggregate_digest=aggregate_digest,
        is_complete=is_complete,
    )


def _reject_nondeterministic(value: Any) -> Any:
    """Refuse to serialise a type whose text form is not guaranteed stable.

    The aggregate digest is only meaningful if it reproduces byte-for-byte for
    an unchanged case. A permissive ``default=str`` would quietly accept an
    object whose ``repr`` embeds a memory address, reintroducing exactly the
    non-determinism the sorting and key-canonicalisation above exist to remove
    -- and it would do so silently, showing up as a spurious stale_inputs
    refusal. Fail loudly instead, so the offending type is fixed at its source.
    """
    raise TypeError(
        f"provenance payload contains a value of type {type(value).__name__!r} "
        "with no deterministic serialisation; add explicit handling rather than "
        "relying on str()"
    )


def _compute_aggregate_digest(
    components: tuple[ProvenanceComponent, ...],
    workflow_digest: str,
    plugin_identity: dict[str, Any],
) -> str:
    sorted_payloads = sorted(
        (_component_digest_payload(component) for component in components),
        key=lambda payload: (payload["kind"], payload["path"]),
    )
    payload = {
        "schema_version": SCHEMA_VERSION,
        "workflow_digest": workflow_digest,
        "plugin_identity": plugin_identity,
        "components": sorted_payloads,
    }
    blob = json.dumps(payload, sort_keys=True, default=_reject_nondeterministic).encode("utf-8")
    return "sha256:" + hashlib.sha256(blob).hexdigest()


def compare(
    before: ProvenanceSnapshot, after: ProvenanceSnapshot
) -> tuple[ProvenanceDiff, ...]:
    """Diff two snapshots into actionable entries: which path was added,
    removed, or modified, plus synthetic ``workflow`` / ``plugin`` entries
    for the non-file scalars so the aggregate digest can never move without
    ``compare()`` reporting something an agent can act on.

    File components are matched by ``(kind, path)``. Identity ignores raw
    mtime (see module docstring). Order follows ``after.components`` for
    added/modified, then ``before.components`` for removed, then the two
    synthetic scalar diffs -- source order throughout, never sorted.
    """
    diffs: list[ProvenanceDiff] = []

    before_index = {(component.kind, component.path): component for component in before.components}
    after_index = {(component.kind, component.path): component for component in after.components}

    for component in after.components:
        key = (component.kind, component.path)
        previous = before_index.get(key)
        if previous is None:
            diffs.append(
                ProvenanceDiff(kind=component.kind, path=component.path, change="added", role=component.role)
            )
        elif _component_identity(previous) != _component_identity(component):
            diffs.append(
                ProvenanceDiff(kind=component.kind, path=component.path, change="modified", role=component.role)
            )

    for component in before.components:
        key = (component.kind, component.path)
        if key not in after_index:
            diffs.append(
                ProvenanceDiff(kind=component.kind, path=component.path, change="removed", role=component.role)
            )

    if before.workflow_digest != after.workflow_digest:
        diffs.append(
            ProvenanceDiff(kind="workflow", path="<workflow_dag>", change="modified", role=None)
        )

    if before.plugin_identity != after.plugin_identity:
        diffs.append(
            ProvenanceDiff(kind="plugin", path="<plugin_identity>", change="modified", role=None)
        )

    return tuple(diffs)
