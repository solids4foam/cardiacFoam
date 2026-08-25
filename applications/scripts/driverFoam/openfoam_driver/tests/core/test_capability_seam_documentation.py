"""The capability seam contract must stay documented, accurate, and rendered.

A plain "has a docstring" check is satisfied by ``\"\"\"TODO.\"\"\"``, so these
tests assert the four structured fields resolve to things that actually exist:
``:adapts:`` must name a real plugin member, ``:fallback:`` a real
``compatibility.py`` function, ``:consumed-by:`` a module that really touches
the capability. A stale or invented reference fails, not just an empty one.

The same field blocks are the source of the ARCHITECTURE.md seam table, so the
last test keeps the rendered table from drifting away from the code.
"""

from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path

import pytest

from openfoam_driver.core import (
    capability_seams,
    compatibility,
    plugin_capabilities,
    plugin_interface,
)

DRIVER_ROOT = Path(plugin_capabilities.__file__).resolve().parents[2]
GENERATOR = DRIVER_ROOT / "scripts" / "export-capability-seams.py"

REQUIRED_FIELDS = ("adapts", "consumed-by", "fallback", "status")
VALID_STATUSES = {"mandatory", "optional", "mixed"}

CAPABILITY_FIELDS = tuple(plugin_capabilities.PluginCapabilities.__annotations__)


def _protocol_for(field: str):
    annotation = plugin_capabilities.PluginCapabilities.__annotations__[field]
    name = annotation if isinstance(annotation, str) else annotation.__name__
    protocol = getattr(plugin_capabilities, name, None)
    assert protocol is not None, f"{field}: no Protocol named {name!r}"
    return name, protocol


def _fields(field: str) -> dict[str, str]:
    _, protocol = _protocol_for(field)
    return capability_seams.parse_fields(protocol.__doc__)


def _plugin_members() -> set[str]:
    """Every member a plugin may legitimately expose, across both protocols."""
    members: set[str] = set()
    for protocol in (
        plugin_interface.SolverPlugin,
        plugin_interface.SolverPluginOptionalHooks,
    ):
        members |= {name for name in dir(protocol) if not name.startswith("_")}
    return members


def _compatibility_functions() -> set[str]:
    return {name for name in dir(compatibility) if name.startswith("legacy_")}


@pytest.mark.parametrize("field", CAPABILITY_FIELDS)
def test_capability_documents_all_four_fields(field: str) -> None:
    name, protocol = _protocol_for(field)
    doc = protocol.__doc__
    assert doc and doc.strip(), f"{name} has no docstring"

    parsed = _fields(field)
    missing = [key for key in REQUIRED_FIELDS if not parsed.get(key, "").strip()]
    assert not missing, f"{name} is missing or has empty fields: {missing}"

    # Prose, not just a field block: the fields describe the wiring, the prose
    # has to say why the seam exists.
    prose = doc.split(":adapts:")[0].strip()
    assert len(prose) > 80, f"{name} has a field block but no substantive prose"


@pytest.mark.parametrize("field", CAPABILITY_FIELDS)
def test_adapts_names_real_plugin_members(field: str) -> None:
    name, _ = _protocol_for(field)
    declared = _fields(field)["adapts"]
    if declared.strip() == "none":
        return
    members = _plugin_members()
    for member in (item.strip() for item in declared.split(",")):
        assert member in members, (
            f"{name} :adapts: names {member!r}, which is not a member of "
            "SolverPlugin or SolverPluginOptionalHooks"
        )


@pytest.mark.parametrize("field", CAPABILITY_FIELDS)
def test_fallback_names_real_compatibility_functions(field: str) -> None:
    name, _ = _protocol_for(field)
    declared = _fields(field)["fallback"]
    if declared.strip() == "none":
        return
    known = _compatibility_functions()
    for fn in (item.strip() for item in declared.split(",")):
        assert fn in known, (
            f"{name} :fallback: names {fn!r}, which is not a function in "
            "core/compatibility.py"
        )


@pytest.mark.parametrize("field", CAPABILITY_FIELDS)
def test_consumed_by_names_modules_that_touch_the_capability(field: str) -> None:
    name, _ = _protocol_for(field)
    declared = _fields(field)["consumed-by"]
    if declared.startswith("none"):
        return
    # Subset semantics: every listed module must really touch this capability,
    # but the list need not be exhaustive -- adding a consumer must not break
    # the build.
    for relpath in (item.strip() for item in declared.split(",")):
        path = DRIVER_ROOT / relpath
        assert path.is_file(), f"{name} :consumed-by: names missing file {relpath}"
        assert f"capabilities.{field}" in path.read_text(), (
            f"{name} :consumed-by: names {relpath}, which does not reference "
            f"capabilities.{field}"
        )


@pytest.mark.parametrize("field", CAPABILITY_FIELDS)
def test_status_is_a_known_value(field: str) -> None:
    name, _ = _protocol_for(field)
    status = _fields(field)["status"].strip()
    assert status in VALID_STATUSES, f"{name} :status: is {status!r}"


def test_every_ungated_cardiac_fallback_is_accounted_for() -> None:
    """No capability fallback may reach cardiac code without checking plugin_id.

    The import-level boundary test allows compatibility.py to import the
    cardiac plugin; it cannot see whether a given fallback checks *which*
    plugin it is answering for. Six fallbacks did not, so a non-cardiac plugin
    silently inherited cardiac semantics. This test keeps that from returning.
    """
    source = Path(compatibility.__file__).read_text()
    tree = ast.parse(source)
    offenders = []
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef):
            continue
        segment = ast.get_source_segment(source, node) or ""
        if "plugins.cardiacfoam" not in segment:
            continue
        if "org.cardiacfoam" in segment:
            continue
        offenders.append(node.name)

    # legacy_default_driver_context selects the default plugin when none is
    # given (a product decision, not a leak); legacy_generic_case_mutation
    # serves direct callers of core make_spec and is a documented Plan-2 seam.
    allowed = {"legacy_default_driver_context", "legacy_generic_case_mutation"}
    assert set(offenders) <= allowed, (
        f"ungated cardiac fallbacks: {sorted(set(offenders) - allowed)}"
    )


def test_architecture_seam_table_is_up_to_date() -> None:
    result = subprocess.run(
        [sys.executable, str(GENERATOR), "--check"],
        capture_output=True,
        text=True,
        cwd=DRIVER_ROOT,
    )
    assert result.returncode == 0, result.stderr or result.stdout


def test_every_probed_hook_is_declared_somewhere() -> None:
    """A hook an adapter probes must be findable in the public contract.

    Fourteen were not, so a plugin author reading plugin_interface.py could not
    discover the extension points existed -- while not implementing one routed
    them into a compatibility fallback. SolverPluginOptionalHooks closed that;
    this keeps the next hook from reopening it.
    """
    source = Path(plugin_capabilities.__file__).read_text()
    probed = set(re.findall(r'getattr\(\s*self\.plugin,\s*"([a-z_]+)"', source))
    undeclared = sorted(probed - _plugin_members())
    assert not undeclared, (
        "adapters probe hooks that no plugin protocol declares: "
        f"{undeclared}. Add them to SolverPluginOptionalHooks."
    )
