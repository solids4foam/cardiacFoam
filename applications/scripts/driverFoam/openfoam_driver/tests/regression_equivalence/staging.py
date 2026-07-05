"""Agent addressability of regression cases (solver-free).

Two driving paths:

- strict: resolve by the registered entry name (only for mapped cases).
- generic: resolve by the case-folder path, which is not a SPEC_FACTORIES key,
  so ``resolve_entry`` falls through to the generic case-folder branch. This is
  how the agent reasons about cases that have no bespoke spec.
"""
from __future__ import annotations

from typing import Any

from openfoam_driver.core.runtime.registry import resolve_entry
from openfoam_driver.tests.regression_equivalence.registry import RegressionCase


def resolve_generic(case: RegressionCase) -> dict[str, Any]:
    """Resolve a case by its folder path (the generic, non-registered branch).

    Resolving by path skips the registered-name lookup and lands in the
    case-folder branch. Some mapped cases (e.g. niederer) exist under the same
    path as both a registered_tutorial and a case_folder entry, which is
    ambiguous; there we pin the case_folder kind. Cases with a single entry at
    their path resolve without a kind filter.
    """
    try:
        return resolve_entry(case.case_dir)
    except KeyError as exc:
        if "ambiguous" not in str(exc):
            raise
        return resolve_entry(case.case_dir, entry_kind="case_folder")


def resolve_strict(case: RegressionCase) -> dict[str, Any]:
    if not case.mapped:
        raise ValueError(f"{case.case_dir} has no registered entry")
    return resolve_entry(case.entry_name)
