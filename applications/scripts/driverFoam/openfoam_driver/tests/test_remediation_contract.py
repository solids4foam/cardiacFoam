from openfoam_driver.core.runtime.remediation import (
    STATIC_REMEDIATION_HINTS,
    RemediationHint,
)
from openfoam_driver.dict_entries import (
    CONTROL_DICT_ENTRIES,
    ELECTRO_PROPERTY_ENTRY_GROUPS,
    PHYSICS_PROPERTY_ENTRIES,
)

_PREFIX = "$ELECTRO_MODEL_COEFFS."


def _addressable_leaves() -> set[str]:
    leaves: set[str] = set()
    for e in CONTROL_DICT_ENTRIES:
        leaves.add(e.driver_path)                       # e.g. "deltaT"
    for e in PHYSICS_PROPERTY_ENTRIES:
        leaves.add(e.driver_path)
    for group in ELECTRO_PROPERTY_ENTRY_GROUPS.values():
        for e in group:
            dp = e.driver_path
            leaves.add(dp[len(_PREFIX):] if dp.startswith(_PREFIX) else dp)
    return leaves


def _all_hints() -> list[RemediationHint]:
    hints: list[RemediationHint] = []
    for group in STATIC_REMEDIATION_HINTS.values():
        hints.extend(group)
    return hints


def test_every_hint_driver_path_is_catalog_addressable():
    leaves = _addressable_leaves()
    for h in _all_hints():
        if not h.driver_path:        # advisory hints are exempt
            continue
        assert h.driver_path in leaves, (
            f"hint for {h.diagnostic_code or h.source!r} targets non-catalog "
            f"driver_path {h.driver_path!r}"
        )
