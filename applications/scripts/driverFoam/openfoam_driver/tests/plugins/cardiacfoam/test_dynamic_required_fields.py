"""Required-field enforcement for every dynamic block, not just one.

``specs/validation.py``'s generic required-field pass skips ``dynamic_path``
entries outright ("concrete required leaves are the user's responsibility
when those blocks are actually configured"), and
``_evaluate_dynamic_required_fields`` was written to discharge that
responsibility -- but only for ``conductionNetworkDomains.<name>.*``.

Ten dynamic catalogue entries carry a required rule; that hardcoding left
six of them unenforced. Three are inert anyway (the ``ecgDomains.<name>.
sampling.*`` trio carry ``typical_value``s, so the builder satisfies them
before validation ever runs). The other three are real: every leaf under
``domainCouplings.<name>`` has an empty ``typical_value``, so a configured
coupling block genuinely can go to the solver missing a required key.

Not in scope: ``_evaluate_pvj_resistance_requirement``. That is a
graph-file-aware check the catalogue cannot express -- "required unless a
FILE provides pvjResistances" -- and it stays hand-written.
"""

from __future__ import annotations

from openfoam_driver.plugins.cardiacfoam.validation import (
    _evaluate_dynamic_required_fields,
)


def _fields(errors) -> list[str]:
    return sorted(e.field for e in errors)


def test_configured_domain_coupling_missing_required_leaf_is_reported():
    """A declared domainCouplings block must have its required leaves."""
    context = {
        "myocardiumSolver": "monodomainSolver",
        "domainCouplings.couplingA.electroDomainCoupler": "eikonalMonodomainPvjCoupler",
        # conductionNetworkDomain and couplingMode are required and absent
    }
    fields = _fields(_evaluate_dynamic_required_fields(context))
    assert "domainCouplings.couplingA.conductionNetworkDomain" in fields
    assert "domainCouplings.couplingA.couplingMode" in fields


def test_a_complete_domain_coupling_is_silent():
    context = {
        "myocardiumSolver": "monodomainSolver",
        "domainCouplings.couplingA.electroDomainCoupler": "monodomainPvjCoupler",
        "domainCouplings.couplingA.conductionNetworkDomain": "purkinjeNetwork",
        "domainCouplings.couplingA.couplingMode": "twoWay",
    }
    assert _evaluate_dynamic_required_fields(context) == []


def test_an_undeclared_block_is_not_invented():
    """No block configured means nothing to require."""
    assert _evaluate_dynamic_required_fields({"myocardiumSolver": "monodomainSolver"}) == []


def test_sibling_instances_are_scoped_independently():
    """One coupling's value must never satisfy another's requirement."""
    context = {
        "myocardiumSolver": "monodomainSolver",
        "domainCouplings.couplingA.electroDomainCoupler": "monodomainPvjCoupler",
        "domainCouplings.couplingA.conductionNetworkDomain": "netA",
        "domainCouplings.couplingA.couplingMode": "twoWay",
        "domainCouplings.couplingB.electroDomainCoupler": "monodomainPvjCoupler",
        # couplingB is missing both of its own required leaves
    }
    fields = _fields(_evaluate_dynamic_required_fields(context))
    assert all(f.startswith("domainCouplings.couplingB.") for f in fields), fields
    assert "domainCouplings.couplingB.couplingMode" in fields


def test_conduction_network_enforcement_is_preserved():
    """The block that already worked must keep working."""
    context = {
        "myocardiumSolver": "monodomainSolver",
        "conductionNetworkDomains.netA.conductionSystemDomain": "purkinjeGraphModel",
    }
    fields = _fields(_evaluate_dynamic_required_fields(context))
    assert all(f.startswith("conductionNetworkDomains.netA.") for f in fields), fields
