"""Compatibility wrapper for the core generic case factory.

The generic case-folder execution path is now owned by
``openfoam_driver.core.runtime.generic_case`` so the main driver remains
solver-agnostic. This module stays in place only to preserve legacy imports.
"""

from openfoam_driver.core.runtime.generic_case import make_spec as _core_make_spec
from openfoam_driver.plugins.cardiacfoam.generic_case_mutation import apply_case_mutation


def make_spec(**kwargs):
    kwargs.setdefault("_apply_case_mutation", apply_case_mutation)
    return _core_make_spec(**kwargs)


def make_generic_case_spec(**kwargs):
    return make_spec(**kwargs)

__all__ = ["make_spec", "make_generic_case_spec"]
