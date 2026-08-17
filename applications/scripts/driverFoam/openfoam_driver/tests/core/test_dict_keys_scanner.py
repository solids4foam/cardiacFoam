"""Unit coverage for the C++/catalog drift scanner's path normalization.

_parse_path's scope-token stripping is a syntactic transform (recognize and
strip a "$SCOPE_TOKEN." shape) -- it must not hardcode the one token the
built-in cardiac plugin happens to declare, since a future plugin can
register its own scope token under the same $TOKEN. convention.
"""
from __future__ import annotations

from openfoam_driver.scripts._dict_keys_scanner import _parse_path


def test_strips_the_cardiac_scope_token() -> None:
    path = _parse_path("$ELECTRO_MODEL_COEFFS.myocardiumSolver", is_dynamic=False)
    assert path.normalised == "myocardiumSolver"
    assert path.leaf == "myocardiumSolver"


def test_strips_an_arbitrary_scope_token_of_the_same_shape() -> None:
    path = _parse_path("$SOME_OTHER_PLUGIN_COEFFS.foo.bar", is_dynamic=False)
    assert path.normalised == "foo.bar"
    assert path.leaf == "bar"
    assert path.parents == ("foo",)


def test_leaves_an_unprefixed_path_unchanged() -> None:
    path = _parse_path("myocardiumSolver", is_dynamic=False)
    assert path.normalised == "myocardiumSolver"


def test_leaves_a_dollar_sign_not_matching_the_scope_token_shape_unchanged() -> None:
    path = _parse_path("$notAToken", is_dynamic=False)
    assert path.normalised == "$notAToken"
