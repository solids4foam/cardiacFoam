"""Round-trip stability of the agent's dict layer (solver-free).

The agent's dict builders synthesize normalized dicts rather than copying the
hand-authored tutorial files, so byte/key-set equality against the committed
dicts is not a valid bar (see memory
``project_driverfoam_dict_synthesizer_not_copier``). What *is* a valid, cheap
invariant: ``build ∘ parse`` must be idempotent on its own generated text.

    build(parse(committed))  ==  build(parse(build(parse(committed))))

i.e. once the agent has synthesized a dict, re-parsing and re-building it must
reproduce the same bytes. If it does not, the agent's generate/parse pair is
internally inconsistent — a real defect independent of any tutorial's
hand-formatting.

Note: comparing the *parsed* selectors/overrides directly is NOT a valid test —
``parse_electro_properties`` deliberately omits values equal to a catalog
``typical_value``, so a default-valued override (e.g. ``groundElectrode yes``)
appears only on the first parse. The generated *text* is the stable object.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

from openfoam_driver.specs.common import tutorials_root_default
from openfoam_driver.specs.dict_builder import (
    build_electro_properties,
    parse_electro_properties,
)
from openfoam_driver.tests.regression_equivalence.registry import RegressionCase


def _build_from_parse(path: Path) -> str:
    parsed = parse_electro_properties(path)
    return build_electro_properties(parsed["selectors"], overrides=parsed["overrides"])


def _write_temp(text: str) -> Path:
    with tempfile.NamedTemporaryFile(
        "w", suffix=".electroProperties", delete=False
    ) as fh:
        fh.write(text)
        return Path(fh.name)


def electro_build_parse_fixpoint(case: RegressionCase) -> tuple[str, str]:
    """Return (once, twice) generated electroProperties texts.

    ``once == twice`` proves ``build ∘ parse`` is idempotent for this case.
    """
    committed = tutorials_root_default() / case.case_dir / "constant/electroProperties"
    once = _build_from_parse(committed)
    tmp = _write_temp(once)
    try:
        twice = _build_from_parse(tmp)
    finally:
        tmp.unlink(missing_ok=True)
    return once, twice
