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
#     ionic_catalog_verification
#
# Description
#     Verifies IONIC_MODEL_CATALOG against the built solver by running
#     listCellModelsVariables and diffing its report. The only durable check:
#     constant naming follows no static rule.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Verify the static ionic catalog against what the solver actually exposes.

Why this cannot be done statically
----------------------------------
An agent writes ``ionicConstantOverrides`` using names from
``IONIC_MODEL_CATALOG``. An unknown name is a **FatalError** at solver startup
(``src/genericWriter/ionicModelIO.C:245-255``), so a stale catalog means a run
that dies -- or a name that silently resolves to something else.

A 2026-08-11 audit established that no static rule describes the naming:

* ``AC_``-prefixed: AlievPanfilov, Courtemanche, Fabbri, Gaur, Grandi,
  PerisYague, Stewart, ToRORd_dynCl, Trovato.
* **Unprefixed**: BuenoOrovio, TNNP, and the three FDA manufactured
  verification models.
* **Mixed within one model**: TWorld carries 273 ``AC_*`` constants and exactly
  one that is not -- ``gnalTissueScale``
  (``src/ionicModels/TWorld/TWorld_2024.H:729``).

The split is by *provenance*: CellML-generated constants get ``AC_``,
hand-added ones do not. So it grows whenever someone hand-adds a constant, and
any rule written down is wrong the next time that happens. The only durable
answer is to ask the built solver.

What makes this authoritative
-----------------------------
``listCellModelsVariables.C:161`` prints ``model.constantVariableNames()``,
which is a thin wrapper (``src/ionicModels/ionicModel/ionicModel.H:518``) over
``ioConstantNames()`` (``:526``) -- the same array handed to the override
matcher at ``src/ionicModels/ionicModel/ionicModel.C:98``. The utility's output
*is* the override vocabulary, not a parallel description of it.

Running it
----------
The live check needs OpenFOAM sourced and the utility built::

    source /path/to/OpenFOAM/etc/bashrc
    (cd applications/utilities/listCellModelsVariables && wmake)
    cd applications/scripts/driverFoam
    uv run pytest openfoam_driver/tests/plugins/cardiacfoam/ \\
        test_ionic_catalog_live_verification.py -v

Without the utility on ``PATH`` every model reports ``skipped`` and
``all_match`` is ``False``. **A skip is never a pass.**

Supersedes ``openfoam_driver/scripts/generate_catalog.py``, which did the same
comparison but took a hand-supplied report file, was referenced by nothing, and
has been broken since the catalog moved to ``plugins/cardiacfoam/`` (it still
imports ``from ionic_model_catalog import ...`` against ``openfoam_driver/``).
Its parser logic is correct and is lifted here.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

_UTILITY = "listCellModelsVariables"

# Sections the report exposes, in the order the utility prints them.
_FIELDS = ("constants", "states", "algebraic")


@dataclass(frozen=True)
class ModelVerificationResult:
    """One model's catalog-vs-runtime comparison.

    ``missing_from_catalog`` is the dangerous direction: the solver exposes a
    name the catalog does not, so an agent never learns the override exists.
    ``extra_in_catalog`` is the loud one: the catalog advertises a name the
    solver lacks, so an agent writes an override that fatals at startup.

    Only non-empty fields appear in either mapping.
    """

    model: str
    status: str  # "match" | "mismatch" | "skipped" | "error"
    missing_from_catalog: dict[str, tuple[str, ...]] = field(default_factory=dict)
    extra_in_catalog: dict[str, tuple[str, ...]] = field(default_factory=dict)
    reason: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "model": self.model,
            "status": self.status,
            "missing_from_catalog": {k: list(v) for k, v in self.missing_from_catalog.items()},
            "extra_in_catalog": {k: list(v) for k, v in self.extra_in_catalog.items()},
            "reason": self.reason,
        }


@dataclass(frozen=True)
class VerificationResult:
    """Whole-run result. ``all_match`` is True only when every model was
    actually compared and agreed -- a skipped model makes it False, because a
    check that did not run is not a check that passed."""

    utility_available: bool
    results: dict[str, ModelVerificationResult]

    @property
    def all_match(self) -> bool:
        return bool(self.results) and all(
            result.status == "match" for result in self.results.values()
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "utility_available": self.utility_available,
            "all_match": self.all_match,
            "results": {name: r.to_json() for name, r in self.results.items()},
        }


def parse_report_text(text: str) -> dict[str, Any]:
    """Parse a listCellModelsVariables report into name lists.

    Lifted from the (now-broken) ``scripts/generate_catalog.py``; the regexes
    match the real output format at ``listCellModelsVariables.C:168-198``. Note
    ``algebraic`` entries carry no ``-->`` value, unlike constants and states.

    Unparseable text yields empty lists rather than raising -- the caller
    surfaces that as a mismatch, so garbage can never read as agreement.
    """
    result: dict[str, Any] = {
        "ionic_model": None,
        "states": [],
        "algebraic": [],
        "constants": [],
    }

    model_match = re.search(r"Selected ionicModel:\s*(\w+)", text)
    if model_match:
        result["ionic_model"] = model_match.group(1)

    for section, value_suffix in (
        ("constants", r"\s*-->"),
        ("states", r"\s*-->"),
        ("algebraic", ""),
    ):
        header = "Ionic " + section
        block = re.search(
            rf"{header}\s*\(\d+\)[^\n]*\n((?:[ \t]+{section}[ \t]*\[\d+\][^\n]*\n)*)",
            text,
        )
        if not block:
            continue
        for line in block.group(1).split("\n"):
            match = re.search(rf"{section}\s*\[\d+\]\s+(\w+){value_suffix}", line)
            if match:
                result[section].append(match.group(1))

    return result


def diff_report(model: str, parsed: dict[str, Any], entry: Any) -> ModelVerificationResult:
    """Compare one parsed report against its catalog entry.

    Order is preserved from the source rather than sorted, so a diff reads in
    the same order as the file it came from.
    """
    missing: dict[str, tuple[str, ...]] = {}
    extra: dict[str, tuple[str, ...]] = {}

    for section in _FIELDS:
        runtime = list(parsed.get(section) or [])
        catalogued = list(getattr(entry, section, ()) or ())
        runtime_set, catalog_set = set(runtime), set(catalogued)

        only_runtime = tuple(n for n in runtime if n not in catalog_set)
        only_catalog = tuple(n for n in catalogued if n not in runtime_set)
        if only_runtime:
            missing[section] = only_runtime
        if only_catalog:
            extra[section] = only_catalog

    status = "match" if not missing and not extra else "mismatch"
    return ModelVerificationResult(
        model=model,
        status=status,
        missing_from_catalog=missing,
        extra_in_catalog=extra,
    )


def find_listCellModelsVariables_binary() -> Path | None:
    """Locate the utility on PATH, or None when OpenFOAM is not sourced."""
    found = shutil.which(_UTILITY)
    return Path(found) if found else None


def _synthesize_case(case_dir: Path, model: str, entry: Any) -> None:
    """Write a minimal single-cell case for one ionic model.

    Uses the driver's own dictionary synthesis rather than hand-written text,
    so this exercises the same path an agent would use to configure a run.
    """
    from openfoam_driver.plugins.cardiacfoam.dict_builder import (
        build_electro_properties,
        build_physics_properties,
    )

    tissues = tuple(getattr(entry, "compatible_tissues", ()) or ("myocyte",))
    (case_dir / "constant").mkdir(parents=True, exist_ok=True)
    (case_dir / "system").mkdir(parents=True, exist_ok=True)

    # The *compactBatched models declare batchedIntegrator required
    # (dict_entries_catalog.py:305-311), so a bare selector set is rejected by
    # the validator. Supply the documented default rather than failing to
    # verify every batched model.
    overrides: dict[str, str] = {}
    if "Batched" in model:
        overrides["$ELECTRO_MODEL_COEFFS.batchedIntegrator"] = "rushLarsen"
    # The FDA manufactured models select their analytical solution by
    # dimensionality and fatal without it (ionicSelector.C:73). Any valid
    # choice exposes the same variable set, so 3D is arbitrary but sufficient.
    # Found by the first live run.
    if "Manufactured" in model:
        # Plain value: dict_builder quotes tokens OpenFOAM cannot lex bare
        # (see _openfoam_value_token). Any valid choice exposes the same
        # variable set, so 3D is arbitrary but sufficient.
        overrides["$ELECTRO_MODEL_COEFFS.dimension"] = "3D"

    (case_dir / "constant" / "electroProperties").write_text(
        build_electro_properties(
            selectors={
                "myocardiumSolver": "singleCellSolver",
                "ionicModel": model,
                "tissue": tissues[0],
            },
            overrides=overrides or None,
        )
    )
    (case_dir / "constant" / "physicsProperties").write_text(
        build_physics_properties(selectors={"type": "electroModel"})
    )

    # Every OpenFOAM application reads system/controlDict via createTime.H
    # before anything else, so the utility fatals without it even though it
    # needs no mesh. Found by the first live run, not by inspection.
    from openfoam_driver.plugins.cardiacfoam.system_templates import build_control_dict

    (case_dir / "system" / "controlDict").write_text(build_control_dict())

    from openfoam_driver.plugins.cardiacfoam.mesh_provisioning import provision_mesh

    # Keyword-only. singleCellSolver copies the bundled 1-cell polyMesh
    # directly, so no blockMesh run is needed before the utility.
    provision_mesh(case_dir=case_dir, myocardium_solver="singleCellSolver")


def _verify_one(model: str, entry: Any, binary: Path, case_root: Path) -> ModelVerificationResult:
    """Run the utility for one model. A per-model failure is reported as that
    model's status, never raised -- one awkward model must not blind the rest."""
    case_dir = case_root / model
    try:
        _synthesize_case(case_dir, model, entry)
    except Exception as exc:  # noqa: BLE001 - reported, not swallowed
        return ModelVerificationResult(
            model=model, status="error", reason=f"could not synthesize a case: {exc}"
        )

    try:
        completed = subprocess.run(
            [str(binary), "-case", str(case_dir)],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return ModelVerificationResult(
            model=model, status="error", reason=f"{_UTILITY} failed to run: {exc}"
        )

    report = completed.stdout
    # The utility also writes postProcessing/listCellModelsVariables.txt; prefer
    # it when stdout was redirected or truncated.
    written = case_dir / "postProcessing" / f"{_UTILITY}.txt"
    if written.is_file():
        report = written.read_text()

    parsed = parse_report_text(report)
    if not any(parsed[section] for section in _FIELDS):
        return ModelVerificationResult(
            model=model,
            status="error",
            reason=(
                f"{_UTILITY} produced no parseable report "
                f"(exit {completed.returncode}); stderr: {completed.stderr[:400]}"
            ),
        )
    return diff_report(model, parsed, entry)


def verify_ionic_catalog(
    model: str | None = None,
    *,
    case_dir: str | Path | None = None,
) -> VerificationResult:
    """Verify one model, or every catalogued model, against the built solver.

    Raises ``KeyError`` for an unknown model name -- that is a caller bug and
    must not be reported as a skip.
    """
    from openfoam_driver.plugins.cardiacfoam.ionic_model_catalog import (
        IONIC_MODEL_CATALOG,
    )

    if model is not None and model not in IONIC_MODEL_CATALOG:
        raise KeyError(
            f"unknown ionic model {model!r}; known models: "
            + ", ".join(sorted(IONIC_MODEL_CATALOG))
        )
    targets = {model: IONIC_MODEL_CATALOG[model]} if model else dict(IONIC_MODEL_CATALOG)

    binary = find_listCellModelsVariables_binary()
    if binary is None:
        return VerificationResult(
            utility_available=False,
            results={
                name: ModelVerificationResult(
                    model=name,
                    status="skipped",
                    reason=(
                        f"{_UTILITY} not on PATH -- source the OpenFOAM bashrc and "
                        "build it. A skip is not a pass."
                    ),
                )
                for name in targets
            },
        )

    if case_dir is not None:
        root = Path(case_dir)
        root.mkdir(parents=True, exist_ok=True)
        results = {
            name: _verify_one(name, entry, binary, root)
            for name, entry in targets.items()
        }
    else:
        with tempfile.TemporaryDirectory(prefix="verify_ionic_catalog_") as tmp:
            root = Path(tmp)
            results = {
                name: _verify_one(name, entry, binary, root)
                for name, entry in targets.items()
            }

    return VerificationResult(utility_available=True, results=results)
