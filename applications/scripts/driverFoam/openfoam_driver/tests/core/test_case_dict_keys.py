"""Warn-level sweep of a case dictionary's keys against the plugin catalogue.

This is the direction the C++ scanner cannot see. ``unmatched_cxx_reads``
asks "what keys does the solver accept?", which lives only in C++ source.
This asks "what keys did the user actually write?", which lives only in the
case file -- and a key nobody catalogued is silently ignored by OpenFOAM.

Measured motivation: misspelling ``activeTensionModel`` in singleCell
produces exit 0, a clean solver log, and an output set quietly missing the
active-tension trace, because the read site is guarded by ``found()``.

Warn, never error: cardiacFoam does not own every key that may legitimately
appear in these dictionaries.
"""

from __future__ import annotations

from pathlib import Path

from openfoam_driver.core.specs.case_dict_keys import case_dict_key_diagnostics

_SINGLE_CELL = (
    Path(__file__).resolve().parents[6]
    / "tutorials" / "electrophysiologyProtocols" / "singleCell"
)

_HEADER = """\
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      electroProperties;
}
"""


def _write(case_root: Path, body: str) -> None:
    d = case_root / "constant"
    d.mkdir(parents=True, exist_ok=True)
    (d / "electroProperties").write_text(_HEADER + body)


def test_misspelled_key_is_reported(tmp_path):
    _write(tmp_path, """
myocardiumSolver singleCellSolver;
singleCellSolverCoeffs
{
    ionicModel      TWorld;
    activeTensionModl LandNiederer;
}
""")
    diags = case_dict_key_diagnostics(
        tmp_path,
        # Catalogue paths are SCOPE-RELATIVE: _parse_path strips the
        # "$ELECTRO_MODEL_COEFFS." prefix, so they describe the inside of the
        # <model>Coeffs dictionary, not the file root.
        catalogued_paths=("myocardiumSolver", "ionicModel", "activeTensionModel"),
        dict_relpaths=("constant/electroProperties",),
    )
    assert [d.field for d in diags] == ["activeTensionModl"]
    assert diags[0].level == "warning"
    assert diags[0].code == "uncatalogued_case_dict_key"


# ---------------------------------------------------------------------------
# Wiring into the strict plan: reported, but never fatal
# ---------------------------------------------------------------------------


def test_runtime_selection_coeffs_dict_is_not_reported(tmp_path):
    """``<model>Coeffs`` is OpenFOAM's runtime-selection convention.

    The catalogue writes these paths with a ``$SCOPE_TOKEN.`` prefix that
    parsing strips, so the literal name never enters the known set. It is
    OpenFOAM's convention, not a cardiacFoam key, and warning about it would
    fire on every single case.
    """
    _write(tmp_path, """
myocardiumSolver singleCellSolver;
singleCellSolverCoeffs
{
    ionicModel TWorld;
}
""")
    diags = case_dict_key_diagnostics(
        tmp_path,
        catalogued_paths=("myocardiumSolver", "ionicModel"),
        dict_relpaths=("constant/electroProperties",),
    )
    assert [d.field for d in diags] == []


def test_a_misspelled_coeffs_dict_is_still_reported(tmp_path):
    """The convention allowance must not swallow an actual typo."""
    _write(tmp_path, """
myocardiumSolver singleCellSolver;
singleCellSolverCoefs
{
    ionicModel TWorld;
}
""")
    diags = case_dict_key_diagnostics(
        tmp_path,
        catalogued_paths=("myocardiumSolver", "ionicModel"),
        dict_relpaths=("constant/electroProperties",),
    )
    assert [d.field for d in diags] == ["singleCellSolverCoefs"]


def test_strict_plan_reports_a_misspelled_key_without_failing(tmp_path):
    """A misspelled key must surface as a warning and leave the plan valid.

    cardiacFoam does not own every key that may appear in these dicts, so an
    unmatched key can never be allowed to fail a plan.
    """
    import shutil

    from openfoam_driver.core.strict_planning import strict_plan

    tutorials_root = tmp_path / "tutorials"
    case = tutorials_root / "case"
    shutil.copytree(_SINGLE_CELL, case)
    ep = case / "constant" / "electroProperties"
    ep.write_text(
        ep.read_text().replace(
            "activeTensionModel LandNiederer;",
            "activeTensionModl LandNiederer;",
        )
    )

    report = strict_plan(
        "case",
        entry_kind="case_folder",
        overrides={"tutorials_root": str(tutorials_root)},
        openfoam_bashrc="/no/such/openfoam/bashrc",
    )
    payload = report.to_json()

    warnings = [
        item
        for item in payload["case_dict_key_diagnostics"]
        if item["code"] == "uncatalogued_case_dict_key"
    ]
    assert [w["field"] for w in warnings] == ["activeTensionModl"]
    assert all(w["level"] == "warning" for w in warnings)

    # Warn-only: the key diagnostics must not appear in any error bucket.
    for bucket in ("validation_diagnostics", "catalog_coverage_errors"):
        assert not [
            item
            for item in payload[bucket]
            if item["code"] == "uncatalogued_case_dict_key"
        ]


# ---------------------------------------------------------------------------
# One matching set, shared with the C++ scanner
# ---------------------------------------------------------------------------


def test_catalogued_names_covers_wildcard_leaves_and_containers():
    """Both directions must agree on what "the catalogue knows" means.

    A second, subtly different copy of this set is what produced the 71%
    false-positive rate on the C++ side.
    """
    from openfoam_driver.core.contracts.dictionary import DictEntry
    from openfoam_driver.scripts._dict_keys_scanner import catalogued_names

    entries = [
        DictEntry(
            driver_path="$C.ecgDomains.<name>.sigmaExtracellular",
            description="",
            phases=frozenset({"physics"}),
            dynamic_path=True,
        ),
        DictEntry(
            driver_path="$C.singleCellSolverCoeffs.ionicModel",
            description="",
            phases=frozenset({"physics"}),
        ),
    ]
    names = catalogued_names(entries)

    assert "sigmaExtracellular" in names, "wildcard-path leaf must count as known"
    assert "ecgDomains" in names, "container segment must count as known"
    assert "singleCellSolverCoeffs" in names
    assert "ionicModel" in names
    assert "<name>" not in names, "wildcard placeholders are not real key names"


# ---------------------------------------------------------------------------
# Position-aware matching: a <placeholder> segment matches any instance name
# ---------------------------------------------------------------------------


def test_user_chosen_instance_name_under_a_wildcard_is_not_reported(tmp_path):
    """``ecgDomains { ECG { ... } }`` -- "ECG" is the author's own label.

    The catalogue models it as ``ecgDomains.<name>.sigmaExtracellular``. A
    flat name-set cannot tell "ECG" from a typo, but the key's position can:
    it sits exactly where the catalogue expects a ``<name>``.
    """
    _write(tmp_path, """
ecgDomains
{
    ECG
    {
        sigmaExtracellular 0.2;
    }
}
""")
    diags = case_dict_key_diagnostics(
        tmp_path,
        catalogued_paths=("ecgDomains.<name>.sigmaExtracellular",),
        dict_relpaths=("constant/electroProperties",),
    )
    assert [d.field for d in diags] == []


def test_a_typo_beneath_a_wildcard_is_still_reported(tmp_path):
    """Wildcard tolerance must not extend to the leaves underneath it."""
    _write(tmp_path, """
ecgDomains
{
    ECG
    {
        sigmaExtracellulr 0.2;
    }
}
""")
    diags = case_dict_key_diagnostics(
        tmp_path,
        catalogued_paths=("ecgDomains.<name>.sigmaExtracellular",),
        dict_relpaths=("constant/electroProperties",),
    )
    assert [d.field for d in diags] == ["sigmaExtracellulr"]


# ---------------------------------------------------------------------------
# Why required_when can never catch this, and the warning is the only signal
# ---------------------------------------------------------------------------


def test_a_misspelled_key_is_silently_replaced_by_the_catalogue_default(tmp_path):
    """The catalogue's own required_when rule cannot catch a misspelling.

    ``stim_amplitude`` is ``required: True`` when ``$singleCellStimulus_present``,
    and that virtual token IS inferred on the read path. The rule is live. It
    still cannot fire, because the builder fills every entry's
    ``typical_value`` BEFORE validation runs -- so requiredness is satisfied
    by construction, and any key carrying a typical_value is structurally
    immune to the required-field check.

    The visible consequence: an author writing ``stim_amplitud 25`` gets
    ``stim_amplitude 60`` -- their value discarded, the catalogue default
    silently substituted. Nothing in the required-field machinery objects.
    The uncatalogued-key warning is the only signal anywhere in the system,
    which is precisely why it is worth having at warn level.
    """
    import shutil

    from openfoam_driver.plugins.cardiacfoam import dict_builder as DB
    from openfoam_driver.core.strict_planning import strict_plan

    tutorials_root = tmp_path / "tutorials"
    case = tutorials_root / "case"
    shutil.copytree(_SINGLE_CELL, case)
    ep = case / "constant" / "electroProperties"
    ep.write_text(ep.read_text().replace("stim_amplitude  60;", "stim_amplitud  25;"))

    parsed = DB.parse_electro_properties(ep)
    rebuilt = str(
        DB.build_electro_properties(
            selectors=parsed["selectors"], overrides=parsed["overrides"]
        )
    )
    amplitude = [l.strip() for l in rebuilt.splitlines() if "stim_amplitude" in l]
    assert amplitude == ["stim_amplitude 60;"], (
        f"expected the catalogue default to be substituted, got {amplitude}"
    )
    assert "stim_amplitud " not in rebuilt, "the misspelled key is dropped entirely"

    payload = strict_plan(
        "case",
        entry_kind="case_folder",
        overrides={"tutorials_root": str(tutorials_root)},
        openfoam_bashrc="/no/such/openfoam/bashrc",
    ).to_json()

    # The plan is valid -- the built dict really is complete and correct.
    assert payload["status"] == "ok"
    assert not [
        item
        for item in payload["validation_diagnostics"]
        if "stim_amplitude" in item["message"]
    ], "required_when cannot fire here; if it starts to, this test should change"

    # ...and the warning is the only thing that noticed.
    assert [
        item["field"] for item in payload["case_dict_key_diagnostics"]
    ] == ["stim_amplitud"]
