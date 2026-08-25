# driverFOAM scripts/ audit — what was left undone (2026-08-25)

Written at the end of the `applications/scripts/driverFoam/scripts/` audit
(commits `98e87db9`..`8528041c` plus the sweep-contract fix). Everything the
audit *did* is in those commits and their messages; this file records only
what was deliberately **not** done, and why, so nobody re-derives it.

Status of the tree at time of writing: suite green at 1701 passed, working
tree clean apart from the pre-existing `modules/solids4foam` submodule.

---

## 1. `ecgModelIO::loadSurface` is dead code

**Where:** `src/genericWriter/ecgModelIO.C:71` (declared at `ecgModelIO.H:88`).

**What:** It performs a real `dict.get<fileName>("torsoSurface")` read — an
optional key pointing at an STL torso surface — but the function is
**declared, defined, and called from nowhere in `src/`**. Verified by grep
across the whole tree, excluding `lnInclude/`.

**Consequence today:** `torsoSurface` shows up in the C++ dictionary-key
scanner as an unmatched read, and is waived in
`openfoam_driver/plugins/cardiacfoam/dict_key_allowlist.json` under
`unmatched_cxx_reads`, tagged group (d) DEAD CODE.

**Decision needed:** either delete `loadSurface` and drop the waiver, or wire
it up — in which case `torsoSurface` becomes a genuine catalogue gap and the
waiver **must** be removed and a `DictEntry` added. The allowlist description
already carries that instruction.

**Cost:** small either way. This is a C++ decision, not a driverFOAM one.

---

## 2. The C++ dict-key scanner cannot see reads through a call receiver

**Where:** `applications/scripts/driverFoam/openfoam_driver/scripts/_dict_keys_scanner.py`

**What:** The scanner regex matches `<identifier>.lookup("key")` and friends.
It does **not** match a receiver that is itself a call expression:

```cpp
electroMechanicalProperties().getOrDefault<scalar>("someKey", 0.5)   // INVISIBLE
const dictionary& d = electroMechanicalProperties();
d.getOrDefault<scalar>("someKey", 0.5);                              // seen
```

Confirmed by injecting both forms into a real solver file: the first was not
detected, the second was, and `plan --strict` failed on the second with
`plugin_dict_key_unmatched_cxx_reads`.

**Consequence:** this is the concrete reason the scanner's own docstring says
accuracy is ~80%, and the reason `--strict` is allowlist-backed rather than
trusted outright. A new C++ tuning knob added through a call receiver will
not be flagged as uncatalogued.

**Fix:** a real C++ parser (libclang), not a bigger regex. Regex cannot
resolve what a call expression returns. This is a genuine project, not a
patch — weigh it against how often keys are actually added this way.

---

## 3. Catalogue `required` cannot catch a misspelling — by construction

**Not a bug to fix. Recorded so it is not rediscovered as one.**

`stim_amplitude` is `required: True` when `$singleCellStimulus_present`, and
that virtual token **is** correctly inferred on the read path (an earlier
diagnosis that it was dormant was wrong). The rule is live. It still cannot
fire, because the builder fills every entry's `typical_value` **before**
validation runs, so requiredness is satisfied by construction.

**Any catalogue key carrying a `typical_value` is structurally immune to the
required-field check.**

Demonstrated: writing `stim_amplitud 25` makes driverFOAM build
`stim_amplitude 60` — the author's value discarded, the catalogue default
silently substituted, zero required-field diagnostics. Pinned by
`test_a_misspelled_key_is_silently_replaced_by_the_catalogue_default`.

This is precisely why the warn-level uncatalogued-key sweep
(`core/specs/case_dict_keys.py`) exists: it is the **only** mechanism in the
stack that can see this class. If someone later "fixes" required_when to fire
here, that test should be revisited deliberately, not deleted.

---

## 4. Explicitly out of scope: C++ ↔ catalogue drift checking

Considered and **declined by the maintainer** on 2026-08-25.

The idea was a test asserting every hand-written C++ requirement (e.g. the
four-key completeness guard at `src/genericWriter/stimulusIO.C:160-175`)
corresponds to a catalogue `required` / `required_when`. Those two really are
the same rule written twice in two languages.

Rejected because the catalogue should stay **as agnostic as possible to logic
living inside C++**. Do not resurrect this without re-asking. The C++ guards
stay as they are; they also protect people invoking `cardiacFoam` directly
without driverFOAM, so they cannot simply be deleted either.

---

## 5. Housekeeping

- **73 commits ahead of `github/ep-work-onto-main`**, nothing pushed.
- A **concurrent session** was editing this repo throughout the audit
  (`CITATION.cff`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`,
  `VERSION_POLICY.md`, the `monodomain1D3D` tutorial READMEs, the
  `buildAndTest.yml` foamlib rename, and a `pytest-cov` dev-dependency
  group). All committed here with attribution, in `7f547301` and `3ca67bcb`.
  Check with whoever drove that session before pushing.
- Two cells in the new `tutorials/README.md` `monodomain1D3D` row were
  completing the documentation contract incorrectly and left two tests red;
  corrected in `7f547301`.
