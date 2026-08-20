"""Differential harness: post-change mutators vs. the tier-1-only reference.

The reference module and the two tests comparing against it
(`test_tier1_output_is_byte_identical`, `test_tier2_case_reference_fails_and_
new_code_succeeds`) were a one-time regression check run against the
pre-migration `mutators.py` (captured from git history as
`_mutators_reference.py`) across the real tutorial corpus, then deleted once
the migration was verified clean. `test_no_directive_is_evaluated` stands
alone and is worth keeping permanently.
"""

import shutil

from openfoam_driver.core.runtime import mutators


def test_no_directive_is_evaluated(tmp_path):
    """A #codeStream entry must survive as text, never execute."""
    path = tmp_path / "d"
    path.write_text(
        "FoamFile { version 2.0; class dictionary; object d; }\n"
        "pwned  #codeStream { code #{ os << system(\"touch /tmp/PWNED_DIFF\"); #}; };\n"
        "sigma  0.2;\n"
    )
    mutators.update_foam_entry(path, "sigma", 0.35)
    assert "#codeStream" in path.read_text()
    assert not shutil.os.path.exists("/tmp/PWNED_DIFF")
