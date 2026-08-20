"""Differential harness: post-change mutators vs. the tier-1-only reference.

The one-time byte-identity and semantic-preservation comparisons that used to
live here (against ``_mutators_reference.py``, the pre-migration
implementation) have been run and verified across the real ``tutorials/``
corpus and removed per plan. What remains is the standing security
regression: a directive embedded in a dict value must never be evaluated.
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
