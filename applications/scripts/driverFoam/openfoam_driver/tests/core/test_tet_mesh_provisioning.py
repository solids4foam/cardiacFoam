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
#     test_tet_mesh_provisioning
#
# Description
#     Tests __LC__ rendering for tet-mesh sweeps. Never calls gmsh/gmshToFoam
#     -- those run later as workflow_dag steps, not during materialization.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path

import pytest

from openfoam_driver.core.specs.tet_mesh_provisioning import render_tet_geo


def _write_template(root: Path, *, relpath: str = "setup/studies/tetConvergence/box.geo.template", body: str = "lc = __LC__;\nBox(1) = {0, 0, 0, 1, 1, 1};\n") -> None:
    path = root / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(body)


def test_renders_lc_from_n(tmp_path):
    _write_template(tmp_path)
    geo_path = render_tet_geo(tmp_path, 10)
    text = geo_path.read_text()
    assert "lc = 0.1;" in text
    assert "__LC__" not in text


def test_returns_the_written_geo_path(tmp_path):
    _write_template(tmp_path)
    geo_path = render_tet_geo(tmp_path, 20)
    assert geo_path == tmp_path / "setup/studies/tetConvergence/box.geo"
    assert geo_path.exists()


def test_rejects_non_positive_n(tmp_path):
    _write_template(tmp_path)
    with pytest.raises(ValueError, match="positive"):
        render_tet_geo(tmp_path, 0)
    with pytest.raises(ValueError, match="positive"):
        render_tet_geo(tmp_path, -5)


def test_rejects_non_integer_n(tmp_path):
    _write_template(tmp_path)
    with pytest.raises(TypeError, match="integer"):
        render_tet_geo(tmp_path, 10.5)


def test_missing_template_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        render_tet_geo(tmp_path, 10)


def test_rejects_template_with_no_lc_placeholder(tmp_path):
    _write_template(tmp_path, body="Box(1) = {0, 0, 0, 1, 1, 1};\n")
    with pytest.raises(ValueError, match="__LC__"):
        render_tet_geo(tmp_path, 10)


def test_rejects_template_with_multiple_lc_placeholders(tmp_path):
    _write_template(tmp_path, body="lc = __LC__;\nlc2 = __LC__;\n")
    with pytest.raises(ValueError, match="__LC__"):
        render_tet_geo(tmp_path, 10)


def test_comment_mentioning_the_placeholder_by_name_does_not_count(tmp_path):
    # The real box.geo.template files document the substitution mechanism in
    # a `//` comment that names the placeholder literally ("__LC__ is
    # substituted by run_corrector_study.sh") -- found via real, non-mocked
    # verification against an actual bidomain scratch copy. A naive
    # text.count("__LC__") treats that as a second occurrence and rejects a
    # perfectly valid template; only the real substitution site should count,
    # and the comment's own textual mention must survive untouched (a blind
    # str.replace would otherwise mangle it into nonsense with the number).
    _write_template(
        tmp_path,
        body="// The characteristic length placeholder __LC__ is substituted here.\nlc = __LC__;\n",
    )
    geo_path = render_tet_geo(tmp_path, 10)
    text = geo_path.read_text()
    assert "lc = 0.1;" in text
    assert "// The characteristic length placeholder __LC__ is substituted here." in text
