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
#     test_fresh
#
# Description
#     Tests the --fresh safety guards and deletion in core/runtime/fresh.py.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from pathlib import Path

from openfoam_driver.core.runtime.fresh import (
    check_fresh_deletion_allowed,
    ensure_fresh_output_dir,
)


def test_allows_deletion_of_nonexistent_directory(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    assert check_fresh_deletion_allowed(target, allowed_root=None) is None


def test_allows_deletion_of_directory_with_top_level_marker(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    target.mkdir(parents=True)
    (target / "workflow_state.json").write_text("{}")
    assert check_fresh_deletion_allowed(target, allowed_root=None) is None


def test_allows_deletion_when_marker_is_one_level_down(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    case_dir = target / "TNNP"
    case_dir.mkdir(parents=True)
    (case_dir / "run_document.json").write_text("{}")
    assert check_fresh_deletion_allowed(target, allowed_root=None) is None


def test_refuses_directory_with_no_driverfoam_marker(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    target.mkdir(parents=True)
    (target / "notes.txt").write_text("do not delete me")
    error = check_fresh_deletion_allowed(target, allowed_root=None)
    assert error is not None
    assert "no recognizable driverFOAM artifact" in error


def test_refuses_filesystem_root():
    error = check_fresh_deletion_allowed(Path("/"), allowed_root=None)
    assert error is not None
    assert "filesystem root" in error


def test_refuses_home_directory():
    error = check_fresh_deletion_allowed(Path.home(), allowed_root=None)
    assert error is not None
    assert "home directory" in error


def test_refuses_path_with_two_segments():
    error = check_fresh_deletion_allowed(Path("/a/b"), allowed_root=None)
    assert error is not None
    assert "too shallow" in error


def test_allows_path_with_three_segments():
    assert check_fresh_deletion_allowed(Path("/a/b/c"), allowed_root=None) is None


def test_refuses_path_outside_allowed_root(tmp_path):
    allowed_root = tmp_path / "allowed"
    allowed_root.mkdir()
    outside = tmp_path / "elsewhere" / "case" / "output"
    outside.mkdir(parents=True)
    (outside / "workflow_state.json").write_text("{}")
    error = check_fresh_deletion_allowed(outside, allowed_root=allowed_root)
    assert error is not None
    assert "DRIVERFOAM_ALLOWED_RUNS_ROOT" in error


def test_allows_path_inside_allowed_root(tmp_path):
    allowed_root = tmp_path / "allowed"
    inside = allowed_root / "case" / "output"
    inside.mkdir(parents=True)
    (inside / "workflow_state.json").write_text("{}")
    assert check_fresh_deletion_allowed(inside, allowed_root=allowed_root) is None


def test_ensure_fresh_output_dir_deletes_when_allowed(tmp_path, capsys):
    target = tmp_path / "a" / "b" / "c"
    target.mkdir(parents=True)
    (target / "workflow_state.json").write_text('{"status": "completed"}')
    (target / "stale_result.txt").write_text("old")

    error = ensure_fresh_output_dir(target, fresh=True, allowed_root=None)

    assert error is None
    assert not target.exists()
    captured = capsys.readouterr()
    assert "deleting" in captured.err
    assert captured.out == ""


def test_ensure_fresh_output_dir_noop_when_fresh_false(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    target.mkdir(parents=True)
    (target / "workflow_state.json").write_text("{}")

    error = ensure_fresh_output_dir(target, fresh=False, allowed_root=None)

    assert error is None
    assert target.exists()


def test_ensure_fresh_output_dir_leaves_directory_untouched_on_refusal(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    target.mkdir(parents=True)
    (target / "notes.txt").write_text("keep me")

    error = ensure_fresh_output_dir(target, fresh=True, allowed_root=None)

    assert error is not None
    assert target.exists()
    assert (target / "notes.txt").read_text() == "keep me"


def test_ensure_fresh_output_dir_noop_when_directory_does_not_exist(tmp_path):
    target = tmp_path / "a" / "b" / "c"
    assert ensure_fresh_output_dir(target, fresh=True, allowed_root=None) is None
    assert not target.exists()
