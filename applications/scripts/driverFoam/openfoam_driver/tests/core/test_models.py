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
#     test_models
#
# Description
#     Tests for core runtime model dataclasses (TutorialSpec, etc.).
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

from __future__ import annotations

import dataclasses

from openfoam_driver.core.runtime.models import TutorialSpec


def test_tutorial_spec_no_longer_has_a_collect_outputs_field():
    # collect_outputs was a TutorialSpec field nothing in the runtime ever
    # called (confirmed: zero call sites for .collect_outputs anywhere,
    # zero test references) -- removed together with every tutorial's
    # per-case assignment to it.
    field_names = {f.name for f in dataclasses.fields(TutorialSpec)}
    assert "collect_outputs" not in field_names
