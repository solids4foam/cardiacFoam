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
#     provenance_dependencies
#
# Description
#     Turns a plugin-declared RuntimeDependency into a fingerprinted
#     ProvenanceComponent, reusing Task 1's provenance model (component_for_
#     path's content-hash / degrade-to-metadata / unavailable policy)
#     instead of building a parallel one.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Compose ``RuntimeDependency`` declarations into ``ProvenanceComponent``s.

A runtime dependency's resolved path is not under any case root -- a library
sits in ``$FOAM_USER_LIBBIN``, a case-local ``.so`` sits under the case's own
``platforms/`` tree, and the solver binary sits wherever ``PATH`` says.
``component_for_path`` fingerprints relative to a root the path is *under*,
so this module fingerprints each dependency relative to its own parent
directory and then restores the declared dependency name as the component's
identity -- a library found via a different search directory, or under a
different platform extension, must still compare as the same dependency
rather than reading as added+removed.

Every dependency always produces a component, present or not. An
unresolved dependency (``path is None``) always degrades to
``strength="unavailable"``, whether or not it was ``required``: reporting
absence, never omitting it, is the fix for the incident this module exists
to prevent (I3) -- a required-but-missing dependency must surface, not
silently vanish the way the old ``tuple[Path, ...]`` contract allowed. An
unavailable *optional* dependency (e.g. libelectroMechanicalModels absent in
the maintainer's lightweight default build) is not an error, but it still
shows up and still makes any snapshot built from it partial per I5 -- "the
honest outcome" the plan calls for. What a resume does with that partiality
is a later task's resume-policy question, not this module's.
"""

from __future__ import annotations

from dataclasses import replace

from ..plugin_capabilities import RuntimeDependency
from .provenance import ProvenanceComponent, component_for_path


def component_for_runtime_dependency(dependency: RuntimeDependency) -> ProvenanceComponent:
    """Fingerprint one runtime dependency, present or not."""
    if dependency.path is None:
        return ProvenanceComponent(
            kind="runtime_dependency",
            path=dependency.name,
            role="required_input",
            method="unavailable",
            strength="unavailable",
        )
    component = component_for_path(
        dependency.path,
        kind="runtime_dependency",
        relative_to=dependency.path.parent,
    )
    return replace(component, path=dependency.name)
