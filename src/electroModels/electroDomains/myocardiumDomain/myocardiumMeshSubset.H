/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    The heart sub-mesh: the myocardium cellZone of the support mesh and the
    patch its exposed faces join. Shared by the myocardium domains and the
    lead-field conductor, so that all build the same heart.

SourceFiles
    myocardiumMeshSubset.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef myocardiumMeshSubset_H
#define myocardiumMeshSubset_H

#include "fvMeshSubset.H"
#include "dictionary.H"
#include "autoPtr.H"

namespace Foam
{

//- Support-mesh patch that receives the faces a myocardium cellZone exposes.
//  electroProperties' exposedFacesPatch names it; otherwise the first patch
//  that is neither empty nor coupled. OpenFOAM's own default for exposed
//  faces is an empty patch, which gradients and Laplacian corrections skip.
label myocardiumExposedFacesPatch
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties
);

//- Sub-mesh of the myocardium cellZone named by cellZone in
//  electroProperties (null without cellZone)
autoPtr<fvMeshSubset> createMyocardiumMeshSubset
(
    const fvMesh& supportMesh,
    const dictionary& electroProperties
);

} // End namespace Foam

#endif

// ************************************************************************* //
