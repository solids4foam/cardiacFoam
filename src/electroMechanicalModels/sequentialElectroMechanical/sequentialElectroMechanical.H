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

Class
    sequentialElectroMechanical

Description
    Sequential (weak) electro-mechanical coupling, where no outer correctors
    are performed each time-step. The electro model is solved first, then the
    solid model.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    sequentialElectroMechanical.C

\*---------------------------------------------------------------------------*/

#ifndef sequentialElectroMechanical_H
#define sequentialElectroMechanical_H

#include "electroMechanicalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroMechanicalModels
{

/*---------------------------------------------------------------------------*\
               Class sequentialElectroMechanical Declaration
\*---------------------------------------------------------------------------*/

class sequentialElectroMechanical
:
    public electroMechanicalModel
{
    // Private data

        //- Active tension field registered on the solid mesh.
        //  Populated each time step from the electro model's Cai field
        //  and picked up by electroMechanicalLaw via objectRegistry lookup.
        volScalarField Ta_;

        //- Scaling factor: Ta = kTa * max(Cai - CaiThreshold, 0)
        //  Units: [Pa / Cai-unit]. The ionic model returns Cai in its
        //  native units (e.g. mM for TNNP). This factor converts to stress.
        //  This simple linear model is a placeholder; it will be replaced by
        //  a dedicated runtime-selectable active tension model in the future.
        const dimensionedScalar kTa_;

        //- Threshold Cai below which no active tension is generated
        const dimensionedScalar CaiThreshold_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sequentialElectroMechanical(const sequentialElectroMechanical&);

        //- Disallow default bitwise assignment
        void operator=(const sequentialElectroMechanical&);


public:

    //- Runtime type information
    TypeName("sequentialElectroMechanical");


    // Constructors

        //- Construct from components
        sequentialElectroMechanical
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~sequentialElectroMechanical()
        {}


    // Member Functions

        // Edit

            //- Evolve the electro-mechanical model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroMechanicalModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
