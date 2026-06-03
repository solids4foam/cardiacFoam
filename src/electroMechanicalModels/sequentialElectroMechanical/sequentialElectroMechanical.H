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
#include "activeTensionModel.H"
#include "overrideTypeName.H"

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
        //  Populated each time step by activeTensionModel::calculateTension()
        //  and picked up by electroMechanicalLaw via objectRegistry lookup.
        volScalarField Ta_;

        //- Per-cell stretch field passed to calculateTension().
        //  Phase 1: uniform 1.0 (no solid feedback).
        //  Phase 2: lambda = sqrt(f0.C.f0) from solid deformation gradient.
        scalarField lambdaField_;

        //- Runtime-selectable active tension model
        autoPtr<activeTensionModel> activeTensionModel_;

        //- Flag: print solid objectRegistry probe on first evolve() call only
        bool firstTimeStep_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sequentialElectroMechanical(const sequentialElectroMechanical&);

        //- Disallow default bitwise assignment
        void operator=(const sequentialElectroMechanical&);


public:

    //- Runtime type information
    OverrideTypeName("sequentialElectroMechanical");


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
            virtual bool evolve() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroMechanicalModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
