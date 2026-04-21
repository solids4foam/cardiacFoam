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
    electroMechanicalModel

Description
    Virtual base class for electro-mechanical coupling.

    Follows the fluidSolidInterface pattern from solids4foam: owns an
    electroModel and a solidModel in separate regions, and delegates to
    derived classes for the coupling scheme.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    electroMechanicalModel.C
\*---------------------------------------------------------------------------*/

#ifndef electroMechanicalModel_H
#define electroMechanicalModel_H

#include "physicsModel.H"
#include "IOdictionary.H"
#include "dynamicFvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "electroModel.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class electroMechanicalModel Declaration
\*---------------------------------------------------------------------------*/

class electroMechanicalModel
:
    public physicsModel,
    public IOdictionary
{
    // Private data

        //- Electro-mechanical properties dictionary
        dictionary& electroMechanicalProperties_;

        //- Electro solver
        autoPtr<electroModel> electro_;

        //- Solid solver
        autoPtr<solidModel> solid_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        electroMechanicalModel(const electroMechanicalModel&);

        //- Disallow default bitwise assignment
        void operator=(const electroMechanicalModel&);


protected:

    // Protected Member Functions

        //- Return sub-dictionary
        dictionary& electroMechanicalProperties()
        {
            return electroMechanicalProperties_;
        }


public:

    //- Runtime type information
    TypeName("electroMechanicalModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            electroMechanicalModel,
            dictionary,
            (
                Time& runTime,
                const word& region
            ),
            (runTime, region)
        );


    // Constructors

        //- Construct from components
        electroMechanicalModel
        (
            const word& type,
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Selectors

        //- Select constructed from time
        static autoPtr<electroMechanicalModel> New
        (
            Time& runTime,
            const word& region
        );


    // Destructor

        virtual ~electroMechanicalModel();


    // Member Functions

        // Access

            //- Return time
            const Time& runTime() const
            {
                return electro_().runTime();
            }

            //- Return sub-dictionary
            const dictionary& electroMechanicalProperties() const
            {
                return electroMechanicalProperties_;
            }

            //- Return electro model
            const electroModel& electro() const
            {
                return electro_();
            }

            //- Return electro model
            electroModel& electro()
            {
                return electro_();
            }

            //- Return solid model
            const solidModel& solid() const
            {
                return solid_();
            }

            //- Return solid model
            solidModel& solid()
            {
                return solid_();
            }

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the electro-mechanical model
            virtual bool evolve() = 0;

            //- Write fields
            virtual void writeFields(const Time& runTime);

            //- This function is called at the end of the simulation
            virtual void end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
