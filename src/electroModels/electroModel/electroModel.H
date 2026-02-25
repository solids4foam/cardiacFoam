/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    electroModel

Description
    Virtual base class for electrophysiology solvers.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    electroModel.C

\*---------------------------------------------------------------------------*/

#ifndef electroModel_H
#define electroModel_H

#include "physicsModel.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pimpleControl.H"
#include "dimVoltage.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class electroModel Declaration
\*---------------------------------------------------------------------------*/

class electroModel
:
    public physicsModel,
    public IOdictionary
{
    // Private data

        //- Mesh
        autoPtr<dynamicFvMesh> meshPtr_;

        //- Electro properties dictionary
        dictionary electroProperties_;

        //- Pimple control for monitoring the convergence of the outer loop
        mutable autoPtr<pimpleControl> pimplePtr_;


    // Private Member Functions

        //- Make the pimple control
        void makePimpleControl() const;

        //- Disallow default bitwise copy construct
        electroModel(const electroModel&);

        //- Disallow default bitwise assignment
        void operator=(const electroModel&);


protected:

    // Protected member functions

        //- Return electro properties dictionary
        dictionary& electroProperties()
        {
            return electroProperties_;
        }

        //- Return electro properties dictionary
        const dictionary& electroProperties() const
        {
            return electroProperties_;
        }

        //- Return non-const acess to the pimple control
        pimpleControl& pimple();

public:

    //- Runtime type information
    TypeName("electroModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            electroModel,
            dictionary,
            (
                Time& runTime,
                const word& region
            ),
            (runTime, region)
        );


    // Constructors

        //- Construct from components
        electroModel
        (
            const word& type,
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<electroModel> New
        (
            Time& runTime, const word& region = "region0"
        );


    // Destructor

        virtual ~electroModel();


    // Member Functions

        // Access

            //- Return reference to the mesh
            dynamicFvMesh& mesh()
            {
                return meshPtr_();
            }

            //- Return const reference to the mesh
            const dynamicFvMesh& mesh() const
            {
                return meshPtr_();
            }

            //- Return time
            const Time& runTime() const
            {
                return meshPtr_().time();
            }

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime)
            {
                physicsModel::setDeltaT(runTime);
            }

            //- Evolve the electro model
            virtual bool evolve() = 0;

            //- Read electroProperties dictionary
            virtual bool read();

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
