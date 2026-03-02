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
    singleCellElectro

Description
    Single-cell electrophysiology solver based on OpenFOAM.

    This solver advances a single cardiac cell model in time using a
    run-time selectable ionic model. In contrast to electroActivationFoam,
    no spatial voltage PDE is solved; instead, the transmembrane voltage
    Vm is integrated directly as part of the ionic model ODE system.

    The solver is primarily intended for:
      - single-cell action potential simulations,
      - testing and validation of ionic models,
      - debugging and development of electrophysiology model components.

    External stimulation is prescribed via a stimulus protocol dictionary.
    Output is written as time series of states, algebraic variables and rates
    using ionicModelIO.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    singleCellElectro.C

\*---------------------------------------------------------------------------*/

#ifndef singleCellElectro_H
#define singleCellElectro_H

#include "electroModel.H"
#include "ionicModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

/*---------------------------------------------------------------------------* \
                        Class singleCellElectro Declaration
\*---------------------------------------------------------------------------*/

class singleCellElectro
:
    public electroModel
{
    // Private data

        //- Cardiac properties dictionary
        const dictionary& cardiacProperties_;

        //- Run-time selectable ionic model
        autoPtr<ionicModel> ionicModelPtr_;

        //- Output file pointer
        autoPtr<OFstream> outputPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        singleCellElectro(const singleCellElectro&);

        //- Disallow default bitwise assignment
        void operator=(const singleCellElectro&);


public:

    //- Runtime type information
    TypeName("singleCellElectro");

    // Constructors

        //- Construct from components
        singleCellElectro
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~singleCellElectro()
        {}


    // Member Functions

        // Edit

            //- Evolve the electro model
            virtual bool evolve();

            //- This function is called at the end of the simulation
            virtual void end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
