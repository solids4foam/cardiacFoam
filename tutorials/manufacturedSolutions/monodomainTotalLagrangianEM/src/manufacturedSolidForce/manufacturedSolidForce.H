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
    Foam::fv::manufacturedSolidForce

Description
    fvOption that applies a manufactured body force to the solid momentum
    equation for electromechanics manufactured-solution verification.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedSolidForce_H
#define manufacturedSolidForce_H

#include "fvOption.H"

namespace Foam
{
namespace fv
{

class manufacturedSolidForce
:
    public option
{
    vector amplitude_;
    scalar Tmax_;
    scalar V0_;
    scalar gamma_;
        //- Passive mechanical parameters
        scalar E_;
        scalar nu_;

        //- Density
        scalar rho0_;

public:

    TypeName("manufacturedSolidForce");

    manufacturedSolidForce
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~manufacturedSolidForce();

    virtual bool isActive() { return true; }

    virtual void addSup
    (
        fvMatrix<vector>& eqn,
        const label fieldi
    );

    virtual bool read(const dictionary& dict);
};

}
}

#endif
