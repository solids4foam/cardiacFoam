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

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  NOT_COMPILED
  This file is not included in Make/files. Bath ECG coupling is not currently
  wired into any active orchestration path. Future extension of BathDomain.
\*---------------------------------------------------------------------------*/

#ifndef bidomainBathECGSolver_H
#define bidomainBathECGSolver_H

#include "Switch.H"
#include "bathECGSolver.H"
#include "overrideTypeName.H"

namespace Foam
{

class bidomainBathECGSolver
:
    public BathECGSolver
{
    label phiEReferenceCell_;
    scalar phiEReferenceValue_;
    Switch reportBathExtrema_;

public:

    OverrideTypeName("bidomainBathECG");

    explicit bidomainBathECGSolver(const dictionary& dict);

    virtual void solve
    (
        BathDomain& domain,
        scalar t0,
        scalar dt
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
