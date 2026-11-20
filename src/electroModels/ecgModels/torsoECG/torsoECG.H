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
    torsoECG

Description
    ECG of the electrodes in the torso, the region outside the heart
    (the torso block of electroProperties), by torsoCoupling:
    - bath: the phiE of the bath domain, solved with the heart, in the
      electrode cells;
    - leadField: the lead fields of buildLeadFields contracted with the
      heart source of the myocardium's Vm (leadFieldECG).

    Every electrode lies in the mesh: an electrode outside it is fatal
    (ecgElectrodeCells).

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef torsoECG_H
#define torsoECG_H

#include "ecgSolver.H"
#include "leadFieldECG.H"
#include "labelList.H"
#include "overrideTypeName.H"

namespace Foam
{

class torsoECG
:
    public ecgSolver
{
    //- Lead-field evaluation (torsoCoupling leadField), else null
    autoPtr<leadFieldECG> leadFieldPtr_;

    //- Electrode cells on the phiE mesh (torsoCoupling bath)
    labelList electrodeCells_;
    bool cellsBuilt_;

    //- phiE of the bath domain in the electrode cells
    void sample(const ecgDomain& domain, scalarField& values);

public:

    OverrideTypeName("torsoECG");

    //- Construct from the ECG domain entry of the ecgDomains block of
    //  electroProperties, whose torsoCoupling selects the evaluation
    explicit torsoECG(const dictionary& dict);

    //- The lead fields it evaluates (torsoCoupling leadField; loaded at the
    //  first solve), else null
    const leadFieldECG* leadField() const
    {
        return leadFieldPtr_.get();
    }

    //- The bath's phiE with torsoCoupling bath, else the myocardium's Vm
    virtual stateSource source() const override
    {
        return leadFieldPtr_ ? stateSource::myocardium : stateSource::bath;
    }

    virtual void solve
    (
        ecgDomain& domain,
        scalar t0,
        scalar dt,
        scalarField& values
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
