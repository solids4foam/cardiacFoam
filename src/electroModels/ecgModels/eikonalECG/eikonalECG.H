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
    eikonalECG

Description
    Template-voltage surrogate ECG solver for eikonal activation-time fields.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef eikonalECG_H
#define eikonalECG_H

#include "ecgSolver.H"
#include "OFstream.H"
#include "Switch.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "overrideTypeName.H"
#include "volFields.H"
#include "eikonalTemplateSource.H"

namespace Foam
{

class eikonalECG
:
    public ecgSolver
{
    eikonalTemplateSource templates_;
    scalar startTime_;
    scalar endTime_;
    scalar deltaT_;
    scalar sigmaE_;
    scalarField lastValues_;
    autoPtr<OFstream> outputPtr_;
    bool sampled_;
    bool leadVectorsCalculated_;
    List<List<vector>> leadVectors_;
    void calculateLeadVectors(const ecgDomain& domain);

    void calculatePseudoECG
    (
        const ecgDomain& domain,
        const volVectorField& gradVm,
        scalarField& values
    ) const;

public:

    OverrideTypeName("eikonalECG");

    explicit eikonalECG(const dictionary& dict);

    virtual bool writesOwnTimeSeries() const override
    {
        return true;
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
