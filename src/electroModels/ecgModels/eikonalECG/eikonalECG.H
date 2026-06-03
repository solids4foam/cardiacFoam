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
#include "overrideTypeName.H"
#include "volFields.H"

namespace Foam
{

class eikonalECG
:
    public ecgSolver
{
    scalarField templateTimes_;
    scalarField templateValues_;
    Switch useManufacturedTemplate_;

    scalar startTime_;
    scalar endTime_;
    scalar deltaT_;
    Switch report_;

    bool written_;
    scalarField lastValues_;
    autoPtr<OFstream> outputPtr_;
    autoPtr<volScalarField> VmPtr_;

    // Optimization: precomputed transfer vectors for pseudo-ECG
    bool leadVectorsCalculated_;
    List<List<vector>> leadVectors_;

    bool weightsCalculated_;
    scalarField wEndo_;
    scalarField wMid_;
    scalarField wEpi_;

    scalar manufacturedTemplateValue(scalar localTime) const;

    void calculateLeadVectors(const ecgDomain& domain);
    void calculateTransmuralWeights(const ecgDomain& domain);

    void reconstructVm
    (
        scalar sampleTime,
        const volScalarField& activationTime,
        volScalarField& Vm
    ) const;

    void calculatePseudoECG
    (
        const ecgDomain& domain,
        const volScalarField& Vm,
        scalarField& values
    ) const;

    volScalarField& surrogateVm(ecgDomain& domain);

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
