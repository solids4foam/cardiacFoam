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
    leadFieldECG

Description
    ECG from lead fields: the potential of electrode e is

        phi_e = sum over heart cells of z_e V b(Vm) + c_e

    where z_e is the lead field written by buildLeadFields in
    constant/[region/]leadField_<e> (per processor in parallel), b is the
    forward heart source of passiveConductor for the myocardium's Vm, and
    c_e is the constant contribution of ground values, surface currents and
    the reference value.

    The conductor is built on the myocardium's sub-mesh at the first solve.
    Its signature must match constant/[region/]leadFieldSignature.

SourceFiles
    leadFieldECG.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef leadFieldECG_H
#define leadFieldECG_H

#include "ecgSolver.H"
#include "electroStateProvider.H"
#include "leadFieldConductor.H"
#include "overrideTypeName.H"

namespace Foam
{

class leadFieldECG
:
    public ecgSolver
{
    //- Myocardium solver Coeffs holding the ECG domain entry
    const dictionary electroCoeffs_;

    const word ecgDomainName_;

    autoPtr<leadFieldConductor> conductorPtr_;

    //- z_e V on the heart cells, per electrode
    List<scalarField> heartWeights_;

    //- Constant term c_e, per electrode
    scalarList constants_;

public:

    OverrideTypeName("leadFieldECG");

    //- Construct from the ECG domain entry of the myocardium solver Coeffs
    explicit leadFieldECG(const dictionary& dict);

    //- Construct from the myocardium solver Coeffs and the ECG domain name
    leadFieldECG(const dictionary& electroCoeffs, const word& ecgDomainName);

    //- Build the conductor on the provider's myocardium, check the
    //  signature and read the lead fields of the named electrodes
    void load
    (
        const electroStateProvider& provider,
        const wordList& electrodeNames
    );

    bool loaded() const
    {
        return conductorPtr_.valid();
    }

    const leadFieldConductor& conductor() const
    {
        return conductorPtr_();
    }

    //- Channels sum z_e V b(Vm) + c_e for a heart Vm on the conductor's
    //  heart mesh
    void evaluate(const volScalarField& Vm, scalarField& values) const;

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
