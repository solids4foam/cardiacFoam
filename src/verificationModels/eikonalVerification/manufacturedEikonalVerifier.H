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
    manufacturedEikonalVerifier

Description
    Manufactured planar-wave verifier for eikonal activation time.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedEikonalVerifier_H
#define manufacturedEikonalVerifier_H

#include "eikonalVerificationModel.H"
#include "overrideTypeName.H"
#include "Switch.H"
#include "dimensionedScalar.H"
#include "dictionary.H"
#include "volFields.H"

namespace Foam
{

class manufacturedEikonalVerifier : public eikonalVerificationModel
{
    const fvMesh& mesh_;
    const volTensorField& conductivity_;
    scalar chi_;
    scalar Cm_;
    scalar c0_;

    //- Write the cellwise signed error as a volScalarField alongside the
    //  aggregate norms. Off by default, so existing cases are unchanged.
    //  Enables correlating the local error against local mesh quality, which
    //  `checkMesh -writeAllFields` writes onto the same mesh.
    Switch writeErrorField_;

    label dimension_;
    bool errorsReported_;

    manufacturedEikonalVerifier
    (
        const manufacturedEikonalVerifier&
    ) = delete;
    void operator=(const manufacturedEikonalVerifier&) = delete;

public:

    OverrideTypeName("manufacturedEikonalVerifier");

    manufacturedEikonalVerifier
    (
        const dictionary& electroProperties,
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const dimensionedScalar& chi,
        const dimensionedScalar& Cm,
        const dimensionedScalar& c0,
        const Switch& eikonalAdvectionDiffusionApproach
    );

    static autoPtr<manufacturedEikonalVerifier> New
    (
        const dictionary& electroProperties,
        const fvMesh& mesh,
        const volTensorField& conductivity,
        const dimensionedScalar& chi,
        const dimensionedScalar& Cm,
        const dimensionedScalar& c0,
        const Switch& eikonalAdvectionDiffusionApproach
    );

    virtual bool enabled() const override
    {
        return true;
    }

    virtual void applyConstraints(volScalarField& activationTime) const override;

    virtual tmp<volScalarField> sourceTerm() const override;

    virtual bool shouldPostProcess() const override
    {
        return !errorsReported_;
    }

    virtual void postProcess(const volScalarField& activationTime) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
