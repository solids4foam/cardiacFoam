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

#include "Switch.H"
#include "autoPtr.H"
#include "dimensionedScalar.H"
#include "dictionary.H"
#include "volFields.H"

namespace Foam
{

class manufacturedEikonalVerifier
{
    const fvMesh& mesh_;
    const volTensorField& conductivity_;
    scalar chi_;
    scalar Cm_;
    scalar c0_;
    Switch enabled_;
    label dimension_;
    bool errorsReported_;

    manufacturedEikonalVerifier
    (
        const manufacturedEikonalVerifier&
    ) = delete;
    void operator=(const manufacturedEikonalVerifier&) = delete;

public:

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

    bool enabled() const
    {
        return enabled_;
    }

    void applyConstraints(volScalarField& activationTime) const;

    bool shouldPostProcess() const
    {
        return enabled_ && !errorsReported_;
    }

    void postProcess(const volScalarField& activationTime);
};

} // End namespace Foam

#endif

// ************************************************************************* //
