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

#ifndef BiDomainDiffusionSolver_H
#define BiDomainDiffusionSolver_H

#include "reactionDiffusionSolver.H"
#include "fvc.H"
#include "fvm.H"

namespace Foam
{

class BiDomainDiffusionSolver
:
    public reactionDiffusionSolver
{
    volScalarField& phiE_;
    volTensorField Gi_;
    volTensorField Ge_;
    volTensorField GiPlusGe_;
    label phiEReferenceCell_;
    scalar phiEReferenceValue_;

    tmp<volTensorField> initialiseConductivityTensor
    (
        const fvMesh& mesh,
        const word& fieldName,
        const dictionary& dict
    ) const;

public:

    BiDomainDiffusionSolver
    (
        const fvMesh& mesh,
        const dictionary& electroProperties,
        volScalarField& phiE
    );

    virtual const volScalarField* phiEPtr() const override
    {
        return &phiE_;
    }

    virtual const volTensorField* conductivityPtr() const override
    {
        return &Gi_;
    }

    virtual const volTensorField& explicitConductivityField() const override
    {
        return Gi_;
    }

    virtual void solveDiffusionExplicit
    (
        MyocardiumDomain& domain,
        scalar dt
    ) override;

    virtual void solveDiffusionImplicit
    (
        MyocardiumDomain& domain,
        scalar dt,
        pimpleControl& pimple
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
