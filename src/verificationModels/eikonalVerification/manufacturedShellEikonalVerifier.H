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
    manufacturedShellEikonalVerifier

Description
    Manufactured eikonal verifier on the spherical shell r1 < |x| < r2, sealed
    on both spheres: activation time psi with n.M grad(psi) = 0 on each sphere,
    M = G/(chi Cm) with G full and rotating through the wall (written cell-wise
    as 0/Conductivity by setManufacturedFields), and the manufactured source
    that makes psi the solution of
    -div(M grad psi) + c0 sqrt(grad psi . M grad psi) = 1 + S.
    The sealed problem fixes psi up to a constant, so the cells whose centre
    lies within rPin of the stimulus point are pinned to the exact psi (the
    nearest cell when none lies within rPin); no patch imposes an exact value.
    With initialiseExact the nonlinear iteration starts from the exact field,
    and pinnedRadius replaces rPin. Fields from
    manufacturedShellEikonalReference.H.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    manufacturedShellEikonalVerifier.C
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellEikonalVerifier_H
#define manufacturedShellEikonalVerifier_H

#include "eikonalVerificationModel.H"
#include "overrideTypeName.H"
#include "Switch.H"
#include "eikonalVerification/manufacturedShellEikonalReference.H"

namespace Foam
{

class manufacturedShellEikonalVerifier
:
    public eikonalVerificationModel
{
    const fvMesh& mesh_;
    const volTensorField& conductivity_;
    shellEikonalMMS::Params p_;
    label nPerDirection_;
    Switch initialiseExact_;
    Switch advectionDiffusion_;
    mutable volScalarField* activationTimePtr_;
    mutable label nPinned_;
    bool errorsReported_;

    //- Check the header functions against their self-test samples
    static void selfTest();

    void validateConfiguration
    (
        const dimensionedScalar& chi,
        const dimensionedScalar& Cm
    ) const;

    //- The conductivity must be the manufactured G at every cell centre
    void validateConductivity() const;

    //- Set a cell to the exact activation time, which must not be negative
    void pinCell
    (
        scalarField& values,
        const label cellI,
        const point& centre
    ) const;

public:

    OverrideTypeName("manufacturedShellEikonalVerifier");

    manufacturedShellEikonalVerifier
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

    //- Pin the cells within rPin of the stimulus point to the exact field
    virtual void applyConstraints(volScalarField& activationTime) const override;

    virtual tmp<volScalarField> sourceTerm() const override;

    virtual bool shouldPostProcess() const override
    {
        return !errorsReported_;
    }

    virtual bool imposesExactValue(const fvPatch&) const override
    {
        return false;
    }

    virtual void writeManufacturedFields(const fvMesh& mesh) const override;

    virtual void postProcess(const volScalarField& activationTime) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
