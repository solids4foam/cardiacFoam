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
    manufacturedRadialEikonalVerifier

Description
    Eikonal verifier with an exact solution and no manufactured source: a held
    ball of radius R inside a sealed sphere of radius Ro = L + R, with a constant
    conductivity sigma Q diag(a^2) Q^T that stretches both into ellipsoids (Q is
    a rotation by angle about (1, 1, 1)). The first-arrival release stays on for
    the deposited junction seeds, so the case exercises the production path.

    Entries of the verificationModel dictionary: wallThickness (L), ballRatio
    (R/L), terminalTime, sigma, axes and angle. An eikonalPvjCoupler in
    domainCouplings seeds the ball at terminalTime with its seedProfile (flat,
    front or curvedFront), and needs a pvjRadius equal to R, no pvjDelay and
    axes (1 1 1). Without a coupler the verifier holds the ellipsoid with its
    own seedProfile: flat (at terminalTime) or exact (at the exact field,
    stopping at half its radius where the field is singular towards the
    centre).

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    manufacturedRadialEikonalVerifier.C
\*---------------------------------------------------------------------------*/

#ifndef manufacturedRadialEikonalVerifier_H
#define manufacturedRadialEikonalVerifier_H

#include "eikonalVerificationModel.H"
#include "overrideTypeName.H"
#include "tensor.H"

namespace Foam
{

class manufacturedRadialEikonalVerifier
:
    public eikonalVerificationModel
{
    const fvMesh& mesh_;
    const volTensorField& conductivity_;
    scalar L_;
    scalar a_;
    scalar terminalTime_;
    scalar sigma_;
    vector axes_;
    tensor Q_;
    bool seededByCoupler_;
    word seedProfile_;
    scalar chiCm_;
    scalar V_;
    scalar Pe_;
    bool errorsReported_;
    mutable bool configurationChecked_;

    //- Check the closed form against its self-test samples
    static void selfTest();

    //- The coordinate in which the solution is radial
    vector stretched(const point& x) const;

    //- The exact activation time for a stretched radius ratio s >= a/2 (s >= a
    //  for curvedFront, whose ball contains the ratio 2/Pe): flat and exact hold
    //  the ball surface at the terminal time, front at the terminal time plus
    //  R/V, and curvedFront follows the solution from the ratio 2/Pe
    scalar exactTime(const scalar s) const;

    //- The exact radial slope w at s >= a
    scalar exactSlope(const scalar s) const;

    //- The conductivity must be the tensor of the solution at every cell, and
    //  the mesh must lie inside the sealed wall
    void validateConfiguration() const;

public:

    OverrideTypeName("manufacturedRadialEikonalVerifier");

    manufacturedRadialEikonalVerifier
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

    //- Hold the ball when no coupler seeds it; deposited seeds are left to the release
    virtual void applyConstraints(volScalarField& activationTime) const override;

    //- No manufactured source: the exact field solves the production equation
    virtual tmp<volScalarField> sourceTerm() const override;

    virtual bool shouldPostProcess() const override
    {
        return !errorsReported_;
    }

    virtual bool imposesExactValue(const fvPatch&) const override
    {
        return false;
    }

    //- A ball held at the exact field stays held; the other seeds follow first arrival
    virtual bool releasesSeeds() const override
    {
        return seedProfile_ != "exact";
    }

    virtual void writeManufacturedFields(const fvMesh& mesh) const override;

    virtual void postProcess(const volScalarField& activationTime) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
