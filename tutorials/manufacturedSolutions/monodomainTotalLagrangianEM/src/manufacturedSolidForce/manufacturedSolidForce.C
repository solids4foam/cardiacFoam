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

#include "manufacturedSolidForce.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "IOdictionary.H"
#include "dimensionedScalar.H"

namespace Foam
{
namespace fv
{

defineTypeNameAndDebug(manufacturedSolidForce, 0);
addToRunTimeSelectionTable(option, manufacturedSolidForce, dictionary);

manufacturedSolidForce::manufacturedSolidForce
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    amplitude_(vector::zero),
    Tmax_(1.0),
    V0_(1.0),
    gamma_(1.0),
    E_(10000.0),
    nu_(0.3),
    rho0_(1060.0)
{
    read(dict);
}

manufacturedSolidForce::~manufacturedSolidForce()
{}

bool manufacturedSolidForce::read(const dictionary& dict)
{
    if (!option::read(dict))
    {
        return false;
    }

    // fv::option leaves fieldNames_/applied_ empty for options that derive
    // directly from the base (only cellSetOption populates them). Without this,
    // applyToField("D") returns -1 and addSup() is never invoked. Populate them
    // here so the body force is actually applied to the D equation.
    dict.readEntry("fieldNames", fieldNames_);
    applied_.setSize(fieldNames_.size(), false);

    // Active-tension parameters are deduced from the (top-level)
    // electroMechanicalProperties dictionary so that the manufactured body
    // force can never drift from the active-tension model the solver runs.
    IOdictionary emProps
    (
        IOobject
        (
            "electroMechanicalProperties",
            mesh_.time().constant(),
            mesh_.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word emModel(emProps.get<word>("electroMechanicalModel"));
    const dictionary& emCoeffs = emProps.subDict(emModel + "Coeffs");
    const dictionary& constants = emCoeffs.subDict("constants");

    Tmax_ = constants.get<scalar>("Tmax");
    V0_ = constants.get<scalar>("V0");
    gamma_ = constants.get<scalar>("gamma");

    // The coupler multiplies Ta by TaScale before it enters the solid stress,
    // so the manufactured body force must use the same effective amplitude.
    const scalar TaScale = emCoeffs.getOrDefault<scalar>("TaScale", 1e3);
    Tmax_ *= TaScale;

    // Manufactured displacement amplitude. Read from the same entry the
    // displacement BC uses so the two impose the same manufactured solution.
    amplitude_ =
        emCoeffs.subDict("electromechanicalVerificationModel")
       .get<vector>("amplitude");

    // Passive material parameters are deduced from the solid-region
    // mechanicalProperties passive law.
    IOdictionary mechProps
    (
        IOobject
        (
            "mechanicalProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const PtrList<entry> lawEntries(mechProps.lookup("mechanical"));
    const dictionary& passiveLaw =
        lawEntries[0].dict().subDict("passiveMechanicalLaw");

    rho0_ = dimensionedScalar("rho", dimDensity, lawEntries[0].dict()).value();
    E_ = dimensionedScalar("E", dimPressure, passiveLaw).value();
    nu_ = dimensionedScalar("nu", dimless, passiveLaw).value();

    Info<< "manufacturedSolidForce: deduced parameters" << nl
        << "    amplitude             = " << amplitude_ << nl
        << "    Tmax (TaScale-scaled) = " << Tmax_ << nl
        << "    V0                    = " << V0_ << nl
        << "    gamma                 = " << gamma_ << nl
        << "    rho0                  = " << rho0_ << nl
        << "    E                     = " << E_ << nl
        << "    nu                    = " << nu_ << endl;

    return true;
}

void manufacturedSolidForce::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const scalar t = mesh_.time().value();
    const volVectorField& C = mesh_.C();

    const scalar Ax = amplitude_.x();
    const scalar Ay = amplitude_.y();
    const scalar Az = amplitude_.z();

    const scalar mu = E_ / (2.0 * (1.0 + nu_));
    const scalar K = (nu_ * E_ / ((1.0 + nu_) * (1.0 - 2.0 * nu_))) + (2.0 / 3.0) * mu;

    const scalar rho0 = rho0_;
    const scalar Tmax = Tmax_;
    const scalar V0 = V0_;
    const scalar gamma = gamma_;

    const volVectorField* Dptr = mesh_.findObject<volVectorField>("D");

    forAll(eqn.source(), celli)
    {
        scalar X = C[celli].x();
        scalar Y = C[celli].y();
        scalar Z = C[celli].z();

        if (Dptr)
        {
            X -= (*Dptr)[celli].x();
            Y -= (*Dptr)[celli].y();
            Z -= (*Dptr)[celli].z();
        }

        #include "B_expr.H"

        eqn.source()[celli] += vector(Bx, By, Bz) * mesh_.V()[celli];
    }
}

}
}
