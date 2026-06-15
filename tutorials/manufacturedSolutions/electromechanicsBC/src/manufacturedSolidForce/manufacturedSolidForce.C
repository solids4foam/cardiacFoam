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
    nu_(0.3)
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

    E_ = dimensionedScalar("E", dimPressure, passiveLaw).value();
    nu_ = dimensionedScalar("nu", dimless, passiveLaw).value();

    Info<< "manufacturedSolidForce: deduced parameters" << nl
        << "    amplitude             = " << amplitude_ << nl
        << "    Tmax (TaScale-scaled) = " << Tmax_ << nl
        << "    V0                    = " << V0_ << nl
        << "    gamma                 = " << gamma_ << nl
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

    const scalar Tmax = Tmax_;
    const scalar V0 = V0_;
    const scalar gamma = gamma_;
    const scalar pi = constant::mathematical::pi;

    forAll(eqn.source(), celli)
    {
        const scalar X = C[celli].x();
        const scalar Y = C[celli].y();
        const scalar Z = C[celli].z();

        #include "B_expr.H"

        eqn.source()[celli] += vector(Bx, By, Bz) * mesh_.V()[celli];
    }
}

}
}
