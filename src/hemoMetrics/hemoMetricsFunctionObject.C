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

#include "hemoMetricsFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;
using namespace Foam::functionObjects;

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(hemoMetricsFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, hemoMetricsFunctionObject, dictionary);
}
}

bool hemoMetricsFunctionObject::isWallPatch(const fvPatch& pp) const
{
    return isA<wallFvPatch>(pp);
}

hemoMetricsFunctionObject::hemoMetricsFunctionObject
(
    const word& name,
    const Time& rt,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, rt, dict),
    rho_(1060.0),
    mu_(3.5e-3),      // default 3.5 mPa·s
    UName_("U"),
    cyclePeriod_(1.0),
    startPhase_(0.0),
    cycleStartTime_(0.0)
{
    read(dict);
}

void hemoMetricsFunctionObject::ensureFields()
{
    const fvMesh& m = mesh_;

    if (!wss_.valid())
    {
        wss_.reset(new volVectorField
        (
            IOobject("wss", time().timeName(), m,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            m,
            dimensionedVector("zero", dimPressure, vector::zero)
        ));
    }

    if (!magWSS_.valid())
    {
        magWSS_.reset(new volScalarField
        (
            IOobject("magWSS", time().timeName(), m,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            m,
            dimensionedScalar("zero", dimPressure, 0.0)
        ));
    }

    if (!intMagWSS_.valid())
    {
        intMagWSS_.reset(new volScalarField
        (
            IOobject("magWSSInt", time().timeName(), m,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            m,
            dimensionedScalar("zero", dimPressure, 0.0)
        ));
    }

    if (!intWSS_.valid())
    {
        intWSS_.reset(new volVectorField
        (
            IOobject("wssInt", time().timeName(), m,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            m,
            dimensionedVector("zero", dimPressure, vector::zero)
        ));
    }

    if (!TAWSS_.valid())
    {
        TAWSS_.reset(new volScalarField
        (
            IOobject("TAWSS", time().timeName(), m,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            m,
            dimensionedScalar("zero", dimPressure, 0.0)
        ));
    }

    if (!OSI_.valid())
    {
        OSI_.reset(new volScalarField
        (
            IOobject("OSI", time().timeName(), m,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            m,
            dimensionedScalar("zero", dimless, 0.0)
        ));
    }

    if (!RRT_.valid())
    {
        RRT_.reset(new volScalarField
        (
            IOobject("RRT", time().timeName(), m,
                     IOobject::NO_READ, IOobject::NO_WRITE),
            m,
            dimensionedScalar("zero", dimless/dimPressure, 0.0) // corrected units
        ));
    }
}

void hemoMetricsFunctionObject::computeWSS()
{
    const fvMesh& m = mesh_;
    const volVectorField& U = m.lookupObject<volVectorField>(UName_);

    // Symmetric rate-of-strain tensor S = 0.5*(gradU + gradU^T)
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU();

    volSymmTensorField S
    (
        IOobject("S_tmp", time().timeName(), m, IOobject::NO_READ, IOobject::NO_WRITE),
        symm(gradU)
    );

    wss_->primitiveFieldRef() = vector::zero;
    magWSS_->primitiveFieldRef() = 0.0;

    const fvBoundaryMesh& bnd = m.boundary();

    forAll(bnd, patchI)
    {
        const fvPatch& fvp = bnd[patchI];
        if (!isWallPatch(fvp)) continue;
        if (patchNames_.size() && !patchNames_.found(fvp.name())) continue;

        const vectorField& Sf = m.Sf().boundaryField()[patchI];
        const scalarField& magSf = m.magSf().boundaryField()[patchI];
        tmp<vectorField> tn = -Sf/magSf;
        const vectorField& n = tn();  // bind to underlying Field

        const symmTensorField& Sp = S.boundaryField()[patchI];

        vectorField traction(Sp.size(), vector::zero);
        forAll(Sp, i)
        {
            vector Sn = vector(Sp[i] & n[i]);      // S·n
            traction[i] = 2.0*mu_*Sn;              // 2*mu*S·n
        }

        vectorField tauw(traction.size(), vector::zero);
        forAll(traction, i)
        {
            const vector& ni = n[i];
            const vector& ti = traction[i];
            tauw[i] = ti - ni*(ti & ni);           // project tangentially
        }

        wss_->boundaryFieldRef()[patchI] = tauw;
        magWSS_->boundaryFieldRef()[patchI] = mag(tauw);
    }
}

void hemoMetricsFunctionObject::accumulate()
{
    const scalar dt = time().deltaTValue();

    intMagWSS_->primitiveFieldRef() += magWSS_->primitiveField()*dt;
    intWSS_->primitiveFieldRef()    += wss_->primitiveField()*dt;

    forAll(mesh_.boundary(), patchI)
    {
        if (!isWallPatch(mesh_.boundary()[patchI])) continue;

        intMagWSS_->boundaryFieldRef()[patchI] += magWSS_->boundaryField()[patchI]*dt;
        intWSS_->boundaryFieldRef()[patchI]    += wss_->boundaryField()[patchI]*dt;
    }
}

void hemoMetricsFunctionObject::finalizeCycleAndWrite()
{
    const scalar T = max(SMALL, cyclePeriod_);

    // TAWSS = (1/T) ∫ |WSS| dt
    TAWSS_->primitiveFieldRef() = intMagWSS_->primitiveField()/T;

    // OSI = 0.5 * (1 - |∫ WSS dt| / ∫ |WSS| dt)
    OSI_->primitiveFieldRef() = 0.0;
    forAll(OSI_->primitiveField(), i)
    {
        const scalar denom = max(SMALL, intMagWSS_->primitiveField()[i]);
        const scalar num   = mag(intWSS_->primitiveField()[i]);
        OSI_->primitiveFieldRef()[i] = 0.5*(1.0 - num/denom);
    }

    // RRT = 1 / ( (1 - 2*OSI) * TAWSS )
    RRT_->primitiveFieldRef() = 0.0;
    forAll(RRT_->primitiveField(), i)
    {
        const scalar g = (1.0 - 2.0*OSI_->primitiveField()[i]) * TAWSS_->primitiveField()[i];
        RRT_->primitiveFieldRef()[i] = (mag(g) > SMALL ? 1.0/g : 0.0);
    }

    // boundary fields
    forAll(mesh_.boundary(), patchI)
    {
        if (!isWallPatch(mesh_.boundary()[patchI])) continue;

        const scalarField& intMag = intMagWSS_->boundaryField()[patchI];
        const vectorField& intVec = intWSS_->boundaryField()[patchI];

        scalarField TAW(intMag.size(), 0.0);
        scalarField osi(intMag.size(), 0.0);
        scalarField rrt(intMag.size(), 0.0);

        forAll(TAW, i)
        {
            TAW[i] = intMag[i]/T;
            const scalar denom = max(SMALL, intMag[i]);
            osi[i] = 0.5*(1.0 - mag(intVec[i])/denom);
            const scalar g = (1.0 - 2.0*osi[i]) * TAW[i];
            rrt[i] = (mag(g) > SMALL ? 1.0/g : 0.0);
        }

        TAWSS_->boundaryFieldRef()[patchI] = TAW;
        OSI_->boundaryFieldRef()[patchI]   = osi;
        RRT_->boundaryFieldRef()[patchI]   = rrt;
    }

    // Write only at solver write steps
    if (time().outputTime())
    {
        TAWSS_->write();
        OSI_->write();
        RRT_->write();
    }

    // reset integrals
    intMagWSS_->primitiveFieldRef() = 0.0;
    intWSS_->primitiveFieldRef()    = vector::zero;
    forAll(mesh_.boundary(), patchI)
    {
        if (!isWallPatch(mesh_.boundary()[patchI])) continue;
        intMagWSS_->boundaryFieldRef()[patchI] = 0.0;
        intWSS_->boundaryFieldRef()[patchI]    = vector::zero;
    }
}

bool hemoMetricsFunctionObject::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("patches",      patchNames_);
    dict.readIfPresent("rho",          rho_);
    dict.readIfPresent("mu",           mu_);
    dict.readIfPresent("U",            UName_);
    dict.readIfPresent("cyclePeriod",  cyclePeriod_);
    dict.readIfPresent("startPhase",   startPhase_);

    return true;
}

void hemoMetricsFunctionObject::start()
{
    ensureFields();
    cycleStartTime_ = time().value() + startPhase_;
}

bool hemoMetricsFunctionObject::execute()
{
    ensureFields();
    computeWSS();
    accumulate();

    const scalar t = time().value();

    // Only finalize exactly at cycle end (within tolerance)
    if (mag(t - (cycleStartTime_ + cyclePeriod_)) < SMALL)
    {
        finalizeCycleAndWrite();
        cycleStartTime_ += cyclePeriod_;
    }

    return true;
}

bool hemoMetricsFunctionObject::end()
{
    finalizeCycleAndWrite();
    return true;
}

bool hemoMetricsFunctionObject::write()
{
    if (time().outputTime())   // respect solver write interval
    {
        wss_->write();
        magWSS_->write();
    }
    return true;
}
