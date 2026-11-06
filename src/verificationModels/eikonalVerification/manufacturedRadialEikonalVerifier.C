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

#include "manufacturedRadialEikonalVerifier.H"
#include "manufacturedRadialEikonalReference.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "PstreamReduceOps.H"
#include "verificationUtils.H"
#include "dictionaryKeyCheck.H"
#include "eikonalPvjCoupler.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedRadialEikonalVerifier, 0);
addToRunTimeSelectionTable
(
    eikonalVerificationModel,
    manufacturedRadialEikonalVerifier,
    dictionary
);

namespace
{

const dimensionSet conductivityDim(-1, -3, 3, 0, 0, 2, 0);

//- Rotation by angle about (1, 1, 1)
tensor rotationAboutDiagonal(const scalar angle)
{
    const scalar k = 1/Foam::sqrt(scalar(3));
    const tensor K(0, -k, k, k, 0, -k, -k, k, 0);
    return tensor(I) + Foam::sin(angle)*K + (1 - Foam::cos(angle))*(K & K);
}


//- The eikonalPvjCoupler block of domainCouplings, or nullptr when there is none
const dictionary* seedingCoupler(const dictionary& coeffs)
{
    const dictionary* coupler = nullptr;
    if (const dictionary* couplings = coeffs.findDict("domainCouplings"))
    {
        for (const entry& e : *couplings)
        {
            if
            (
                e.isDict()
             && e.dict().lookupOrDefault<word>("electroDomainCoupler", word::null)
             == eikonalPvjCoupler::typeName
            )
            {
                if (coupler)
                {
                    FatalIOErrorInFunction(*couplings)
                        << "The radial verifier takes one eikonalPvjCoupler."
                        << exit(FatalIOError);
                }
                coupler = &e.dict();
            }
        }
    }
    return coupler;
}


//- The seed profile of the ball: the coupler's when one seeds the ball,
//  otherwise the verifier's own, flat or exact, with which it holds the ball
word readSeedProfile(const dictionary& coeffs, const scalar R, const vector& axes)
{
    const dictionary& verifierDict = coeffs.subDict("verificationModel");
    const dictionary* coupler = seedingCoupler(coeffs);
    if (!coupler)
    {
        const word profile(verifierDict.get<word>("seedProfile"));
        dictionaryKeyCheck::checkValue
        (
            verifierDict, "seedProfile", profile, wordList({"flat", "exact"})
        );
        return profile;
    }

    if (verifierDict.found("seedProfile"))
    {
        FatalIOErrorInFunction(verifierDict)
            << "seedProfile is set by the eikonalPvjCoupler that seeds the ball."
            << exit(FatalIOError);
    }
    if (mag(axes - vector::one) > SMALL)
    {
        FatalIOErrorInFunction(verifierDict)
            << "A ball seeded by the eikonalPvjCoupler needs axes (1 1 1): "
            << "the coupler seeds a sphere." << exit(FatalIOError);
    }
    const scalar pvjRadius = coupler->get<scalar>("pvjRadius");
    if (mag(pvjRadius - R) > 1e-9*R)
    {
        FatalIOErrorInFunction(*coupler)
            << "pvjRadius " << pvjRadius << " differs from the ball radius "
            << "ballRatio*wallThickness = " << R << "." << exit(FatalIOError);
    }
    if (coupler->lookupOrDefault<scalar>("pvjDelay", 0) != 0)
    {
        FatalIOErrorInFunction(*coupler)
            << "The ball is seeded at terminalTime: pvjDelay must be 0."
            << exit(FatalIOError);
    }
    return eikonalPvjCoupler::seedProfile(*coupler);
}

} // End anonymous namespace


void manufacturedRadialEikonalVerifier::selfTest()
{
    for (const radialEikonalMMS::SelfTestSample& t : radialEikonalMMS::selfTestSamples)
    {
        if
        (
            mag(radialEikonalMMS::slope(t.s, t.Pe, t.a) - t.slope) > 1e-11
         || mag(radialEikonalMMS::scaledTime(t.s, t.Pe, t.a) - t.scaledTime) > 1e-10
        )
        {
            FatalErrorInFunction
                << "Radial eikonal self-test failed at Pe = " << t.Pe
                << ", a = " << t.a << ", s = " << t.s << exit(FatalError);
        }
    }
}


manufacturedRadialEikonalVerifier::manufacturedRadialEikonalVerifier
(
    const dictionary& electroProperties,
    const fvMesh& mesh,
    const volTensorField& conductivity,
    const dimensionedScalar& chi,
    const dimensionedScalar& Cm,
    const dimensionedScalar& c0,
    const Switch& eikonalAdvectionDiffusionApproach
)
:
    eikonalVerificationModel(electroProperties),
    mesh_(mesh),
    conductivity_(conductivity),
    L_(electroProperties.subDict("verificationModel").get<scalar>("wallThickness")),
    a_(electroProperties.subDict("verificationModel").get<scalar>("ballRatio")),
    terminalTime_(electroProperties.subDict("verificationModel").get<scalar>("terminalTime")),
    sigma_(electroProperties.subDict("verificationModel").get<scalar>("sigma")),
    axes_(electroProperties.subDict("verificationModel").get<vector>("axes")),
    Q_(rotationAboutDiagonal(electroProperties.subDict("verificationModel").get<scalar>("angle"))),
    seededByCoupler_(seedingCoupler(electroProperties) != nullptr),
    seedProfile_(readSeedProfile(electroProperties, a_*L_, axes_)),
    chiCm_(chi.value()*Cm.value()),
    V_(c0.value()*Foam::sqrt(sigma_/chiCm_)),
    Pe_(c0.value()*L_/Foam::sqrt(sigma_/chiCm_)),
    errorsReported_(false),
    configurationChecked_(false)
{
    selfTest();

    if (!eikonalAdvectionDiffusionApproach)
    {
        FatalIOErrorInFunction(electroProperties)
            << typeName << " needs eikonalAdvectionDiffusionApproach true."
            << exit(FatalIOError);
    }
    if
    (
        L_ <= 0 || a_ <= 0 || sigma_ <= 0
     || axes_.x() <= 0 || axes_.y() <= 0 || axes_.z() <= 0 || Pe_ <= 0
    )
    {
        FatalIOErrorInFunction(electroProperties.subDict("verificationModel"))
            << "wallThickness, ballRatio, sigma, axes and c0 must be positive."
            << exit(FatalIOError);
    }
    if (seedProfile_ == "curvedFront" && 2/Pe_ > a_)
    {
        FatalIOErrorInFunction(electroProperties.subDict("verificationModel"))
            << "seedProfile curvedFront needs ballRatio at least 2/Pe = " << 2/Pe_
            << ": the exact field holds outside the critical radius."
            << exit(FatalIOError);
    }
}


vector manufacturedRadialEikonalVerifier::stretched(const point& x) const
{
    const vector local(Q_.T() & x);
    return vector(local.x()/axes_.x(), local.y()/axes_.y(), local.z()/axes_.z());
}


scalar manufacturedRadialEikonalVerifier::exactTime(const scalar s) const
{
    if (seedProfile_ == "curvedFront")
    {
        return
            terminalTime_
          + (L_/V_)
           *(
                radialEikonalMMS::scaledTime(s, Pe_, a_)
              - radialEikonalMMS::scaledTime(2/Pe_, Pe_, a_)
            );
    }
    const scalar edge = seedProfile_ == "front" ? a_ : scalar(0);
    return terminalTime_ + (L_/V_)*(edge + radialEikonalMMS::scaledTime(s, Pe_, a_));
}


scalar manufacturedRadialEikonalVerifier::exactSlope(const scalar s) const
{
    return radialEikonalMMS::slope(s, Pe_, a_);
}


void manufacturedRadialEikonalVerifier::validateConfiguration() const
{
    const tensor expected
    (
        sigma_*(Q_ & diagTensor(sqr(axes_.x()), sqr(axes_.y()), sqr(axes_.z())) & Q_.T())
    );
    scalar worst = 0;
    forAll(conductivity_, cellI)
    {
        worst = max(worst, mag(conductivity_[cellI] - expected)/mag(expected));
    }
    reduce(worst, maxOp<scalar>());
    if (worst > 1e-9)
    {
        FatalErrorInFunction
            << typeName << " needs the conductivity written by "
            << "setManufacturedFields; the largest relative difference is "
            << worst << "." << exit(FatalError);
    }

    label nOutside = 0;
    for (const point& pt : mesh_.points())
    {
        if (mag(stretched(pt)) > (1 + a_)*L_*(1 + 1e-8))
        {
            ++nOutside;
        }
    }
    reduce(nOutside, sumOp<label>());
    if (nOutside)
    {
        FatalErrorInFunction
            << typeName << ": " << nOutside
            << " mesh points lie outside the sealed wall." << exit(FatalError);
    }
}


void manufacturedRadialEikonalVerifier::applyConstraints
(
    volScalarField& activationTime
) const
{
    if (!configurationChecked_)
    {
        validateConfiguration();
        configurationChecked_ = true;
    }

    if (seededByCoupler_)
    {
        return;
    }

    const vectorField& C = mesh_.C().primitiveField();
    scalarField& values = activationTime.primitiveFieldRef();
    forAll(C, cellI)
    {
        const scalar s = mag(stretched(C[cellI]))/L_;
        if (s <= a_)
        {
            values[cellI] =
                seedProfile_ == "exact" ? exactTime(max(s, 0.5*a_)) : terminalTime_;
            if (values[cellI] < 0)
            {
                FatalErrorInFunction
                    << "terminalTime is too small: a held value is negative."
                    << exit(FatalError);
            }
        }
    }
}


tmp<volScalarField> manufacturedRadialEikonalVerifier::sourceTerm() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Smms",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("Smms", dimless, 0.0)
        )
    );
}


void manufacturedRadialEikonalVerifier::postProcess
(
    const volScalarField& activationTime
)
{
    if (!shouldPostProcess())
    {
        return;
    }

    const vectorField& C = mesh_.C().primitiveField();
    const scalarField& V = mesh_.V();
    const scalarField& T = activationTime.primitiveField();
    const scalar span = (L_/V_)*radialEikonalMMS::scaledTime(a_ + 1, Pe_, a_);

    const auto exactAt = [&](const point& x)
    {
        return exactTime(max(mag(stretched(x))/L_, a_));
    };

    // The gradient of the exact field in x from its radial form in the stretched coordinate
    vectorField gradExact(C.size(), Zero);
    boolList isFree(C.size(), false);
    scalarField error(C.size(), 0);
    forAll(C, cellI)
    {
        const vector xi(stretched(C[cellI]));
        const scalar s = mag(xi)/L_;
        isFree[cellI] = s > a_;
        if (isFree[cellI])
        {
            error[cellI] = (T[cellI] - exactTime(s))/span;
            const vector dXi((exactSlope(s)/V_)*xi/mag(xi));
            gradExact[cellI] =
                Q_ & vector(dXi.x()/axes_.x(), dXi.y()/axes_.y(), dXi.z()/axes_.z());
        }
    }

    const boolList isWall(boundaryAdjacentCells(mesh_));
    const boolList besideBall(interfaceBandCells(mesh_, isFree));

    const tmp<volVectorField> tgrad(fvc::grad(activationTime));
    const tmp<volVectorField> tgradRecon
    (
        exactFieldGradient(mesh_, activationTime.name(), exactAt)
    );

    scalar gradScale = 0;
    forAll(C, cellI)
    {
        gradScale += V[cellI]*magSqr(gradExact[cellI]);
    }
    const scalar total = gSum(V);
    gradScale = Foam::sqrt(returnReduce(gradScale, sumOp<scalar>())/total);

    //- L1, L2 and Linf over the selected cells; every part uses the whole volume
    const auto normsOf = [&](const auto& select, const auto& value)
    {
        scalar l1 = 0, l2 = 0, lInf = 0;
        forAll(C, cellI)
        {
            if (select(cellI))
            {
                const scalar e = value(cellI);
                l1 += V[cellI]*e;
                l2 += V[cellI]*e*e;
                lInf = max(lInf, e);
            }
        }
        reduce(l1, sumOp<scalar>());
        reduce(l2, sumOp<scalar>());
        reduce(lInf, maxOp<scalar>());
        return Tuple2<Tuple2<scalar, scalar>, scalar>
        (
            Tuple2<scalar, scalar>(l1/total, Foam::sqrt(l2/total)),
            lInf
        );
    };

    const auto errorAt = [&](const label c) { return mag(error[c]); };
    const auto isBulk = [&](const label c)
    {
        return isFree[c] && !isWall[c] && !besideBall[c];
    };
    const auto gradAt = [&](const label c)
    {
        return mag(tgrad().primitiveField()[c] - gradExact[c])/gradScale;
    };
    const auto gradReconAt = [&](const label c)
    {
        return mag(tgradRecon().primitiveField()[c] - gradExact[c])/gradScale;
    };

    const auto eps = normsOf([&](const label c) { return isFree[c]; }, errorAt);
    const auto epsBulk = normsOf(isBulk, errorAt);
    const auto epsWall = normsOf([&](const label c) { return isFree[c] && isWall[c]; }, errorAt);
    const auto epsBall = normsOf([&](const label c) { return isFree[c] && besideBall[c]; }, errorAt);
    const auto gradBulk = normsOf(isBulk, gradAt);
    const auto gradReconBulk = normsOf(isBulk, gradReconAt);

    const label nCells = returnReduce(mesh_.nCells(), sumOp<label>());
    const scalar h = Foam::cbrt(total/max(nCells, label(1)));
    const fileName outputDir(mesh_.time().globalPath()/"postProcessing");

    if (Pstream::master())
    {
        mkDir(outputDir);
        OFstream os(outputDir/"radial_errors.dat");
        os.precision(12);
        os  << "# Radial eikonal errors, relative to the activation span and the RMS exact gradient\n"
            << "# Every row uses the whole-domain volume, so the parts add in quadrature\n"
            << "# PeL " << Pe_ << "\n"
            << "# ballRatio " << a_ << "\n"
            << "# seedProfile " << seedProfile_ << "\n"
            << "# hoverL " << h/L_ << "\n"
            << "# cells " << nCells << "\n"
            << "# span " << span << "\n"
            << "# field L1 L2 Linf\n";
        writeNormRow(os, "eps", eps);
        writeNormRow(os, "eps_bulk", epsBulk);
        writeNormRow(os, "eps_wall", epsWall);
        writeNormRow(os, "eps_ball", epsBall);
        writeNormRow(os, "grad_bulk", gradBulk);
        writeNormRow(os, "gradRecon_bulk", gradReconBulk);
    }
    Info<< nl << "Radial eikonal error summary" << nl;
    writeNormRow(Info, "eps", eps);
    writeNormRow(Info, "grad_bulk", gradBulk);

    errorsReported_ = true;
}


void manufacturedRadialEikonalVerifier::writeManufacturedFields
(
    const fvMesh& mesh
) const
{
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << typeName << " writes the manufactured conductivity of the "
            << "undecomposed case; run it in serial." << exit(FatalError);
    }

    const tensor D
    (
        sqr(axes_.x()), 0, 0,
        0, sqr(axes_.y()), 0,
        0, 0, sqr(axes_.z())
    );
    const unsigned int oldPrecision = IOstream::defaultPrecision(17);

    volSymmTensorField G
    (
        IOobject("Conductivity", "0", mesh, IOobject::NO_READ, IOobject::NO_WRITE),
        mesh,
        dimensionedSymmTensor(conductivityDim, symm(sigma_*(Q_ & D & Q_.T()))),
        "zeroGradient"
    );
    G.write();

    IOstream::defaultPrecision(oldPrecision);

    Info<< typeName << ": wrote Conductivity in 0." << endl;
}

} // End namespace Foam

// ************************************************************************* //
