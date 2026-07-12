/*---------------------------------------------------------------------------*\
Description
    Coefficient-level tests for extracellular face interpolation.
\*---------------------------------------------------------------------------*/

#include "extracellularFaceConductivity.H"
#include "IOstreams.H"

using namespace Foam;
using namespace Foam::extracellularFaceConductivity;

namespace
{

bool close(const scalar actual, const scalar expected)
{
    return mag(actual - expected) <= 1e-12*max(scalar(1), mag(expected));
}


void require(const bool condition, const char* message)
{
    if (!condition)
    {
        FatalErrorInFunction << message << exit(FatalError);
    }
}

}


int main()
{
    const scalar sigmaE = 0.04;
    const scalar sigmaB = 0.02;

    const scalar equalUnweighted = unweightedHarmonic(sigmaE, sigmaB);
    const scalar equalWeighted =
        distanceWeightedHarmonic(sigmaE, sigmaB, 0.5, 0.5);

    require
    (
        close(equalUnweighted, equalWeighted),
        "Equal-distance harmonic coefficients must agree."
    );

    const scalar distanceHeart = 0.25;
    const scalar distanceBath = 0.75;
    const scalar exactResistance =
        distanceHeart/sigmaE + distanceBath/sigmaB;
    const scalar exactCoefficient =
        (distanceHeart + distanceBath)/exactResistance;
    const scalar weighted = distanceWeightedHarmonic
    (
        sigmaE,
        sigmaB,
        distanceHeart,
        distanceBath
    );

    require
    (
        close(weighted, exactCoefficient),
        "Distance-weighted harmonic coefficient must reproduce exact series resistance."
    );
    require
    (
        !close(equalUnweighted, exactCoefficient),
        "Unweighted harmonic coefficient must expose the unequal-distance mismatch."
    );

    const scalar sigmaI = 0.12;
    const tensor heartSigmaTotal((sigmaI + sigmaE)*tensor::I);
    const tensor heartSigmaI(sigmaI*tensor::I);
    const tensor bathSigma(sigmaB*tensor::I);
    const tensor zeroTensor(tensor::zero);
    const tensor physicsCorrect = unweightedExtracellularFaceTensor
    (
        heartSigmaTotal,
        bathSigma,
        heartSigmaI,
        zeroTensor
    );
    const scalar naiveLinear = 0.5*(sigmaI + sigmaE) + 0.5*sigmaB;
    const tensor weightedPhysicsCorrect =
        distanceWeightedExtracellularFaceTensor
        (
            heartSigmaTotal,
            bathSigma,
            heartSigmaI,
            zeroTensor,
            distanceBath/(distanceHeart + distanceBath)
        );

    require
    (
        close(physicsCorrect.xx(), equalUnweighted),
        "Heart-bath coefficient must exclude intracellular conductivity."
    );
    require
    (
        !close(physicsCorrect.xx(), naiveLinear),
        "Naive sigmaTotal interpolation must differ when Gi contaminates the interface."
    );
    require
    (
        close(weightedPhysicsCorrect.xx(), exactCoefficient),
        "Weighted heart-bath tensor coefficient must reproduce series resistance."
    );

    Info<< "equalDistanceCoefficient=" << equalWeighted << nl
        << "unequalDistanceExactCoefficient=" << exactCoefficient << nl
        << "unequalDistanceUnweightedCoefficient=" << equalUnweighted << nl
        << "physicsCorrectInterfaceCoefficient=" << physicsCorrect.xx() << nl
        << "naiveLinearSigmaTotalCoefficient=" << naiveLinear << nl
        << "PASS" << endl;

    return 0;
}

// ************************************************************************* //
