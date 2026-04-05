/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "BuenoOrovioModelInfo.H"
#include "BuenoOrovio_2008.H"

namespace
{
    Foam::scalar buenoOrovioTransformedVm(const Foam::scalarField& S)
    {
        return S[0] * 85.7 - 84.0;
    }
}

const Foam::ionicModelFamilyInfo& Foam::BuenoOrovioFamilyInfo()
{
    static const ionicModelFamilyInfo info
    {
        NUM_CONSTANTS,
        NUM_STATES,
        NUM_ALGEBRAIC,
        BuenoOrovioCONSTANTS_NAMES,
        BuenoOrovioSTATES_NAMES,
        BuenoOrovioALGEBRAIC_NAMES,
        u,
        1000.0,
        1000.0/85.7,
        84.0/85.7,
        &buenoOrovioTransformedVm
    };

    return info;
}

// ************************************************************************* //
