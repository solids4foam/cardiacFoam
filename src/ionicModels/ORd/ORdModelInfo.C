/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "ORdModelInfo.H"
#include "ORd_2011.H"

const Foam::ionicModelFamilyInfo& Foam::ORdFamilyInfo()
{
    static const ionicModelFamilyInfo info
    {
        NUM_CONSTANTS,
        NUM_STATES,
        NUM_ALGEBRAIC,
        ORdCONSTANTS_NAMES,
        ORdSTATES_NAMES,
        ORdALGEBRAIC_NAMES,
        membrane_V,
        1000.0,
        1000.0,
        0.0,
        nullptr
    };

    return info;
}

// ************************************************************************* //
