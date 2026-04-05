/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.
\*---------------------------------------------------------------------------*/

#include "TNNPModelInfo.H"
#include "TNNP_2004.H"

const Foam::ionicModelFamilyInfo& Foam::TNNPFamilyInfo()
{
    static const ionicModelFamilyInfo info
    {
        NUM_CONSTANTS,
        NUM_STATES,
        NUM_ALGEBRAIC,
        TNNP_CONSTANTS_NAMES,
        TNNP_STATES_NAMES,
        TNNP_ALGEBRAIC_NAMES,
        V,
        1000.0,
        1000.0,
        0.0,
        nullptr
    };

    return info;
}

// ************************************************************************* //
