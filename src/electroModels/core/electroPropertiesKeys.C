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

#include "electroPropertiesKeys.H"
#include "dictionaryKeyCheck.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

using namespace Foam;

// ODE integrator settings, flat in the block that selects the ionic model
const wordList odeKeys
({
    "solver", "absTol", "relTol", "maxSteps",
    "safeScale", "alphaInc", "alphaDec", "minScale", "maxScale",
    "ODESolver"
});

// Read by ionicModel from the block that selects it
const wordList ionicKeys
({
    "ionicModel", "tissue", "sex", "dimension", "ionicHeterogeneity",
    "ionicConstantOverrides", "outputVariables", "utilities",
    "initSampleCell", "batchedIntegrator", "batchedSubsteps",
    "batchedParallelCells", "batchedParallelMinCells",
    "storeBatchedAlgebraics"
});

// Myocardium tissue solvers (monodomainSolver, bidomainSolver, eikonalSolver)
const wordList tissueKeys
({
    "chi", "cm", "conductivity", "conductivitySource",
    "conductivityIntracellular", "conductivityExtracellular", "cellZone",
    "sealedHeartBoundary", "exposedFacesPatch", "solutionAlgorithm",
    "timeCouplingScheme", "activationThreshold", "externalStimulus",
    "singleCellStimulus", "verificationModel", "electrophysicsAdvanceScheme",
    "bathPredictorCorrector", "bathPotentialDomain", "bathHeartPhiETrace",
    "phiERefPoint", "phiEReferenceValue", "ecgDomains",
    "conductionNetworkDomains", "domainCouplings", "c0",
    "eikonalAdvectionDiffusionApproach", "stimulusLocationMin",
    "stimulusLocationMax", "useGraphPrePopulation"
});

// singleCellSolver, including its optional active-tension model
const wordList singleCellKeys
({
    "singleCellStimulus", "writeAfterTime", "writeFrequency",
    "verificationModel", "activeTensionModel", "constants", "initialStates",
    "couplingSignal", "preconditioningTime"
});

const wordList singleCellStimulusKeys
({
    "stim_start", "stim_period_S1", "stim_duration", "stim_amplitude",
    "nstim1", "stim_period_S2", "nstim2"
});


wordRes join(std::initializer_list<const wordList*> lists)
{
    label n = 0;
    for (const wordList* l : lists)
    {
        n += l->size();
    }

    wordRes result(n);
    n = 0;
    for (const wordList* l : lists)
    {
        for (const word& w : *l)
        {
            result[n++] = wordRe(w);
        }
    }

    return result;
}


wordRes keys(const wordList& names)
{
    return join({&names});
}


void checkSubDict
(
    const dictionary& dict,
    const word& name,
    const wordRes& allowed
)
{
    if (const dictionary* subPtr = dict.findDict(name))
    {
        dictionaryKeyCheck::checkKeys(*subPtr, allowed);
    }
}


void checkOutputVariables(const dictionary& coeffs)
{
    if (const dictionary* ovPtr = coeffs.findDict("outputVariables"))
    {
        dictionaryKeyCheck::checkKeys
        (
            *ovPtr,
            keys({"ionic", "activeTension"})
        );
        checkSubDict(*ovPtr, "ionic", keys({"export", "debug"}));
        checkSubDict(*ovPtr, "activeTension", keys({"export", "debug"}));
    }
}


void checkIonicHeterogeneity(const dictionary& coeffs)
{
    const dictionary* hetPtr = coeffs.findDict("ionicHeterogeneity");

    if (!hetPtr)
    {
        return;
    }

    dictionaryKeyCheck::checkKeys
    (
        *hetPtr,
        keys
        ({
            "mode", "field", "regions", "gradientAxes", "smoothing",
            "transitionMode", "transitionWidth"
        })
    );

    if (const dictionary* regionsPtr = hetPtr->findDict("regions"))
    {
        dictionaryKeyCheck::checkNamedBlocks
        (
            *regionsPtr,
            keys({"baseline", "cellZone", "range"})
        );
    }

    if (const dictionary* axesPtr = hetPtr->findDict("gradientAxes"))
    {
        dictionaryKeyCheck::checkNamedBlocks
        (
            *axesPtr,
            keys({"field", "beta", "scalingMin", "scalingMax", "variables"})
        );
    }
}


// The ionic model's own block and the blocks under it
void checkIonicBlock(const dictionary& coeffs)
{
    checkIonicHeterogeneity(coeffs);
    checkSubDict(coeffs, "singleCellStimulus", keys(singleCellStimulusKeys));
}


void checkConductionNetworks(const dictionary& networks)
{
    dictionaryKeyCheck::checkNamedBlocks
    (
        networks,
        keys({"conductionSystemDomain", "purkinjeGraphModelCoeffs"})
    );

    const wordList graphKeys
    ({
        "conductionSystemSolver", "graphFile", "chi", "cm",
        "purkinjeConductivity", "purkinjeFibreRadius", "purkinjeCV",
        "vm1DRest", "rootStimulus", "verificationModel",
        "activationVmThreshold", "apdNominal", "minimumDI90",
        "escapeInterval", "referenceConductance", "useEdgeConductance",
        "stimulus"
    });

    for (const entry& e : networks)
    {
        if (!e.isDict())
        {
            continue;
        }

        const dictionary* graphPtr =
            e.dict().findDict("purkinjeGraphModelCoeffs");

        if (!graphPtr)
        {
            continue;
        }

        dictionaryKeyCheck::checkKeys
        (
            *graphPtr,
            join({&graphKeys, &ionicKeys, &odeKeys})
        );

        checkIonicHeterogeneity(*graphPtr);

        // The network writes its outputs directly, without ionic{} or
        // activeTension{} groups
        checkSubDict(*graphPtr, "outputVariables", keys({"export", "debug"}));
        checkSubDict
        (
            *graphPtr,
            "rootStimulus",
            keys({"node", "startTime", "startTimeList", "duration", "intensity"})
        );
    }
}


void checkECGDomains(const dictionary& ecgDomains)
{
    dictionaryKeyCheck::checkNamedBlocks
    (
        ecgDomains,
        keys
        ({
            "ecgSolver", "electrodePositions", "sampling",
            "sigmaExtracellular", "torsoSurface", "verificationModel",
            "ecgVerificationModel", "personalizedTemplates", "coupling"
        }),
        wordRes(),
        wordList({"electrodePositions"})
    );

    for (const entry& e : ecgDomains)
    {
        if (!e.isDict() || e.keyword() == "electrodePositions")
        {
            continue;
        }

        checkSubDict(e.dict(), "sampling", keys({"start", "end", "deltaT"}));

        const dictionary* templatesPtr =
            e.dict().findDict("personalizedTemplates");

        if (!templatesPtr)
        {
            continue;
        }

        dictionaryKeyCheck::checkKeys
        (
            *templatesPtr,
            keys({"dt", "duration", "nBeats", "ionicModelConfig"})
        );

        const wordList stimulusKey({"singleCellStimulus"});

        if
        (
            const dictionary* configPtr =
                templatesPtr->findDict("ionicModelConfig")
        )
        {
            dictionaryKeyCheck::checkKeys
            (
                *configPtr,
                join({&ionicKeys, &odeKeys, &stimulusKey})
            );
            checkIonicBlock(*configPtr);
        }
    }
}


void checkTissueCoeffs(const dictionary& coeffs)
{
    checkSubDict
    (
        coeffs,
        "externalStimulus",
        keys
        ({
            "stimulusLocationMin", "stimulusLocationMax",
            "stimulusStartTime", "stimulusDuration", "stimulusIntensity",
            "stimulusLocationMinList", "stimulusLocationMaxList",
            "stimulusStartTimeList", "stimulusDurationList",
            "stimulusIntensityList"
        })
    );

    checkSubDict
    (
        coeffs,
        "bathPotentialDomain",
        keys
        ({
            "bathCellZones", "bathConductivityField", "groundPatches",
            "surfaceCurrentPatches", "interfaceConductivityInterpolation",
            "phiERefPoint", "phiEReferenceValue"
        })
    );

    if (const dictionary* couplingsPtr = coeffs.findDict("domainCouplings"))
    {
        dictionaryKeyCheck::checkNamedBlocks
        (
            *couplingsPtr,
            keys
            ({
                "electroDomainCoupler", "conductionNetworkDomain",
                "couplingMode", "pvjRadius", "pvjKernel",
                "pvjCouplingScheme", "rPvj", "debugCoupling",
                "verificationModel"
            })
        );
    }

    if
    (
        const dictionary* networksPtr =
            coeffs.findDict("conductionNetworkDomains")
    )
    {
        checkConductionNetworks(*networksPtr);
    }

    if (const dictionary* ecgPtr = coeffs.findDict("ecgDomains"))
    {
        checkECGDomains(*ecgPtr);
    }
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::electroPropertiesKeys::check
(
    const dictionary& electroProperties,
    const word& solverType
)
{
    // Coeffs blocks of solvers not selected may stay in the file
    wordRes topKeys(keys({"FoamFile", "myocardiumSolver"}));
    topKeys.append(wordRe(".*Coeffs", wordRe::REGEX));
    dictionaryKeyCheck::checkKeys(electroProperties, topKeys);

    const dictionary& coeffs = electroProperties.subDict(solverType + "Coeffs");

    if (solverType == "singleCellSolver")
    {
        dictionaryKeyCheck::checkKeys
        (
            coeffs,
            join({&singleCellKeys, &ionicKeys, &odeKeys})
        );
    }
    else
    {
        dictionaryKeyCheck::checkKeys
        (
            coeffs,
            join({&tissueKeys, &ionicKeys, &odeKeys})
        );
        checkTissueCoeffs(coeffs);
    }

    checkIonicBlock(coeffs);
    checkOutputVariables(coeffs);
}


// ************************************************************************* //
