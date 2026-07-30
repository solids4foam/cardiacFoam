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

#include "coupled1D3DMonodomainVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "conductionSystemDomain.H"
#include "electroDomainCouplingEndpoints.H"
#include "electroDomainInterface.H"
#include "electroVolumeFieldDomain.H"
#include "ionicModel.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "pvjMapper.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(coupled1D3DMonodomainVerifier, 0);
addToRunTimeSelectionTable
(
    couplingVerificationModel,
    coupled1D3DMonodomainVerifier,
    dictionary
);

namespace
{

struct SourceStats
{
    scalar minValue = GREAT;
    scalar maxValue = -GREAT;
    scalar l1 = 0.0;
    scalar l2 = 0.0;
    scalar total = 0.0;
    scalar totalAbs = 0.0;
    label nonZero = 0;
};


SourceStats computeSourceStats
(
    const fvMesh& mesh,
    const scalarField& values
)
{
    const scalarField& volumes = mesh.V();

    SourceStats stats;
    scalar weightedAbs = 0.0;
    scalar weightedSq = 0.0;
    scalar totalVolume = 0.0;

    forAll(values, i)
    {
        const scalar source = values[i];
        stats.minValue = min(stats.minValue, source);
        stats.maxValue = max(stats.maxValue, source);
        weightedAbs += volumes[i]*mag(source);
        weightedSq += volumes[i]*source*source;
        totalVolume += volumes[i];
        stats.total += source*volumes[i];
        stats.totalAbs += mag(source)*volumes[i];

        if (mag(source) > SMALL)
        {
            ++stats.nonZero;
        }
    }

    reduce(stats.minValue, minOp<scalar>());
    reduce(stats.maxValue, maxOp<scalar>());
    reduce(weightedAbs, sumOp<scalar>());
    reduce(weightedSq, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());
    reduce(stats.total, sumOp<scalar>());
    reduce(stats.totalAbs, sumOp<scalar>());
    reduce(stats.nonZero, sumOp<label>());

    if (totalVolume > VSMALL)
    {
        stats.l1 = weightedAbs/totalVolume;
        stats.l2 = Foam::sqrt(weightedSq/totalVolume);
    }

    return stats;
}


networkCouplingEndpoint& requireNetworkDomain
(
    electroDomainInterface& secondaryDomain,
    const char* context
)
{
    auto* networkDomain =
        dynamic_cast<networkCouplingEndpoint*>(&secondaryDomain);

    if (!networkDomain)
    {
        FatalErrorInFunction
            << context << " requires a networkCouplingEndpoint "
            << "secondary domain."
            << exit(FatalError);
    }

    return *networkDomain;
}


electroVolumeFieldDomain& requireVolumeDomain
(
    tissueCouplingEndpoint& primaryDomain,
    const char* context
)
{
    auto* volumeDomain =
        dynamic_cast<electroVolumeFieldDomain*>(&primaryDomain);

    if (!volumeDomain)
    {
        FatalErrorInFunction
            << context << " requires an electroVolumeFieldDomain "
            << "primary domain."
            << exit(FatalError);
    }

    return *volumeDomain;
}


conductionSystemDomain* conductionDomainPtr
(
    electroDomainInterface& secondaryDomain
)
{
    return dynamic_cast<conductionSystemDomain*>(&secondaryDomain);
}


scalarField terminalX(const pointField& terminalLocations)
{
    scalarField x(terminalLocations.size());

    forAll(terminalLocations, i)
    {
        x[i] = terminalLocations[i].x();
    }

    return x;
}


scalarField terminalResistances
(
    const networkCouplingEndpoint& networkDomain,
    const dictionary& couplingDict
)
{
    if (const scalarField* pRes = networkDomain.terminalResistances())
    {
        return *pRes;
    }

    return scalarField
    (
        networkDomain.terminalLocations().size(),
        couplingDict.get<scalar>("rPvj")
    );
}


void computeExactVmField
(
    const fvMesh& mesh,
    const scalar timeValue,
    scalarField& VmExact
)
{
    const vectorField& centres = mesh.C().primitiveField();
    scalarField X(centres.component(vector::X));
    scalarField Y(centres.component(vector::Y));
    scalarField Z(centres.component(vector::Z));

    computeManufacturedV(VmExact, X, Y, Z, timeValue, 3);
}


void computeExactTerminalVm
(
    const pointField& terminalLocations,
    const scalar timeValue,
    scalarField& VmExact
)
{
    scalarField X(terminalX(terminalLocations));
    scalarField zeroY(X.size(), 0.0);
    scalarField zeroZ(X.size(), 0.0);

    computeManufacturedV(VmExact, X, zeroY, zeroZ, timeValue, 1);
}


void computeExactTerminalCurrent
(
    const fvMesh& mesh,
    const pvjMapper& mapper,
    const pointField& terminalLocations,
    const scalarField& R_pvj,
    scalar primaryTime,
    scalar secondaryTime,
    scalarField& terminalCurrent
)
{
    scalarField VmExact3D;
    computeExactVmField(mesh, primaryTime, VmExact3D);

    volScalarField VmExactField
    (
        IOobject
        (
            "coupled1D3DMonodomainVerifier:VmExactForCurrent",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("VmExactForCurrent", dimless, 0.0)
    );
    VmExactField.primitiveFieldRef() = VmExact3D;

    scalarField VmExact3DAtTerminals;
    mapper.gatherVm3DPvjs(VmExactField, VmExact3DAtTerminals);

    scalarField VmExact1D;
    computeExactTerminalVm(terminalLocations, secondaryTime, VmExact1D);

    terminalCurrent.setSize(VmExact1D.size(), 0.0);
    forAll(terminalCurrent, i)
    {
        terminalCurrent[i] =
            (VmExact1D[i] - VmExact3DAtTerminals[i])/R_pvj[i];
    }
}


} // End anonymous namespace


coupled1D3DMonodomainVerifier::coupled1D3DMonodomainVerifier
(
    const dictionary& dict
)
:
    couplingVerificationModel(dict),
    diagnosticsWritten_(false),
    exactPrimarySource_(),
    exactSecondaryAppliedCurrent_()
{}


void coupled1D3DMonodomainVerifier::preProcess
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain
)
{
    (void)primaryDomain;
    (void)secondaryDomain;
}


void coupled1D3DMonodomainVerifier::updateManufacturedSource
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    scalar primaryTime,
    scalar secondaryTime,
    bool implicitCoupling,
    bool bidirectionalCoupling,
    const word& phaseName
)
{
    networkCouplingEndpoint& networkDomain = requireNetworkDomain
    (
        secondaryDomain,
        "coupled1D3DMonodomainVerifier::updateManufacturedSource"
    );

    const fvMesh& mesh = primaryDomain.mesh();
    const pointField& terminalLocations = networkDomain.terminalLocations();
    pvjMapper exactMapper
    (
        mesh,
        terminalLocations,
        dict().parent().lookupOrDefault<scalar>("pvjRadius", 0.5e-3),
        dict().parent().lookupOrDefault<word>("pvjKernel", "uniform"),
        false
    );

    const scalarField R_pvj =
        terminalResistances(networkDomain, dict().parent());

    if (phaseName == "secondary")
    {
        conductionSystemDomain* graphDomain = conductionDomainPtr(secondaryDomain);
        if (!graphDomain || !graphDomain->ionicModelPtr())
        {
            return;
        }

        exactSecondaryAppliedCurrent_.setSize
        (
            graphDomain->membranePotential().size(),
            0.0
        );
        exactSecondaryAppliedCurrent_ = 0.0;

        if (bidirectionalCoupling)
        {
            scalarField exactCurrent;
            computeExactTerminalCurrent
            (
                mesh,
                exactMapper,
                terminalLocations,
                R_pvj,
                primaryTime,
                secondaryTime,
                exactCurrent
            );

            const labelList& terminalNodes = graphDomain->terminalNodes();
            forAll(terminalNodes, i)
            {
                exactSecondaryAppliedCurrent_[terminalNodes[i]] -= exactCurrent[i];
            }
        }

        graphDomain->ionicModelPtr()->setManufacturedSourceTerm
        (
            exactSecondaryAppliedCurrent_,
            graphDomain->chi(),
            graphDomain->Cm(),
            graphDomain->localStartNode()
        );
        return;
    }

    if (phaseName != "primary")
    {
        return;
    }

    electroVolumeFieldDomain& volumeDomain = requireVolumeDomain
    (
        primaryDomain,
        "coupled1D3DMonodomainVerifier::updateManufacturedSource"
    );

    ionicModel* model = primaryDomain.ionicModelPtr();
    if (!model)
    {
        return;
    }

    const volScalarField& sourceField = primaryDomain.sourceField();
    volScalarField exactSourceField
    (
        IOobject
        (
            "coupled1D3DMonodomainVerifier:manufacturedSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("manufacturedSource", sourceField.dimensions(), 0.0)
    );

    if (implicitCoupling)
    {
        const volScalarField* implicitSourceCoeff =
            primaryDomain.implicitSourceCoeffPtr();
        if (!implicitSourceCoeff)
        {
            FatalErrorInFunction
                << "Implicit manufactured PVJ source requires the primary "
                << "domain to expose an implicit source coefficient field."
                << exit(FatalError);
        }

        volScalarField exactImplicitCoeff
        (
            IOobject
            (
                "coupled1D3DMonodomainVerifier:manufacturedImplicitCoeff",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar
            (
                "manufacturedImplicitCoeff",
                implicitSourceCoeff->dimensions(),
                0.0
            )
        );

        scalarField VmExact1D;
        computeExactTerminalVm(terminalLocations, secondaryTime, VmExact1D);

        exactMapper.depositImplicitCoupling
        (
            VmExact1D,
            R_pvj,
            exactSourceField,
            exactImplicitCoeff
        );

        scalarField VmExact3D;
        computeExactVmField(mesh, secondaryTime, VmExact3D);

        scalarField& exactSource = exactSourceField.primitiveFieldRef();
        const scalarField& exactCoeff =
            exactImplicitCoeff.primitiveField();

        forAll(exactSource, cellI)
        {
            exactSource[cellI] -= exactCoeff[cellI]*VmExact3D[cellI];
        }
    }
    else
    {
        scalarField exactCurrent;
        computeExactTerminalCurrent
        (
            mesh,
            exactMapper,
            terminalLocations,
            R_pvj,
            primaryTime,
            secondaryTime,
            exactCurrent
        );

        exactMapper.depositCoupling(exactCurrent, exactSourceField);
    }

    exactPrimarySource_ = exactSourceField.primitiveField();
    model->setManufacturedSourceTerm
    (
        exactPrimarySource_,
        volumeDomain.chi().value(),
        volumeDomain.Cm().value()
    );
}


void coupled1D3DMonodomainVerifier::postProcess
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain
)
{
    if (diagnosticsWritten_)
    {
        return;
    }

    diagnosticsWritten_ = true;

    const fvMesh& mesh = primaryDomain.mesh();
    const volScalarField& sourceField = primaryDomain.sourceField();
    scalarField sourceValues(sourceField.primitiveField());
    bool implicitCoupling = false;

    if (const volScalarField* coeff = primaryDomain.implicitSourceCoeffPtr())
    {
        const scalarField& coeffValues = coeff->primitiveField();
        const scalarField& VmValues = primaryDomain.Vm().primitiveField();
        implicitCoupling = gMax(mag(coeffValues)) > SMALL;

        forAll(sourceValues, i)
        {
            sourceValues[i] -= coeffValues[i]*VmValues[i];
        }
    }

    SourceStats sourceStats = computeSourceStats(mesh, sourceValues);

    scalarField exactSourceValues(sourceValues.size(), 0.0);

    networkCouplingEndpoint& networkDomain = requireNetworkDomain
    (
        secondaryDomain,
        "coupled1D3DMonodomainVerifier::postProcess"
    );

    const pointField& terminalLocations = networkDomain.terminalLocations();
    pvjMapper exactMapper
    (
        mesh,
        terminalLocations,
        dict().parent().lookupOrDefault<scalar>("pvjRadius", 0.5e-3),
        dict().parent().lookupOrDefault<word>("pvjKernel", "uniform"),
        false
    );

    scalarField VmExact3D;
    computeExactVmField(mesh, mesh.time().value(), VmExact3D);

    volScalarField VmExactField
    (
        IOobject
        (
            "coupled1D3DMonodomainVerifier:VmExact",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("VmExactField", dimless, 0.0)
    );
    VmExactField.primitiveFieldRef() = VmExact3D;

    scalarField VmExact1D;
    computeExactTerminalVm(terminalLocations, mesh.time().value(), VmExact1D);

    const scalarField R_pvj = terminalResistances(networkDomain, dict().parent());

    volScalarField exactSourceField
    (
        IOobject
        (
            "coupled1D3DMonodomainVerifier:exactSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("exactSource", sourceField.dimensions(), 0.0)
    );

    if (implicitCoupling)
    {
        volScalarField exactImplicitCoeff
        (
            IOobject
            (
                "coupled1D3DMonodomainVerifier:exactImplicitCoeff",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar
            (
                "exactImplicitCoeff",
                primaryDomain.implicitSourceCoeffPtr()->dimensions(),
                0.0
            )
        );

        exactMapper.depositImplicitCoupling
        (
            VmExact1D,
            R_pvj,
            exactSourceField,
            exactImplicitCoeff
        );

        const scalarField& exactCoeff =
            exactImplicitCoeff.primitiveField();
        scalarField& exactSource = exactSourceField.primitiveFieldRef();

        forAll(exactSource, i)
        {
            exactSource[i] -= exactCoeff[i]*VmExact3D[i];
        }
    }
    else
    {
        scalarField VmExact3DAtTerminals;
        exactMapper.gatherVm3DPvjs(VmExactField, VmExact3DAtTerminals);

        scalarField exactCurrent(VmExact1D.size(), 0.0);
        forAll(exactCurrent, i)
        {
            exactCurrent[i] =
                (VmExact1D[i] - VmExact3DAtTerminals[i])/R_pvj[i];
        }

        exactMapper.depositCoupling(exactCurrent, exactSourceField);
    }

    exactSourceValues = exactSourceField.primitiveField();
    scalarField sourceError(sourceValues);
    sourceError -= exactSourceValues;

    const SourceStats exactStats = computeSourceStats(mesh, exactSourceValues);
    const SourceStats errorStats = computeSourceStats(mesh, sourceError);
    const label totalCells = globalManufacturedCellCount(mesh);

    if (Pstream::master())
    {
        mkDir(mesh.time().globalPath()/"verification");
        autoPtr<OFstream> osPtr
        (
            new OFstream
            (
                mesh.time().globalPath()/"verification"
              / "coupled1D3DMonodomain_diagnostics.csv"
            )
        );

        if (osPtr.valid())
        {
            *osPtr
                << "t,sourceMin,sourceMax,sourceL1,sourceL2,"
                << "totalCurrent,totalAbsCurrent,nonZeroCells,totalCells,"
                << "exactSourceL1,exactSourceL2,sourceErrorL1,"
                << "sourceErrorL2,totalExactAbsCurrent,totalAbsSourceError"
                << nl;

            *osPtr
                << mesh.time().value() << ','
                << sourceStats.minValue << ','
                << sourceStats.maxValue << ','
                << sourceStats.l1 << ','
                << sourceStats.l2 << ','
                << sourceStats.total << ','
                << sourceStats.totalAbs << ','
                << sourceStats.nonZero << ','
                << totalCells << ','
                << exactStats.l1 << ','
                << exactStats.l2 << ','
                << errorStats.l1 << ','
                << errorStats.l2 << ','
                << exactStats.totalAbs << ','
                << errorStats.totalAbs << nl;
        }
    }
}


} // End namespace Foam

// ************************************************************************* //
