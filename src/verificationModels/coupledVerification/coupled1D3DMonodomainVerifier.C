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
#include "mathematicalConstants.H"
#include "boundBox.H"
#include "monodomainVerification/manufacturedShellMonodomainReference.H"
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


//- Network volume pi*rho^2*L/2 summed per node from the graph inputs, not
//  from the network under test, so a wrong volume there fails the MMS
scalarField graphNodeVolumes(const conductionSystemDomain& graphDomain)
{
    const scalar radius(graphDomain.purkinjeFibreRadius());
    const scalarField& edgeLengths = graphDomain.edgeLengths();
    const labelList& edgeStart = graphDomain.edgeStartNodes();
    const labelList& edgeEnd = graphDomain.edgeEndNodes();

    label nNodes = 0;
    forAll(edgeStart, edgeI)
    {
        nNodes = max(nNodes, max(edgeStart[edgeI], edgeEnd[edgeI]) + 1);
    }

    scalarField nodeVolume(nNodes, 0.0);
    forAll(edgeLengths, edgeI)
    {
        const scalar halfEdge =
            0.5*constant::mathematical::pi*sqr(radius)*edgeLengths[edgeI];
        nodeVolume[edgeStart[edgeI]] += halfEdge;
        nodeVolume[edgeEnd[edgeI]] += halfEdge;
    }

    return nodeVolume;
}


} // End anonymous namespace


coupled1D3DMonodomainVerifier::coupled1D3DMonodomainVerifier
(
    const dictionary& dict
)
:
    couplingVerificationModel(dict),
    reference_(dict.lookupOrDefault<word>("reference", "discrete")),
    tissueField_(dict.lookupOrDefault<word>("tissueField", "box")),
    xOrigin_(dict.lookupOrDefault<scalar>("xOrigin", 0)),
    diagnosticKernels_
    (
        dict.lookupOrDefault<wordList>
        (
            "diagnosticKernels",
            wordList({"uniform", "linear", "gaussian", "cosine"})
        )
    ),
    diagnosticsWritten_(false),
    exactPrimarySource_(),
    exactSecondaryAppliedCurrent_(),
    continuumReady_(false),
    R_(0),
    kernel_(),
    Wc_(),
    Fc_(),
    refCells_(),
    refOmega_(),
    refOmegaF_(),
    graphVmOld_(),
    sourceBeforeDeposit_(),
    stepFile_(),
    maxChargeImbalance_(0),
    maxCurrentError_()
{
    if (reference_ != "discrete" && reference_ != "continuum")
    {
        FatalErrorInFunction
            << "reference must be discrete or continuum; found "
            << reference_ << "." << exit(FatalError);
    }

    if (tissueField_ != "box" && tissueField_ != "shell")
    {
        FatalErrorInFunction
            << "tissueField must be box or shell; found "
            << tissueField_ << "." << exit(FatalError);
    }

    if (reference_ == "discrete" && (tissueField_ != "box" || xOrigin_ != 0))
    {
        FatalErrorInFunction
            << "reference discrete supports only tissueField box with "
            << "xOrigin 0." << exit(FatalError);
    }
}


scalar coupled1D3DMonodomainVerifier::tissueF(const point& x) const
{
    if (tissueField_ == "shell")
    {
        static const shellMonodomainMMS::Params p =
            shellMonodomainMMS::defaultParams();
        return shellMonodomainMMS::F(x, p);
    }

    const scalar pi = constant::mathematical::pi;
    return
        Foam::cos(pi*x.x())*Foam::cos(2.0*pi*x.y())*Foam::cos(3.0*pi*x.z());
}


scalar coupled1D3DMonodomainVerifier::graphF(const point& x) const
{
    return Foam::cos(constant::mathematical::pi*(x.x() - xOrigin_));
}


pvjContinuumReference::support
coupled1D3DMonodomainVerifier::supportGeometry
(
    const fvMesh& mesh,
    const point& t,
    scalar R
) const
{
    pvjContinuumReference::support s;
    s.centre = t;
    s.normal = vector::zero;
    s.wallRadius = 0;
    s.kind = pvjContinuumReference::supportKind::interior;

    // Walls closer than this to the terminal count as the terminal's wall
    const scalar tol = 1e-9*max(R, scalar(1));

    if (tissueField_ == "shell")
    {
        const shellMonodomainMMS::Params p =
            shellMonodomainMMS::defaultParams();
        const scalar rt = mag(t);

        if (mag(rt - p.r1) <= tol)
        {
            if (p.r1 + R > p.r2)
            {
                FatalErrorInFunction
                    << "Support of the terminal at " << t << " reaches the "
                    << "outer sphere: r1 + R = " << p.r1 + R << " > r2 = "
                    << p.r2 << "." << exit(FatalError);
            }
            s.kind = pvjContinuumReference::supportKind::cavitySphere;
            s.normal = t/rt;
            s.wallRadius = p.r1;
            return s;
        }

        if (rt - R < p.r1 || rt + R > p.r2)
        {
            FatalErrorInFunction
                << "Terminal at " << t << " is neither on the inner sphere "
                << "nor at least R = " << R << " from both spheres."
                << exit(FatalError);
        }
        return s;
    }

    // Box: the mesh bounding box
    const boundBox bb(mesh.points(), true);
    label nOn = 0;
    for (direction d = 0; d < vector::nComponents; ++d)
    {
        const scalar toMin = t[d] - bb.min()[d];
        const scalar toMax = bb.max()[d] - t[d];

        if (mag(toMin) <= tol || mag(toMax) <= tol)
        {
            ++nOn;
            s.kind = pvjContinuumReference::supportKind::plane;
            s.normal = vector::zero;
            s.normal[d] = mag(toMin) <= tol ? 1 : -1;

            if ((mag(toMin) <= tol ? toMax : toMin) < R)
            {
                FatalErrorInFunction
                    << "Support of the terminal at " << t
                    << " reaches the opposite wall." << exit(FatalError);
            }
        }
        else if (toMin < R || toMax < R)
        {
            FatalErrorInFunction
                << "Support of the terminal at " << t << " is cut by a "
                << "wall the terminal is not on; the continuum reference "
                << "does not cover that geometry." << exit(FatalError);
        }
    }

    if (nOn > 1)
    {
        FatalErrorInFunction
            << "Terminal at " << t << " is on an edge or corner of the box."
            << exit(FatalError);
    }

    return s;
}


void coupled1D3DMonodomainVerifier::initialiseContinuum
(
    const fvMesh& mesh,
    const pointField& terminalLocations,
    bool implicitCoupling
)
{
    if (continuumReady_)
    {
        return;
    }

    const dictionary& couplingDict = dict().parent();
    R_ = couplingDict.lookupOrDefault<scalar>("pvjRadius", 0.5e-3);
    kernel_ = couplingDict.lookupOrDefault<word>("pvjKernel", "uniform");

    const pvjContinuumReference ref(kernel_, R_, dict());
    const pvjContinuumReference::scalarFunction F =
        [this](const point& x) { return tissueF(x); };

    const label nT = terminalLocations.size();
    Wc_.setSize(nT);
    Fc_.setSize(nT);
    refCells_.setSize(nT);
    refOmega_.setSize(nT);
    refOmegaF_.setSize(nT);

    forAll(terminalLocations, i)
    {
        const pvjContinuumReference::support s =
            supportGeometry(mesh, terminalLocations[i], R_);

        Wc_[i] = ref.supportIntegral(s);
        Fc_[i] = ref.supportIntegral(s, F)/Wc_[i];

        // Cell integrals of w F only for the local (implicit) operator
        ref.cellIntegrals
        (
            mesh,
            terminalLocations[i],
            implicitCoupling ? F : pvjContinuumReference::scalarFunction(),
            refCells_[i],
            refOmega_[i],
            refOmegaF_[i]
        );
    }

    maxCurrentError_.setSize(nT, 0.0);
    continuumReady_ = true;

    writeMapDiagnostics(mesh, terminalLocations);

    if (Pstream::master())
    {
        mkDir(mesh.time().globalPath()/"verification");
        stepFile_.reset
        (
            new OFstream
            (
                mesh.time().globalPath()/"verification"
              / "coupled1D3DMonodomain_steps.csv"
            )
        );
        stepFile_() << "t,dt,graphJunctionCharge,tissueJunctionCharge,"
            << "relativeImbalance";
        forAll(terminalLocations, i)
        {
            stepFile_() << ",Ih_" << i << ",Ic_" << i;
        }
        stepFile_() << nl;
    }
}


void coupled1D3DMonodomainVerifier::writeMapDiagnostics
(
    const fvMesh& mesh,
    const pointField& terminalLocations
) const
{
    const scalarField& V = mesh.V();
    const vectorField& C = mesh.C().primitiveField();

    volScalarField Ffield
    (
        IOobject
        (
            "coupled1D3DMonodomainVerifier:F3",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    );
    forAll(C, cellI)
    {
        Ffield.primitiveFieldRef()[cellI] = tissueF(C[cellI]);
    }

    const pvjContinuumReference::scalarFunction F =
        [this](const point& x) { return tissueF(x); };

    autoPtr<OFstream> osPtr;
    if (Pstream::master())
    {
        mkDir(mesh.time().globalPath()/"verification");
        osPtr.reset
        (
            new OFstream
            (
                mesh.time().globalPath()/"verification"
              / "coupled1D3DMonodomain_map.csv"
            )
        );
        osPtr()
            << "kernel,terminal,supportCells,Wh,Wc,normalisationError,"
            << "Fh,Fc,gatherError,depositL1,depositL2,depositL2rel,"
            << "meshIntegratedNormalisationError" << nl;
        osPtr().precision(12);
    }

    scalarField dS(mesh.nCells(), 0.0);
    scalarField Sc(mesh.nCells(), 0.0);

    forAll(diagnosticKernels_, kI)
    {
        const word& kernel = diagnosticKernels_[kI];
        const pvjMapper mapper(mesh, terminalLocations, R_, kernel);
        const pvjContinuumReference ref(kernel, R_, dict());

        scalarField Fh;
        mapper.gatherVm3DPvjs(Ffield, Fh);

        forAll(terminalLocations, i)
        {
            const pvjContinuumReference::support s =
                supportGeometry(mesh, terminalLocations[i], R_);
            const scalar Wc = ref.supportIntegral(s);
            const scalar Fc = ref.supportIntegral(s, F)/Wc;

            labelList cells;
            scalarField omega, omegaF;
            if (kernel == kernel_)
            {
                cells = refCells_[i];
                omega = refOmega_[i];
                omegaF = refOmegaF_[i];
            }
            else
            {
                ref.cellIntegrals
                (
                    mesh,
                    terminalLocations[i],
                    pvjContinuumReference::scalarFunction(),
                    cells,
                    omega,
                    omegaF
                );
            }

            const scalar Wh = mapper.supportNormalisation()[i];
            const labelList& hCells = mapper.terminalCellSets()[i];
            const scalarList& hWeights = mapper.terminalCellWeights()[i];

            // Source per unit current: discrete w/Wh, continuum cell average
            forAll(hCells, cI)
            {
                dS[hCells[cI]] += hWeights[cI]/Wh;
            }
            scalar sumOmega = 0;
            forAll(cells, cI)
            {
                const scalar avg = omega[cI]/(Wc*V[cells[cI]]);
                dS[cells[cI]] -= avg;
                Sc[cells[cI]] += avg;
                sumOmega += omega[cI];
            }

            scalar L1 = 0, L2 = 0, L2c = 0;
            forAll(hCells, cI)
            {
                const label c = hCells[cI];
                L1 += mag(dS[c])*V[c];
                L2 += sqr(dS[c])*V[c];
                L2c += sqr(Sc[c])*V[c];
                dS[c] = 0;
                Sc[c] = 0;
            }
            forAll(cells, cI)
            {
                const label c = cells[cI];
                L1 += mag(dS[c])*V[c];
                L2 += sqr(dS[c])*V[c];
                L2c += sqr(Sc[c])*V[c];
                dS[c] = 0;
                Sc[c] = 0;
            }

            label nCells = hCells.size();
            reduce(nCells, sumOp<label>());
            reduce(L1, sumOp<scalar>());
            reduce(L2, sumOp<scalar>());
            reduce(L2c, sumOp<scalar>());
            reduce(sumOmega, sumOp<scalar>());

            if (osPtr.valid())
            {
                osPtr()
                    << kernel << ',' << i << ',' << nCells << ','
                    << Wh << ',' << Wc << ',' << Wh/Wc - 1 << ','
                    << Fh[i] << ',' << Fc << ',' << Fh[i] - Fc << ','
                    << L1 << ',' << Foam::sqrt(L2) << ','
                    << Foam::sqrt(L2/max(L2c, VSMALL)) << ','
                    << sumOmega/Wc - 1 << nl;
            }
        }
    }
}


scalarField coupled1D3DMonodomainVerifier::continuumCurrent
(
    const pointField& terminalLocations,
    const scalarField& R_pvj,
    scalar t
) const
{
    const scalar timeFactor = Foam::sqrt(1.0 + t);

    scalarField current(terminalLocations.size());
    forAll(current, i)
    {
        current[i] =
            timeFactor*(graphF(terminalLocations[i]) - Fc_[i])/R_pvj[i];
    }

    return current;
}


void coupled1D3DMonodomainVerifier::updateContinuumSource
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    scalar evaluationTime,
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

    initialiseContinuum(mesh, terminalLocations, implicitCoupling);

    const scalarField R_pvj =
        terminalResistances(networkDomain, dict().parent());
    const scalarField Ic =
        continuumCurrent(terminalLocations, R_pvj, evaluationTime);

    if (phaseName == "secondary")
    {
        conductionSystemDomain* graphDomain =
            conductionDomainPtr(secondaryDomain);
        if (!graphDomain || !graphDomain->ionicModelPtr())
        {
            return;
        }

        // Network state before its advance, for the charge balance
        graphVmOld_ = graphDomain->membranePotential();

        exactSecondaryAppliedCurrent_.setSize
        (
            graphDomain->membranePotential().size(),
            0.0
        );
        exactSecondaryAppliedCurrent_ = 0.0;

        if (bidirectionalCoupling)
        {
            const scalarField nodeVolume(graphNodeVolumes(*graphDomain));
            const labelList& terminalNodes = graphDomain->terminalNodes();
            forAll(terminalNodes, i)
            {
                exactSecondaryAppliedCurrent_[terminalNodes[i]] -=
                    Ic[i]/nodeVolume[terminalNodes[i]];
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

    // Cell average of the continuum junction source at evaluationTime
    const scalarField& V = mesh.V();
    const scalar timeFactor = Foam::sqrt(1.0 + evaluationTime);

    exactPrimarySource_.setSize(mesh.nCells());
    exactPrimarySource_ = 0.0;

    forAll(terminalLocations, i)
    {
        const labelList& cells = refCells_[i];
        const scalarField& omega = refOmega_[i];
        const scalarField& omegaF = refOmegaF_[i];

        if (implicitCoupling)
        {
            // Local continuum operator: w (V1 - V3(x))/(R W)
            const scalar V1 = timeFactor*graphF(terminalLocations[i]);
            forAll(cells, cI)
            {
                exactPrimarySource_[cells[cI]] +=
                    (V1*omega[cI] - timeFactor*omegaF[cI])
                  / (R_pvj[i]*Wc_[i]*V[cells[cI]]);
            }
        }
        else
        {
            // Gathered continuum operator: I w/W
            forAll(cells, cI)
            {
                exactPrimarySource_[cells[cI]] +=
                    Ic[i]*omega[cI]/(Wc_[i]*V[cells[cI]]);
            }
        }
    }

    model->setManufacturedSourceTerm
    (
        exactPrimarySource_,
        volumeDomain.chi().value(),
        volumeDomain.Cm().value()
    );

    // The coupler deposits next; recordPrimaryDeposit takes the difference
    sourceBeforeDeposit_ = primaryDomain.sourceField().primitiveField();
}


void coupled1D3DMonodomainVerifier::recordPrimaryDeposit
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain,
    scalar t0,
    scalar dt
)
{
    if (!continuum() || !continuumReady_)
    {
        return;
    }

    networkCouplingEndpoint& networkDomain = requireNetworkDomain
    (
        secondaryDomain,
        "coupled1D3DMonodomainVerifier::recordPrimaryDeposit"
    );

    const fvMesh& mesh = primaryDomain.mesh();
    const pointField& terminalLocations = networkDomain.terminalLocations();
    const scalarField& V = mesh.V();
    const vectorField& C = mesh.C().primitiveField();
    const scalarField& S = primaryDomain.sourceField().primitiveField();

    // Charge entering the tissue: the deposited source, per support
    scalarField tissueCurrent(terminalLocations.size(), 0.0);
    forAll(S, cellI)
    {
        const scalar q = (S[cellI] - sourceBeforeDeposit_[cellI])*V[cellI];
        if (q == 0)
        {
            continue;
        }
        forAll(terminalLocations, i)
        {
            if (mag(C[cellI] - terminalLocations[i]) <= R_)
            {
                tissueCurrent[i] += q;
                break;
            }
        }
    }
    reduce
    (
        tissueCurrent.data(),
        tissueCurrent.size(),
        sumOp<scalar>(),
        UPstream::msgType(),
        UPstream::worldComm
    );
    const scalar tissueCharge = dt*sum(tissueCurrent);

    // Charge leaving the graph, from its own budget over the step:
    // chi Cm A sum l (V' - V) = -dt chi Cm A sum l Iion + dt A sum l Iapp
    //                           - (junction charge)
    scalar graphCharge = 0;
    if (conductionSystemDomain* graphDomain = conductionDomainPtr(secondaryDomain))
    {
        const scalarField& Vnew = graphDomain->membranePotential();
        const scalarField& Iion = graphDomain->ionicCurrent();
        scalarField Iapp(Vnew.size(), 0.0);
        graphDomain->assembleAppliedCurrent(t0, Iapp);

        const scalarField nodeVolume(graphNodeVolumes(*graphDomain));
        const scalar chiCm = graphDomain->chi()*graphDomain->Cm();

        if (graphVmOld_.size() == Vnew.size())
        {
            forAll(Vnew, nodeI)
            {
                graphCharge -=
                    nodeVolume[nodeI]
                   *(
                        chiCm*(Vnew[nodeI] - graphVmOld_[nodeI])
                      + dt*chiCm*Iion[nodeI]
                      - dt*Iapp[nodeI]
                    );
            }
        }
    }

    const scalarField R_pvj =
        terminalResistances(networkDomain, dict().parent());
    const scalarField Ic =
        continuumCurrent(terminalLocations, R_pvj, t0 + dt);

    // Relative to the charge crossing all junctions, not their net: the
    // terminals of a symmetric graph carry opposite currents
    const scalar imbalance =
        (graphCharge - tissueCharge)
       /max(dt*sum(mag(tissueCurrent)), VSMALL);
    maxChargeImbalance_ = max(maxChargeImbalance_, mag(imbalance));
    forAll(Ic, i)
    {
        maxCurrentError_[i] =
            max(maxCurrentError_[i], mag(tissueCurrent[i] - Ic[i]));
    }

    if (stepFile_.valid())
    {
        stepFile_().precision(12);
        stepFile_()
            << t0 + dt << ',' << dt << ',' << graphCharge << ','
            << tissueCharge << ',' << imbalance;
        forAll(Ic, i)
        {
            stepFile_() << ',' << tissueCurrent[i] << ',' << Ic[i];
        }
        stepFile_() << nl;
    }
}


void coupled1D3DMonodomainVerifier::postProcessContinuum
(
    tissueCouplingEndpoint& primaryDomain,
    electroDomainInterface& secondaryDomain
)
{
    networkCouplingEndpoint& networkDomain = requireNetworkDomain
    (
        secondaryDomain,
        "coupled1D3DMonodomainVerifier::postProcess"
    );

    const fvMesh& mesh = primaryDomain.mesh();
    const pointField& terminalLocations = networkDomain.terminalLocations();
    initialiseContinuum(mesh, terminalLocations, false);

    const scalar t = mesh.time().value();
    const scalar timeFactor = Foam::sqrt(1.0 + t);
    const scalarField& V = mesh.V();
    const vectorField& C = mesh.C().primitiveField();
    const scalarField& Vm = primaryDomain.Vm().primitiveField();

    // Final-time tissue error at the cell centres
    scalarField err(Vm.size());
    forAll(err, cellI)
    {
        err[cellI] = Vm[cellI] - timeFactor*tissueF(C[cellI]);
    }

    scalar sumAll = 0, volAll = 0;
    forAll(err, cellI)
    {
        sumAll += sqr(err[cellI])*V[cellI];
        volAll += V[cellI];
    }
    reduce(sumAll, sumOp<scalar>());
    reduce(volAll, sumOp<scalar>());
    const scalar L2all = Foam::sqrt(sumAll/max(volAll, VSMALL));

    // Far field: cells at least farFieldDistance from every terminal
    const scalar farDistance =
        dict().lookupOrDefault<scalar>("farFieldDistance", 0.25);
    scalar sumFar = 0, volFar = 0, maxFar = 0;
    label nFar = 0;
    forAll(err, cellI)
    {
        bool far = true;
        forAll(terminalLocations, i)
        {
            far = far && mag(C[cellI] - terminalLocations[i]) >= farDistance;
        }
        if (far)
        {
            sumFar += sqr(err[cellI])*V[cellI];
            volFar += V[cellI];
            maxFar = max(maxFar, mag(err[cellI]));
            ++nFar;
        }
    }
    reduce(sumFar, sumOp<scalar>());
    reduce(volFar, sumOp<scalar>());
    reduce(maxFar, maxOp<scalar>());
    reduce(nFar, sumOp<label>());
    const scalar L2far = Foam::sqrt(sumFar/max(volFar, VSMALL));

    // The junction's own view of both sides
    const pvjMapper mapper(mesh, terminalLocations, R_, kernel_);
    scalarField tissueGathered;
    mapper.gatherVm3DPvjs(primaryDomain.Vm(), tissueGathered);
    scalarField graphTerminal;
    networkDomain.terminalVm(graphTerminal);

    const scalarField R_pvj =
        terminalResistances(networkDomain, dict().parent());
    const scalarField Ic = continuumCurrent(terminalLocations, R_pvj, t);

    if (Pstream::master())
    {
        mkDir(mesh.time().globalPath()/"verification");
    }
    autoPtr<OFstream> osPtr;
    if (Pstream::master())
    {
        osPtr.reset
        (
            new OFstream
            (
                mesh.time().globalPath()/"verification"
              / "coupled1D3DMonodomain_local.csv"
            )
        );
        osPtr().precision(12);
        osPtr()
            << "t,terminal,kernel,R,rPvj,tissueL2,"
            << "supportCells,supportL2,supportLinf,"
            << "shellCells,shellL2,shellLinf,"
            << "graphTerminalError,tissueGatheredError,currentError,"
            << "Ic,maxStepCurrentError,maxChargeImbalance,"
            << "farFieldDistance,farCells,farL2,farLinf" << nl;
    }

    forAll(terminalLocations, i)
    {
        scalar s1 = 0, v1 = 0, m1 = 0;
        scalar s2 = 0, v2 = 0, m2 = 0;
        label n1 = 0, n2 = 0;

        forAll(err, cellI)
        {
            const scalar r = mag(C[cellI] - terminalLocations[i]);
            if (r <= R_)
            {
                s1 += sqr(err[cellI])*V[cellI];
                v1 += V[cellI];
                m1 = max(m1, mag(err[cellI]));
                ++n1;
            }
            if (r <= 2*R_)
            {
                s2 += sqr(err[cellI])*V[cellI];
                v2 += V[cellI];
                m2 = max(m2, mag(err[cellI]));
                ++n2;
            }
        }
        reduce(s1, sumOp<scalar>());
        reduce(v1, sumOp<scalar>());
        reduce(m1, maxOp<scalar>());
        reduce(n1, sumOp<label>());
        reduce(s2, sumOp<scalar>());
        reduce(v2, sumOp<scalar>());
        reduce(m2, maxOp<scalar>());
        reduce(n2, sumOp<label>());

        const scalar V1exact = timeFactor*graphF(terminalLocations[i]);
        const scalar V3gathered = timeFactor*Fc_[i];
        const scalar Ih = (graphTerminal[i] - tissueGathered[i])/R_pvj[i];

        if (osPtr.valid())
        {
            osPtr()
                << t << ',' << i << ',' << kernel_ << ',' << R_ << ','
                << R_pvj[i] << ',' << L2all << ','
                << n1 << ',' << Foam::sqrt(s1/max(v1, VSMALL)) << ',' << m1
                << ',' << n2 << ',' << Foam::sqrt(s2/max(v2, VSMALL)) << ','
                << m2 << ',' << graphTerminal[i] - V1exact << ','
                << tissueGathered[i] - V3gathered << ',' << Ih - Ic[i] << ','
                << Ic[i] << ',' << maxCurrentError_[i] << ','
                << maxChargeImbalance_ << ',' << farDistance << ','
                << nFar << ',' << L2far << ',' << maxFar << nl;
        }
    }
}


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
    if (continuum())
    {
        // Both sides at the new time level t0 + dt
        updateContinuumSource
        (
            primaryDomain,
            secondaryDomain,
            secondaryTime,
            implicitCoupling,
            bidirectionalCoupling,
            phaseName
        );
        return;
    }

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
        dict().parent().lookupOrDefault<word>("pvjKernel", "uniform")
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

            // The node volume pi*rho^2*L from the inputs, not from the
            // network under test, so a wrong volume there fails the MMS.
            const scalar radius(graphDomain->purkinjeFibreRadius());
            const scalarField& edgeLengths = graphDomain->edgeLengths();
            scalarField nodeVolume(exactSecondaryAppliedCurrent_.size(), 0.0);
            forAll(edgeLengths, edgeI)
            {
                const scalar halfEdge =
                    0.5*constant::mathematical::pi*sqr(radius)
                   *edgeLengths[edgeI];
                nodeVolume[graphDomain->edgeStartNodes()[edgeI]] += halfEdge;
                nodeVolume[graphDomain->edgeEndNodes()[edgeI]] += halfEdge;
            }

            const labelList& terminalNodes = graphDomain->terminalNodes();
            forAll(terminalNodes, i)
            {
                exactSecondaryAppliedCurrent_[terminalNodes[i]] -=
                    exactCurrent[i]/nodeVolume[terminalNodes[i]];
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

    // Support-local errors against the exact fields, and the map against
    // the continuum one, under either reference
    postProcessContinuum(primaryDomain, secondaryDomain);

    if (continuum())
    {
        return;
    }

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
        dict().parent().lookupOrDefault<word>("pvjKernel", "uniform")
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
