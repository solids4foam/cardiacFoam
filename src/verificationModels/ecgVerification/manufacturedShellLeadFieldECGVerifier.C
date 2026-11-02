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

#include "manufacturedShellLeadFieldECGVerifier.H"

#include "addToRunTimeSelectionTable.H"
#include "emptyFvPatch.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(manufacturedShellLeadFieldECGVerifier, 0);
addToRunTimeSelectionTable
(
    ecgVerificationModel,
    manufacturedShellLeadFieldECGVerifier,
    dictionary
);


namespace
{

// sphericalShell: relative tolerance of the case conductivities against the
// solution's at the cell centres
const scalar shellConductivityTol = 1e-10;


// sphericalShell: the zone of the transmission solution for every cell,
// from the cellZones bloodPool, myocardium and bath of the mesh (fatal if a
// cell is in none). A cell centre may lie outside the radial range of its
// zone, where a polyhedral interface cuts the sphere, by at most the cell
// size cbrt(V) (the zone formulas extend smoothly across the spheres).
labelList shellZones(const fvMesh& mesh)
{
    using namespace shellTransmissionMMS;

    const Params p = defaultParams();
    const wordList names({"bloodPool", "myocardium", "bath"});
    const label zoneIds[3] = {bloodPool, myocardium, bath};

    labelList zone(mesh.nCells(), -1);
    forAll(names, i)
    {
        const label zoneI = mesh.cellZones().findZoneID(names[i]);
        if (zoneI < 0)
        {
            FatalErrorInFunction
                << "The sphericalShell solution needs the cellZone "
                << names[i] << "." << exit(FatalError);
        }
        for (const label cellI : mesh.cellZones()[zoneI])
        {
            zone[cellI] = zoneIds[i];
        }
    }

    const vectorField& C = mesh.C().primitiveField();
    const scalarField& V = mesh.V().field();
    label unassigned = 0;
    label misplaced = 0;
    forAll(zone, cellI)
    {
        if (zone[cellI] < 0)
        {
            ++unassigned;
            continue;
        }
        const label z = zone[cellI];
        const scalar r = mag(C[cellI]);
        const scalar rMin =
            z == bloodPool ? 0 : (z == myocardium ? p.r1 : p.r2);
        const scalar rMax =
            z == bloodPool ? p.r1 : (z == myocardium ? p.r2 : GREAT);
        const scalar size = Foam::cbrt(V[cellI]);
        if (r < rMin - size || r > rMax + size)
        {
            ++misplaced;
        }
    }
    reduce(unassigned, sumOp<label>());
    reduce(misplaced, sumOp<label>());
    if (unassigned > 0 || misplaced > 0)
    {
        FatalErrorInFunction
            << unassigned << " cells lie in no sphericalShell cellZone and "
            << misplaced << " cells lie farther than their own size outside "
            << "the radial range of their zone." << exit(FatalError);
    }

    return zone;
}

} // End anonymous namespace


manufacturedShellLeadFieldECGVerifier::manufacturedShellLeadFieldECGVerifier
(
    const electroStateProvider& stateProvider,
    const word& domainName,
    const dictionary& dict,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
:
    leadFieldECGVerifierBase
    (
        stateProvider, domainName, electrodeNames, electrodePositions
    ),
    shellParams_(shellTransmissionMMS::defaultParams()),
    leadFields_(),
    shellZones_()
{
    read(dict);
}


bool manufacturedShellLeadFieldECGVerifier::read(const dictionary& dict)
{
    const dictionary& verifierDict = dict.subDict("verificationModel");

    for (const entry& e : verifierDict)
    {
        if (e.keyword() != "type")
        {
            FatalIOErrorInFunction(verifierDict)
                << type() << " takes no key beyond type (the solution is the "
                << "spherical shell); found " << e.keyword() << "."
                << exit(FatalIOError);
        }
    }

    return true;
}


word manufacturedShellLeadFieldECGVerifier::referenceName() const
{
    return "sphericalShell";
}


void manufacturedShellLeadFieldECGVerifier::writeManufacturedFields
(
    const fvMesh& mesh
) const
{
    using namespace shellTransmissionMMS;

    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "The sphericalShell lead-field conductivities are written on "
            << "the undecomposed case; run in serial." << exit(FatalError);
    }

    const dimensionSet conductivityDim(-1, -3, 3, 0, 0, 2, 0);
    const vectorField& C = mesh.C().primitiveField();
    const unsigned int oldPrecision = IOstream::defaultPrecision(17);

    const auto io = [&mesh](const word& name)
    {
        return IOobject(name, "0", mesh, IOobject::NO_READ, IOobject::NO_WRITE);
    };

    const auto copyInternalToPatches = [&mesh](auto& field)
    {
        forAll(mesh.boundary(), patchI)
        {
            if (polyPatch::constraintType(mesh.boundaryMesh()[patchI].type()))
            {
                continue;
            }
            field.boundaryFieldRef()[patchI] ==
                field.boundaryField()[patchI].patchInternalField();
        }
    };

    for (label i = 0; i < 2; ++i)
    {
        volSymmTensorField G
        (
            io
            (
                i == 0
              ? "ConductivityIntracellular"
              : "ConductivityExtracellular"
            ),
            mesh,
            dimensionedSymmTensor(conductivityDim, Zero),
            "calculated"
        );
        forAll(C, cellI)
        {
            G[cellI] =
                i == 0
              ? Gi(C[cellI], shellParams_)
              : Ge(C[cellI], shellParams_);
        }
        copyInternalToPatches(G);
        G.write();
        Info<< "Wrote " << G.name() << " in 0." << endl;
    }

    volScalarField sigma
    (
        io("bodyAndOrgansConductivity"),
        mesh,
        dimensionedScalar(conductivityDim, Zero),
        "calculated"
    );
    const labelList zone(shellZones(mesh));
    forAll(C, cellI)
    {
        sigma[cellI] = sigmaB
        (
            C[cellI],
            zone[cellI] == myocardium ? label(bath) : zone[cellI],
            shellParams_
        );
    }
    copyInternalToPatches(sigma);
    sigma.write();
    Info<< "Wrote " << sigma.name() << " in 0." << endl;

    IOstream::defaultPrecision(oldPrecision);
}


void manufacturedShellLeadFieldECGVerifier::initialiseReference
(
    const leadFieldConductor& lfc
)
{
    using namespace shellTransmissionMMS;

    const passiveConductor& conductor = lfc.conductor();
    const fvMesh& baseMesh = stateProvider_.baseMesh();

    const passiveConductor::torsoSettings& torso = conductor.torso();

    if
    (
        !torso.groundPatches.empty()
     || !torso.surfaceCurrentPatches.empty()
    )
    {
        FatalErrorInFunction
            << "The sphericalShell reference is pinned and takes its outer "
            << "current from the solution: it needs no ground or "
            << "surface-current patches, but has " << torso.groundPatches
            << " and " << torso.surfaceCurrentPatches << "."
            << exit(FatalError);
    }
    if (conductor.totalWallConductivity())
    {
        FatalErrorInFunction
            << "The sphericalShell reference seals the intracellular "
            << "potential at the heart wall: it needs the bidomain conductor."
            << exit(FatalError);
    }

    // Zones from the cellZones against the conductor's heart cells, and the
    // torso conductivity against sigmaB
    shellZones_ = shellZones(baseMesh);
    const vectorField& C = baseMesh.C().primitiveField();
    const tensorField& sigma = conductor.sigmaTotal().primitiveField();
    label misplaced = 0;
    scalar worst = 0;
    forAll(C, cellI)
    {
        const label zone = shellZones_[cellI];
        if (heartCells_[cellI] != (zone == myocardium))
        {
            ++misplaced;
        }
        else if (!heartCells_[cellI])
        {
            const tensor exact = sigmaB(C[cellI], zone, shellParams_)*tensor::I;
            worst = max(worst, mag(sigma[cellI] - exact)/mag(exact));
        }
    }

    // Gi and Ge of the heart cells
    const vectorField& Ch = lfc.Gi().mesh().C().primitiveField();
    const tensorField& GiI = lfc.Gi().primitiveField();
    const tensorField& GeI = lfc.Ge().primitiveField();
    forAll(Ch, cellI)
    {
        const tensor gi = Gi(Ch[cellI], shellParams_) & tensor::I;
        const tensor ge = Ge(Ch[cellI], shellParams_) & tensor::I;
        worst = max(worst, mag(GiI[cellI] - gi)/mag(gi));
        worst = max(worst, mag(GeI[cellI] - ge)/mag(ge));
    }
    reduce(misplaced, sumOp<label>());
    reduce(worst, maxOp<scalar>());

    if (misplaced > 0)
    {
        FatalErrorInFunction
            << misplaced << " cells of the cellZone myocardium differ from "
            << "the conductor's heart cells." << exit(FatalError);
    }
    if (worst > shellConductivityTol)
    {
        FatalErrorInFunction
            << "The case conductivities differ from the sphericalShell "
            << "solution's at the cell centres by up to " << worst
            << " (relative); write them with setManufacturedFields."
            << exit(FatalError);
    }

    forAll(electrodePositions_, e)
    {
        if (zoneOf(electrodePositions_[e], shellParams_) == myocardium)
        {
            FatalErrorInFunction
                << "Electrode " << electrodeNames_[e] << " at "
                << electrodePositions_[e] << " lies in the myocardium of the "
                << "sphericalShell solution." << exit(FatalError);
        }
    }

    leadFields_.setSize(electrodeNames_.size());
    forAll(electrodeNames_, e)
    {
        leadFields_.set
        (
            e,
            new volScalarField
            (
                IOobject
                (
                    "leadField_" + electrodeNames_[e],
                    baseMesh.time().constant(),
                    baseMesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                baseMesh
            )
        );
    }
}


tmp<volScalarField> manufacturedShellLeadFieldECGVerifier::shellSource
(
    const scalar t,
    const dimensionSet& dims
) const
{
    using namespace shellTransmissionMMS;

    const fvMesh& baseMesh = stateProvider_.baseMesh();

    tmp<volScalarField> tS
    (
        new volScalarField
        (
            IOobject
            (
                "shellPhiESource",
                baseMesh.time().timeName(),
                baseMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            baseMesh,
            dimensionedScalar(dims, 0),
            "zeroGradient"
        )
    );
    scalarField& S = tS.ref().primitiveFieldRef();

    // div(K grad phi) = b + S in every zone
    const vectorField& C = baseMesh.C().primitiveField();
    forAll(C, cellI)
    {
        S[cellI] = phiESource(C[cellI], t, shellZones_[cellI], shellParams_);
    }

    // The exact outward current g on the insulated outer boundary moves into
    // the owner cell's balance: S -= g |Sf|/V
    const scalarField& V = baseMesh.V().field();
    forAll(baseMesh.boundary(), patchI)
    {
        const fvPatch& patch = baseMesh.boundary()[patchI];
        if (patch.coupled() || isA<emptyFvPatch>(patch))
        {
            continue;
        }
        const vectorField nf(patch.nf());
        const vectorField& Cf = patch.Cf();
        const scalarField& magSf = patch.magSf();
        const labelUList& cells = patch.faceCells();
        forAll(cells, faceI)
        {
            const point& pf = Cf[faceI];
            const scalar g =
                sigmaB(pf, bath, shellParams_)
               *(nf[faceI] & gradPhiE(pf, t, bath, shellParams_));
            S[cells[faceI]] -= g*magSf[faceI]/V[cells[faceI]];
        }
    }

    // The pure-Neumann problem needs a source of zero integral
    S -= gSum(V*S)/gSum(V);

    return tS;
}


scalar manufacturedShellLeadFieldECGVerifier::exactVm
(
    const point& p,
    const scalar t
) const
{
    return shellTransmissionMMS::Vm(p, t, shellParams_);
}


scalar manufacturedShellLeadFieldECGVerifier::referencePhi
(
    const point& p,
    const scalar t,
    const label zone
) const
{
    using namespace shellTransmissionMMS;

    return phiE(p, t, zone < 0 ? zoneOf(p, shellParams_) : zone, shellParams_);
}


tmp<volScalarField> manufacturedShellLeadFieldECGVerifier::staticSource
(
    const volScalarField& VmExact,
    const scalar t,
    scalarField& leadFieldStatic
) const
{
    const passiveConductor& conductor = leadField().conductor().conductor();

    tmp<volScalarField> tS
    (
        shellSource(t, conductor.heartSource(VmExact)().dimensions())
    );
    const scalarField& S = tS().primitiveField();
    const scalarField& V = stateProvider_.baseMesh().V().field();
    forAll(leadFieldStatic, e)
    {
        leadFieldStatic[e] += gSum(leadFields_[e].primitiveField()*V*S);
    }

    return tS;
}

} // End namespace Foam

// ************************************************************************* //
