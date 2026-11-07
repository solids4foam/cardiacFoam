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

#include "manufacturedLeadFieldECGVerifier.H"

#include "addToRunTimeSelectionTable.H"
#include "verificationUtils.H"
#include "PstreamReduceOps.H"

namespace Foam
{

defineTypeNameAndDebug(manufacturedLeadFieldECGVerifier, 0);
addToRunTimeSelectionTable
(
    ecgVerificationModel,
    manufacturedLeadFieldECGVerifier,
    dictionary
);


const Enum<manufacturedLeadFieldECGVerifier::referenceType>
manufacturedLeadFieldECGVerifier::referenceTypeNames
({
    { referenceType::slab, "slab" },
    { referenceType::layered, "layered" }
});


namespace
{

// Layered: tolerance of the face normals and the cell extents (unit box)
const scalar layeredGeometryTol = 1e-9;

// Layered: relative spread of current times n_x over the surface-current
// faces and relative net current, both equal by construction up to round-off
const scalar layeredCurrentTol = 1e-12;

// Layered: relative anisotropy of a torso layer conductivity, isotropic by
// construction up to the round-off of the tensor read from the field
const scalar layeredIsotropyTol = 1e-10;

} // End anonymous namespace


manufacturedLeadFieldECGVerifier::manufacturedLeadFieldECGVerifier
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
    reference_(referenceType::slab),
    slabK_(0),
    layers_(),
    alpha_(0)
{
    read(dict);
}


bool manufacturedLeadFieldECGVerifier::read(const dictionary& dict)
{
    const dictionary& verifierDict = dict.subDict("verificationModel");
    const referenceType selected =
        referenceTypeNames.get("reference", verifierDict);

    reference_ = selected;

    if (reference_ == referenceType::layered)
    {
        const FixedList<scalar, 2> x
        (
            verifierDict.get<FixedList<scalar, 2>>("layerInterfaces")
        );
        if (!(x[0] < 0) || !(x[1] > 1))
        {
            FatalIOErrorInFunction(verifierDict)
                << "layerInterfaces " << x << " must be (xL xR) with xL < 0 "
                << "and xR > 1." << exit(FatalIOError);
        }
        layers_.xLeft = x[0];
        layers_.xRight = x[1];
    }
    else if (verifierDict.found("layerInterfaces"))
    {
        FatalIOErrorInFunction(verifierDict)
            << "layerInterfaces applies to reference layered only."
            << exit(FatalIOError);
    }

    return true;
}


word manufacturedLeadFieldECGVerifier::referenceName() const
{
    return referenceTypeNames[reference_];
}


void manufacturedLeadFieldECGVerifier::initialiseReference
(
    const leadFieldConductor& lfc
)
{
    if (reference_ == referenceType::slab)
    {
        const passiveConductor::torsoSettings& torso =
            lfc.conductor().torso();

        if
        (
            !torso.groundPatches.empty()
         || !torso.surfaceCurrentPatches.empty()
        )
        {
            FatalErrorInFunction
                << "The slab reference is insulated: it needs no ground "
                << "or surface-current patches, but has "
                << torso.groundPatches << " and "
                << torso.surfaceCurrentPatches << "."
                << exit(FatalError);
        }
    }

    slabK_ = manufacturedLeadFieldSlabK(lfc.Gi(), lfc.Ge());

    if (reference_ == referenceType::layered)
    {
        initialiseLayered(lfc);
    }
}


void manufacturedLeadFieldECGVerifier::initialiseLayered
(
    const leadFieldConductor& lfc
)
{
    const fvMesh& baseMesh = stateProvider_.baseMesh();

    const passiveConductor::torsoSettings& torso = lfc.conductor().torso();

    if (!torso.pinsReference())
    {
        FatalErrorInFunction
            << "The layered reference is pinned: it needs no ground patches, "
            << "but has " << torso.groundPatches << "." << exit(FatalError);
    }

    // alpha = I n_x on every surface-current face; no net current
    const wordList& patches = torso.surfaceCurrentPatches;
    const scalarList& currents = torso.surfaceCurrentValues;
    scalar alphaMin = GREAT;
    scalar alphaMax = -GREAT;
    scalar netCurrent = 0;
    scalar grossCurrent = 0;
    label nFaces = 0;
    forAll(patches, i)
    {
        const fvPatch& patch = baseMesh.boundary()[patches[i]];
        const vectorField nf(patch.nf());
        forAll(nf, faceI)
        {
            if (mag(mag(nf[faceI].x()) - 1) > layeredGeometryTol)
            {
                FatalErrorInFunction
                    << "The layered reference needs surface-current faces "
                    << "normal to x; patch '" << patches[i] << "' has the "
                    << "face normal " << nf[faceI] << "." << exit(FatalError);
            }
            const scalar a = currents[i]*nf[faceI].x();
            alphaMin = min(alphaMin, a);
            alphaMax = max(alphaMax, a);
            netCurrent += currents[i]*patch.magSf()[faceI];
            grossCurrent += mag(currents[i])*patch.magSf()[faceI];
            ++nFaces;
        }
    }
    reduce(alphaMin, minOp<scalar>());
    reduce(alphaMax, maxOp<scalar>());
    reduce(netCurrent, sumOp<scalar>());
    reduce(grossCurrent, sumOp<scalar>());
    reduce(nFaces, sumOp<label>());

    if (nFaces == 0)
    {
        alpha_ = 0;
    }
    else
    {
        if
        (
            alphaMax - alphaMin
          > layeredCurrentTol*max(mag(alphaMin), mag(alphaMax))
        )
        {
            FatalErrorInFunction
                << "The layered reference needs current times n_x equal on "
                << "every surface-current face; it ranges from " << alphaMin
                << " to " << alphaMax << " on " << patches << "."
                << exit(FatalError);
        }
        if (mag(netCurrent) > layeredCurrentTol*grossCurrent)
        {
            FatalErrorInFunction
                << "The layered reference is insulated apart from the "
                << "current along x: the net surface current " << netCurrent
                << " on " << patches << " must be zero." << exit(FatalError);
        }
        alpha_ = 0.5*(alphaMin + alphaMax);
    }

    // With a current the exact heart slope alpha/(Gi + Ge)_xx at the walls
    // carries the intracellular flux alpha Gi_xx/(Gi + Ge)_xx, which a
    // sealed Vm admits and an insulated bidomain intracellular potential
    // does not
    if (alpha_ != 0 && !lfc.conductor().totalWallConductivity())
    {
        FatalErrorInFunction
            << "The layered reference with a surface current (alpha "
            << alpha_ << ") is a monodomain solution; the bidomain conductor "
            << "insulates the intracellular potential at the heart wall."
            << exit(FatalError);
    }

    // (Gi + Ge)_xx; slabK_ has checked Gi and Ge
    layers_.sigmaHeart =
        verificationUtils::uniformTensor
        (
            lfc.Gi().primitiveField() + lfc.Ge().primitiveField(),
            "layered reference"
        ).xx();

    // Each torso cell lies in one layer; one isotropic value per layer
    const passiveConductor& conductor = lfc.conductor();
    const scalar xL = layers_.xLeft;
    const scalar xR = layers_.xRight;
    const scalar tol = layeredGeometryTol;
    const tensorField& sigma = conductor.sigmaTotal().primitiveField();
    const pointField& points = baseMesh.points();
    List<DynamicList<tensor>> layerSigma(4);

    forAll(sigma, cellI)
    {
        if (heartCells_[cellI])
        {
            continue;
        }

        scalar xMin = GREAT;
        scalar xMax = -GREAT;
        for (const label pointI : baseMesh.cellPoints()[cellI])
        {
            xMin = min(xMin, points[pointI].x());
            xMax = max(xMax, points[pointI].x());
        }

        label layer = -1;
        if (xMax <= xL + tol)
        {
            layer = 0;
        }
        else if (xMin >= xL - tol && xMax <= tol)
        {
            layer = 1;
        }
        else if (xMin >= 1 - tol && xMax <= xR + tol)
        {
            layer = 2;
        }
        else if (xMin >= xR - tol)
        {
            layer = 3;
        }
        else
        {
            FatalErrorInFunction
                << "Torso cell " << cellI << " spans " << xMin << " <= x <= "
                << xMax << ", across a layer interface of the layered "
                << "reference (" << xL << ", 0, 1, " << xR << ")."
                << exit(FatalError);
        }
        layerSigma[layer].append(sigma[cellI]);
    }

    const wordList layerNames
    ({
        "outer left", "inner left", "inner right", "outer right"
    });
    scalarList layerValue(4);
    forAll(layerSigma, layer)
    {
        const tensor s = verificationUtils::uniformTensor
        (
            tensorField(layerSigma[layer]),
            "layered reference, " + layerNames[layer] + " layer"
        );
        if (mag(s - s.xx()*tensor::I) > layeredIsotropyTol*mag(s))
        {
            FatalErrorInFunction
                << "The layered reference needs an isotropic torso "
                << "conductivity; the " << layerNames[layer] << " layer has "
                << s << "." << exit(FatalError);
        }
        layerValue[layer] = s.xx();
    }
    layers_.sigmaOuterLeft = layerValue[0];
    layers_.sigmaInnerLeft = layerValue[1];
    layers_.sigmaInnerRight = layerValue[2];
    layers_.sigmaOuterRight = layerValue[3];
}


scalar manufacturedLeadFieldECGVerifier::exactVm
(
    const point& p,
    const scalar t
) const
{
    return manufacturedLeadFieldSlabVm(p, t);
}


scalar manufacturedLeadFieldECGVerifier::referencePhi
(
    const point& p,
    const scalar t,
    const label
) const
{
    return
        reference_ == referenceType::slab
      ? manufacturedLeadFieldSlabPhi(p, t, slabK_)
      : manufacturedLeadFieldLayeredPhi(p, t, slabK_, alpha_, layers_);
}


void manufacturedLeadFieldECGVerifier::reportParameters(Ostream& os) const
{
    os  << ", k " << slabK_;
    if (reference_ == referenceType::layered)
    {
        os  << ", alpha " << alpha_ << ", interfaces " << layers_.xLeft
            << " " << layers_.xRight << ", sigma outer left "
            << layers_.sigmaOuterLeft << ", inner left "
            << layers_.sigmaInnerLeft << ", heart " << layers_.sigmaHeart
            << ", inner right " << layers_.sigmaInnerRight
            << ", outer right " << layers_.sigmaOuterRight;
    }
}


void manufacturedLeadFieldECGVerifier::writeParameters(Ostream& os) const
{
    os  << "# k " << slabK_ << "\n";
    if (reference_ == referenceType::layered)
    {
        os  << "# alpha " << alpha_ << "\n"
            << "# xLeft " << layers_.xLeft << "\n"
            << "# xRight " << layers_.xRight << "\n"
            << "# sigmaOuterLeft " << layers_.sigmaOuterLeft << "\n"
            << "# sigmaInnerLeft " << layers_.sigmaInnerLeft << "\n"
            << "# sigmaHeart " << layers_.sigmaHeart << "\n"
            << "# sigmaInnerRight " << layers_.sigmaInnerRight << "\n"
            << "# sigmaOuterRight " << layers_.sigmaOuterRight << "\n";
    }
}

} // End namespace Foam

// ************************************************************************* //
