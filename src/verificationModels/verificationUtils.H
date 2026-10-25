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

Description
    Shared utility helpers used by all electroVerificationModel concrete
    classes.  Placing these in a single header eliminates the copy-paste that
    previously existed across monodomainVerifier, bidomainVerifier, and
    bathBidomainVerifier.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef verificationUtils_H
#define verificationUtils_H

#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "fvMesh.H"
#include "fvcGrad.H"
#include "syncTools.H"
#include "volFields.H"
#include "wordList.H"

#include <cctype>

namespace Foam
{
namespace verificationUtils
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return the index of 'name' in 'names', or FatalError.
inline label requireFieldIndex
(
    const wordList& names,
    const word& name,
    const char* phase
)
{
    forAll(names, i)
    {
        if (names[i] == name)
        {
            return i;
        }
    }

    FatalErrorInFunction
        << "Required " << phase << " field '" << name
        << "' is missing from configured hook fields " << names
        << exit(FatalError);

    return -1;
}


//- Return true when the current time step is at (or past) the end time.
//  Uses a half-step tolerance to tolerate floating-point accumulation.
inline bool shouldReportManufacturedErrors(const volScalarField& Vm)
{
    const Time& time = Vm.mesh().time();
    return time.value() + 0.5*time.deltaTValue() >= time.endTime().value();
}


//- Return the global (across all MPI ranks) cell count of 'mesh'.
inline label globalManufacturedCellCount(const fvMesh& mesh)
{
    label totalCells = mesh.nCells();
    reduce(totalCells, sumOp<label>());
    return totalCells;
}


//- For a structured uniform mesh, return the number of cells per direction.
//  Falls back to totalCells for 1-D domains.
inline label structuredCellsPerDirection
(
    const label totalCells,
    const label dimension
)
{
    if (totalCells <= 0 || dimension <= 0)
    {
        return 0;
    }

    if (dimension == 1)
    {
        return totalCells;
    }

    return max
    (
        label(1),
        label(Foam::pow(scalar(totalCells), 1.0/scalar(dimension)) + 0.5)
    );
}


//- Grid spacing for a uniform [0,1]^d mesh with nPerDirection cells.
inline scalar structuredManufacturedDx(const label nPerDirection)
{
    if (nPerDirection <= 0)
    {
        return 0.0;
    }

    return 1.0/scalar(nPerDirection);
}


//- Compute parallel L1, L2, and Linf error norms between two scalar fields.
//
//  Returns Tuple2< Tuple2<L1, L2>, Linf >.
//
//  Guards:
//  - size mismatch → FatalError
//  - zero global cell count → all norms return 0
template<class FieldType1, class FieldType2>
inline Tuple2<Tuple2<scalar, scalar>, scalar> computeNorms
(
    const FieldType1& numeric,
    const FieldType2& exact
)
{
    if (numeric.size() != exact.size())
    {
        FatalErrorInFunction
            << "Cannot compute manufactured norms: field size mismatch ("
            << "numeric=" << numeric.size()
            << ", exact=" << exact.size() << ")."
            << exit(FatalError);
    }

    scalar sumAbs = 0.0;
    scalar sumSq  = 0.0;
    scalar maxAbs = 0.0;

    forAll(numeric, i)
    {
        const scalar diff = Foam::mag(numeric[i] - exact[i]);
        sumAbs += diff;
        sumSq  += diff*diff;
        maxAbs  = max(maxAbs, diff);
    }

    reduce(sumAbs, sumOp<scalar>());
    reduce(sumSq,  sumOp<scalar>());
    reduce(maxAbs, maxOp<scalar>());

    label n = numeric.size();
    reduce(n, sumOp<label>());

    const scalar denom = max(scalar(1), scalar(n));

    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>(sumAbs/denom, Foam::sqrt(sumSq/denom)),
        maxAbs
    );
}


//- Compute parallel, volume-weighted L1 and L2 error norms, together with
//  the (unweighted) Linf norm, for cell-centred finite-volume fields.
//
//  L1 = integral(|e| dV)/integral(dV)
//  L2 = sqrt(integral(e^2 dV)/integral(dV))
//
//  This normalization makes the result an RMS-type field error and leaves
//  the values unchanged on equal-volume meshes.  The mesh-free overload
//  above remains available for graph/node data without control volumes.
template<class FieldType1, class FieldType2>
inline Tuple2<Tuple2<scalar, scalar>, scalar> computeNorms
(
    const fvMesh& mesh,
    const FieldType1& numeric,
    const FieldType2& exact
)
{
    if
    (
        numeric.size() != exact.size()
     || numeric.size() != mesh.nCells()
    )
    {
        FatalErrorInFunction
            << "Cannot compute volume-weighted manufactured norms: size "
            << "mismatch (numeric=" << numeric.size()
            << ", exact=" << exact.size()
            << ", mesh cells=" << mesh.nCells() << ")."
            << exit(FatalError);
    }

    const scalarField& volumes = mesh.V();
    scalar weightedAbs = 0.0;
    scalar weightedSq = 0.0;
    scalar totalVolume = 0.0;
    scalar maxAbs = 0.0;

    forAll(numeric, i)
    {
        const scalar diff = Foam::mag(numeric[i] - exact[i]);
        weightedAbs += volumes[i]*diff;
        weightedSq += volumes[i]*diff*diff;
        totalVolume += volumes[i];
        maxAbs = max(maxAbs, diff);
    }

    reduce(weightedAbs, sumOp<scalar>());
    reduce(weightedSq, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());
    reduce(maxAbs, maxOp<scalar>());

    if (totalVolume <= VSMALL)
    {
        return Tuple2<Tuple2<scalar, scalar>, scalar>
        (
            Tuple2<scalar, scalar>(0.0, 0.0),
            maxAbs
        );
    }

    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>
        (
            weightedAbs/totalVolume,
            Foam::sqrt(weightedSq/totalVolume)
        ),
        maxAbs
    );
}


//- Split the volume-weighted L2 error into bulk and boundary-adjacent
//  contributions.
//
//  A cell is boundary-adjacent if it owns at least one face on a non-empty
//  boundary patch.  Both contributions are normalised by the *total* mesh
//  volume rather than by their own sub-volume, so they are parts of one norm
//  and satisfy
//
//      L2_total^2 = L2_bulk^2 + L2_boundary^2
//
//  The same split is reported for the field and gradient errors of the other
//  manufactured verifiers (bulkBoundaryRows), so solved-field and
//  reconstruction-only errors can be compared level by level on one ladder.
//
//  Returns Tuple2< Tuple2<L2_bulk, L2_boundary>, L2_total >.
template<class FieldType1, class FieldType2>
inline Tuple2<Tuple2<scalar, scalar>, scalar> computeBoundaryBulkNorms
(
    const fvMesh& mesh,
    const FieldType1& numeric,
    const FieldType2& exact
)
{
    if
    (
        numeric.size() != exact.size()
     || numeric.size() != mesh.nCells()
    )
    {
        FatalErrorInFunction
            << "Cannot compute bulk/boundary manufactured norms: size "
            << "mismatch (numeric=" << numeric.size()
            << ", exact=" << exact.size()
            << ", mesh cells=" << mesh.nCells() << ")."
            << exit(FatalError);
    }

    // Mark cells owning at least one physical boundary face. Empty and coupled
    // (processor, cyclic) patches are excluded, so the split does not depend
    // on the domain decomposition.
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    boolList isBoundaryCell(mesh.nCells(), false);

    forAll(bMesh, patchI)
    {
        const polyPatch& pp = bMesh[patchI];

        if (pp.empty() || pp.type() == "empty" || pp.coupled())
        {
            continue;
        }

        const labelUList& faceCells = pp.faceCells();
        forAll(faceCells, i)
        {
            isBoundaryCell[faceCells[i]] = true;
        }
    }

    const scalarField& volumes = mesh.V();

    scalar weightedSqBulk = 0.0;
    scalar weightedSqBound = 0.0;
    scalar totalVolume = 0.0;

    forAll(numeric, i)
    {
        const scalar diff = Foam::mag(numeric[i] - exact[i]);
        const scalar contrib = volumes[i]*diff*diff;

        if (isBoundaryCell[i])
        {
            weightedSqBound += contrib;
        }
        else
        {
            weightedSqBulk += contrib;
        }

        totalVolume += volumes[i];
    }

    reduce(weightedSqBulk, sumOp<scalar>());
    reduce(weightedSqBound, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());

    if (totalVolume <= VSMALL)
    {
        return Tuple2<Tuple2<scalar, scalar>, scalar>
        (
            Tuple2<scalar, scalar>(0.0, 0.0),
            0.0
        );
    }

    const scalar bulkL2 = Foam::sqrt(weightedSqBulk/totalVolume);
    const scalar boundL2 = Foam::sqrt(weightedSqBound/totalVolume);
    const scalar totalL2 =
        Foam::sqrt((weightedSqBulk + weightedSqBound)/totalVolume);

    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>(bulkL2, boundL2),
        totalL2
    );
}


//- Cells owning a face on a non-empty, non-coupled patch of the mesh.
//
//  Processor, cyclic and empty patches never mark a cell, so the set does not
//  depend on the decomposition. On a cellZone sub-mesh the exposed-faces patch
//  is a physical patch, so its cells are the wall cells.
inline boolList boundaryAdjacentCells(const fvMesh& mesh)
{
    boolList isBoundary(mesh.nCells(), false);

    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (pp.empty() || pp.type() == "empty" || pp.coupled())
        {
            continue;
        }

        for (const label cellI : pp.faceCells())
        {
            isBoundary[cellI] = true;
        }
    }

    return isBoundary;
}


//- Volume-weighted L1, L2 and Linf of the error over the masked cells, each
//  normalised by the total volume of the mesh.
//
//  The parts of a partition of the cells therefore add up: L1 = sum of the part
//  L1 values and L2^2 = sum of the part L2^2 values. Linf is the maximum over
//  the masked cells. Returns Tuple2< Tuple2<L1, L2>, Linf >.
template<class FieldType1, class FieldType2>
inline Tuple2<Tuple2<scalar, scalar>, scalar> computeMaskedNorms
(
    const fvMesh& mesh,
    const FieldType1& numeric,
    const FieldType2& exact,
    const boolList& mask
)
{
    if
    (
        numeric.size() != exact.size()
     || numeric.size() != mesh.nCells()
     || mask.size() != mesh.nCells()
    )
    {
        FatalErrorInFunction
            << "Cannot compute masked manufactured norms: size mismatch "
            << "(numeric=" << numeric.size() << ", exact=" << exact.size()
            << ", mask=" << mask.size() << ", mesh cells="
            << mesh.nCells() << ")."
            << exit(FatalError);
    }

    const scalarField& volumes = mesh.V();
    scalar weightedAbs = 0.0;
    scalar weightedSq = 0.0;
    scalar totalVolume = 0.0;
    scalar maxAbs = 0.0;

    forAll(numeric, i)
    {
        totalVolume += volumes[i];

        if (mask[i])
        {
            const scalar diff = Foam::mag(numeric[i] - exact[i]);
            weightedAbs += volumes[i]*diff;
            weightedSq += volumes[i]*diff*diff;
            maxAbs = max(maxAbs, diff);
        }
    }

    reduce(weightedAbs, sumOp<scalar>());
    reduce(weightedSq, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());
    reduce(maxAbs, maxOp<scalar>());

    if (totalVolume <= VSMALL)
    {
        return Tuple2<Tuple2<scalar, scalar>, scalar>
        (
            Tuple2<scalar, scalar>(0.0, 0.0),
            maxAbs
        );
    }

    return Tuple2<Tuple2<scalar, scalar>, scalar>
    (
        Tuple2<scalar, scalar>
        (
            weightedAbs/totalVolume,
            Foam::sqrt(weightedSq/totalVolume)
        ),
        maxAbs
    );
}


//- Cell gradient of an exact field, through the scheme named
//  grad(<fieldName>): exact values in the cells and on the faces of every
//  non-coupled patch, coupled patches evaluated from the neighbouring cells.
//
//  The field is unregistered and has calculated patches, so the result is the
//  reconstruction operator alone, with no boundary condition of the solved field
//  and no solve. exactAt(point) returns the exact value at a cell or face centre.
template<class ExactFn>
inline tmp<volVectorField> exactFieldGradient
(
    const fvMesh& mesh,
    const word& fieldName,
    const ExactFn& exactAt
)
{
    volScalarField q
    (
        IOobject
        (
            "exactForGradient(" + fieldName + ")",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        "calculated"
    );

    const vectorField& C = mesh.C().primitiveField();
    scalarField& qI = q.primitiveFieldRef();
    forAll(C, cellI)
    {
        qI[cellI] = exactAt(C[cellI]);
    }

    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& fp = mesh.boundary()[patchI];
        if (fp.coupled() || fp.size() == 0)
        {
            continue;
        }

        const vectorField& Cf = fp.Cf();
        scalarField values(Cf.size());
        forAll(Cf, faceI)
        {
            values[faceI] = exactAt(Cf[faceI]);
        }
        q.boundaryFieldRef()[patchI] == values;
    }

    q.correctBoundaryConditions();

    return fvc::grad(q, "grad(" + fieldName + ")");
}


//- Exact gradient by a fourth-order central difference of a closed-form
//  value function, for references that provide values only.
//
//  With step 1e-4 the truncation error is below 1e-12 for the smooth
//  manufactured fields; valueAt must use one smooth formula around the point.
template<class ValueFn>
inline vector centralDifferenceGradient(const ValueFn& valueAt, const point& p)
{
    const scalar h = 1e-4;
    vector g(Zero);

    for (direction d = 0; d < vector::nComponents; ++d)
    {
        vector e(Zero);
        e[d] = h;
        g[d] =
        (
            valueAt(p - 2*e) - 8*valueAt(p - e)
          + 8*valueAt(p + e) - valueAt(p + 2*e)
        )/(12*h);
    }

    return g;
}


//- Write one "name L1 L2 Linf" row.
inline void writeNormRow
(
    Ostream& os,
    const word& name,
    const Tuple2<Tuple2<scalar, scalar>, scalar>& norms
)
{
    os  << name << ' ' << norms.first().first() << ' '
        << norms.first().second() << ' ' << norms.second() << nl;
}


//- Rows for one field, computed on every rank: <name>_bulk and <name>_boundary
//  from the field error; grad<Name>, grad<Name>_bulk and grad<Name>_boundary from
//  the solved-field gradient error; the same three gradRecon rows from the
//  exact-field gradient error. Every part is normalised by the total volume of
//  the mesh, so total^2 = bulk^2 + boundary^2 in L2. gradStem names the
//  gradient rows when they differ from the field rows (gauge-shifted fields).
inline List<Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>>
bulkBoundaryRows
(
    const fvMesh& mesh,
    const word& name,
    const scalarField& numeric,
    const scalarField& exact,
    const vectorField& gradNumeric,
    const vectorField& gradRecon,
    const vectorField& gradExact,
    const word& gradStem = word::null
)
{
    const boolList boundary = boundaryAdjacentCells(mesh);
    boolList bulk(boundary.size());
    forAll(bulk, i)
    {
        bulk[i] = !boundary[i];
    }

    word cap(gradStem.empty() ? name : gradStem);
    cap[0] = char(std::toupper(cap[0]));

    return List<Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>>
    ({
        {name + "_bulk", computeMaskedNorms(mesh, numeric, exact, bulk)},
        {name + "_boundary", computeMaskedNorms(mesh, numeric, exact, boundary)},
        {"grad" + cap, computeNorms(mesh, gradNumeric, gradExact)},
        {"grad" + cap + "_bulk", computeMaskedNorms(mesh, gradNumeric, gradExact, bulk)},
        {"grad" + cap + "_boundary", computeMaskedNorms(mesh, gradNumeric, gradExact, boundary)},
        {"gradRecon" + cap, computeNorms(mesh, gradRecon, gradExact)},
        {"gradRecon" + cap + "_bulk", computeMaskedNorms(mesh, gradRecon, gradExact, bulk)},
        {"gradRecon" + cap + "_boundary", computeMaskedNorms(mesh, gradRecon, gradExact, boundary)}
    });
}


//- Cells sharing a face, internal or across a processor boundary, with a cell
//  on the other side of the heart-bath interface.
inline boolList interfaceBandCells(const fvMesh& mesh, const boolList& isHeart)
{
    boolList band(mesh.nCells(), false);
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    forAll(nei, faceI)
    {
        if (isHeart[own[faceI]] != isHeart[nei[faceI]])
        {
            band[own[faceI]] = true;
            band[nei[faceI]] = true;
        }
    }

    boolList nbrIsHeart;
    syncTools::swapBoundaryCellList(mesh, isHeart, nbrIsHeart);
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        if (!pp.coupled())
        {
            continue;
        }

        const labelUList& faceCells = pp.faceCells();
        const label start = pp.start() - mesh.nInternalFaces();
        forAll(faceCells, i)
        {
            if (isHeart[faceCells[i]] != nbrIsHeart[start + i])
            {
                band[faceCells[i]] = true;
            }
        }
    }

    return band;
}


//- Rows for a potential on the heart-and-bath mesh, computed on every rank. The
//  cells are split into the interface band (rows <name>_interface), the outer
//  boundary (boundary-adjacent cells outside the band) and the bulk, each
//  normalised by the total volume.
//  The exact field has a normal-derivative jump across the interface, so the
//  reconstruction-only rows exclude the interface band.
inline List<Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>>
globalPotentialRows
(
    const fvMesh& mesh,
    const boolList& isHeart,
    const word& name,
    const scalarField& numeric,
    const scalarField& exact,
    const vectorField& gradNumeric,
    const vectorField& gradRecon,
    const vectorField& gradExact
)
{
    const boolList band = interfaceBandCells(mesh, isHeart);
    const boolList adjacent = boundaryAdjacentCells(mesh);
    boolList boundary(band.size());
    boolList bulk(band.size());
    forAll(band, i)
    {
        boundary[i] = adjacent[i] && !band[i];
        bulk[i] = !adjacent[i] && !band[i];
    }

    word cap(name);
    cap[0] = char(std::toupper(cap[0]));

    return List<Tuple2<word, Tuple2<Tuple2<scalar, scalar>, scalar>>>
    ({
        {name + "_bulk", computeMaskedNorms(mesh, numeric, exact, bulk)},
        {name + "_boundary", computeMaskedNorms(mesh, numeric, exact, boundary)},
        {name + "_interface", computeMaskedNorms(mesh, numeric, exact, band)},
        {"grad" + cap, computeNorms(mesh, gradNumeric, gradExact)},
        {"grad" + cap + "_bulk", computeMaskedNorms(mesh, gradNumeric, gradExact, bulk)},
        {"grad" + cap + "_boundary", computeMaskedNorms(mesh, gradNumeric, gradExact, boundary)},
        {"grad" + cap + "_interface", computeMaskedNorms(mesh, gradNumeric, gradExact, band)},
        {"gradRecon" + cap + "_bulk", computeMaskedNorms(mesh, gradRecon, gradExact, bulk)},
        {"gradRecon" + cap + "_boundary", computeMaskedNorms(mesh, gradRecon, gradExact, boundary)}
    });
}


//- The gradient the extracellular-potential operator consumes: the jump-aware
//  phiEInterfaceGradient when the bath domain publishes it, else fvc::grad.
inline tmp<volVectorField> consumedPhiEGradient(const volScalarField& phiE)
{
    const fvMesh& mesh = phiE.mesh();
    if (mesh.foundObject<volVectorField>("phiEInterfaceGradient"))
    {
        return tmp<volVectorField>::New
        (
            mesh.lookupObject<volVectorField>("phiEInterfaceGradient")
        );
    }
    return fvc::grad(phiE);
}


//- "tet" when every cell has four faces, "hex" when every cell has six,
//  "mixed" otherwise, over all ranks.
inline word meshCellKind(const fvMesh& mesh)
{
    bool allTet = true;
    bool allHex = true;
    for (const cell& c : mesh.cells())
    {
        allTet = allTet && c.size() == 4;
        allHex = allHex && c.size() == 6;
    }
    reduce(allTet, andOp<bool>());
    reduce(allHex, andOp<bool>());
    return allTet ? "tet" : (allHex ? "hex" : "mixed");
}


//- Human-readable dimension tag ("1D", "2D", "3D", or "unknown").
inline word dimensionName(const label dimension)
{
    return
        dimension == 1 ? "1D"
      : dimension == 2 ? "2D"
      : dimension == 3 ? "3D"
                       : "unknown";
}


inline scalar computeVolumeMeanError
(
    const fvMesh& mesh,
    const scalarField& numeric,
    const scalarField& exact
)
{
    const scalarField& volumes = mesh.V();

    scalar weightedError = 0.0;
    scalar totalVolume = 0.0;

    forAll(numeric, i)
    {
        weightedError += volumes[i]*(numeric[i] - exact[i]);
        totalVolume += volumes[i];
    }

    reduce(weightedError, sumOp<scalar>());
    reduce(totalVolume, sumOp<scalar>());

    if (totalVolume <= VSMALL)
    {
        return 0.0;
    }

    return weightedError/totalVolume;
}




//- Return the global mean of a spatially constant tensor field, or FatalError.
inline tensor uniformTensor(const tensorField& sigma, const word& owner)
{
    tensor sigmaSum(Zero);
    forAll(sigma, cellI)
    {
        sigmaSum += sigma[cellI];
    }
    reduce(sigmaSum, sumOp<tensor>());

    label nCells = sigma.size();
    reduce(nCells, sumOp<label>());

    if (nCells <= 0)
    {
        FatalErrorInFunction
            << owner << ": cannot validate conductivity on an empty mesh."
            << exit(FatalError);
    }

    const tensor meanSigma = sigmaSum/scalar(nCells);
    scalar maxDeviation = 0.0;
    forAll(sigma, cellI)
    {
        maxDeviation = max(maxDeviation, mag(sigma[cellI] - meanSigma));
    }
    reduce(maxDeviation, maxOp<scalar>());

    const scalar tolerance = 1e-10*max(scalar(1), mag(meanSigma));

    if (maxDeviation > tolerance)
    {
        FatalErrorInFunction
            << owner << " assumes a spatially constant conductivity tensor. "
            << "Maximum cellwise deviation from the global mean is "
            << maxDeviation << ", exceeding tolerance " << tolerance << "."
            << exit(FatalError);
    }

    return meanSigma;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace verificationUtils
} // End namespace Foam

#endif

// ************************************************************************* //
