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
#include "volFields.H"
#include "wordList.H"

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




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace verificationUtils
} // End namespace Foam

#endif

// ************************************************************************* //
