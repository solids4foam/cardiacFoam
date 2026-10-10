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

Class
    transmissionManufacturedSolution

Description
    Geometry of a bath-bidomain transmission manufactured solution: the
    exact fields, conductivities and sources of one generated reference
    header with its parameters bound, the zone of a point, the checks of a
    mesh against the analytic geometry, and the manual decomposition that
    setManufacturedFields writes.

    Zones are indexed 0 (bath on the first interface), 1 (myocardium) and 2
    (bath on the second interface). Selected by the verification dictionary
    entry geometry: obliqueSlab or sphericalShell.

SourceFiles
    transmissionManufacturedSolution.C

\*---------------------------------------------------------------------------*/

#ifndef transmissionManufacturedSolution_H
#define transmissionManufacturedSolution_H

#include "autoPtr.H"
#include "dictionary.H"
#include "error.H"
#include "fvMesh.H"
#include "labelList.H"
#include "symmTensor.H"

namespace Foam
{

class transmissionManufacturedSolution
{
public:

    static constexpr label heartZone = 1;

    virtual ~transmissionManufacturedSolution() = default;

    //- Select by the geometry entry, read the parameters and run the self-test
    static autoPtr<transmissionManufacturedSolution> New
    (
        const dictionary& verificationDict
    );

    virtual word geometry() const = 0;

    virtual word zoneName(const label zone) const = 0;

    virtual label zoneOf(const point& pt) const = 0;

    virtual label nSolutionD() const = 0;

    virtual scalar Vm(const point& pt, const scalar t) const = 0;
    virtual scalar phiE(const point& pt, const scalar t, const label zone) const = 0;
    virtual scalar phiI(const point& pt, const scalar t) const = 0;
    virtual vector gradPhiE(const point& pt, const scalar t, const label zone) const = 0;
    virtual vector gradVm(const point& pt, const scalar t) const = 0;
    virtual scalar u1(const point& pt, const scalar t) const = 0;
    virtual scalar u2(const point& pt, const scalar t) const = 0;
    virtual symmTensor Gi(const point& pt) const = 0;
    virtual symmTensor Ge(const point& pt) const = 0;
    virtual scalar sigmaB(const point& pt, const label zone) const = 0;
    virtual scalar divGiGradPhiI(const point& pt, const scalar t) const = 0;
    virtual scalar phiESource(const point& pt, const scalar t, const label zone) const = 0;

    //- n . Ge grad(phiE) of the heart formula
    virtual scalar interfaceFlux
    (
        const point& pt,
        const vector& n,
        const scalar t
    ) const = 0;

    //- Unit normal of the analytic interface nearest to pt, either sign
    virtual vector interfaceNormal(const point& pt) const = 0;

    //- Number of cells whose position contradicts their zone
    virtual label nMisplacedCells
    (
        const fvMesh& mesh,
        const labelUList& zone
    ) const = 0;

    //- Check the heart-bath faces and the heart extent against the analytic geometry
    virtual void validateInterface
    (
        const fvMesh& baseMesh,
        const labelUList& zone,
        const fvMesh& heartMesh
    ) const = 0;

    //- Four-rank manual decomposition with interfaces on internal and processor faces
    virtual labelList decomposition
    (
        const fvMesh& mesh,
        const labelUList& zone
    ) const = 0;

    //- Number of interface parts that carry a partial-current row, 0 for none
    virtual label nInterfaceParts() const
    {
        return 0;
    }

    //- Interface part of a heart-bath face centre
    virtual label interfacePart(const point& Cf) const
    {
        FatalErrorInFunction
            << geometry() << " has no interface parts." << exit(FatalError);
        return -1;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
