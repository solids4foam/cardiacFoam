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
    obliqueSlabTransmission

Description
    Transmission manufactured solution on the oblique slab: a myocardium
    slab between the bath slabs bathLow and bathHigh, with planar interfaces
    s = s1 and s = s2, s = (cos alpha, sin alpha, 0).x, periodic along the
    interfaces. Fields from manufacturedTransmissionBathBidomainReference.H.

SourceFiles
    obliqueSlabTransmission.C

\*---------------------------------------------------------------------------*/

#ifndef obliqueSlabTransmission_H
#define obliqueSlabTransmission_H

#include "transmissionManufacturedSolution.H"
#include "manufacturedTransmissionBathBidomainReference.H"

namespace Foam
{

class obliqueSlabTransmission
:
    public transmissionManufacturedSolution
{
    transmissionMMS::Params p_;

    vector normal() const;

public:

    //- Exact-solution parameters from a verification dictionary.
    static transmissionMMS::Params readParams(const dictionary& dict);

    //- Check the reader against the reference defaults and samples
    static void selfTest();

    explicit obliqueSlabTransmission(const dictionary& dict);

    virtual word geometry() const override
    {
        return "obliqueSlab";
    }

    virtual word zoneName(const label zone) const override;

    virtual label zoneOf(const point& pt) const override
    {
        return transmissionMMS::zoneOf(pt, p_);
    }

    virtual label nSolutionD() const override
    {
        return 2;
    }

    virtual scalar Vm(const point& pt, const scalar t) const override
    {
        return transmissionMMS::Vm(pt, t, p_);
    }

    virtual scalar phiE(const point& pt, const scalar t, const label zone) const override
    {
        return transmissionMMS::phiE(pt, t, zone, p_);
    }

    virtual scalar phiI(const point& pt, const scalar t) const override
    {
        return transmissionMMS::phiI(pt, t, p_);
    }

    virtual vector gradPhiE(const point& pt, const scalar t, const label zone) const override
    {
        return transmissionMMS::gradPhiE(pt, t, zone, p_);
    }

    virtual vector gradVm(const point& pt, const scalar t) const override
    {
        return transmissionMMS::gradVm(pt, t, p_);
    }

    virtual scalar u1(const point& pt, const scalar t) const override
    {
        return transmissionMMS::u1(pt, t, p_);
    }

    virtual scalar u2(const point& pt, const scalar t) const override
    {
        return transmissionMMS::u2(pt, t, p_);
    }

    virtual symmTensor Gi(const point& pt) const override
    {
        return transmissionMMS::Gi(pt, p_);
    }

    virtual symmTensor Ge(const point& pt) const override
    {
        return transmissionMMS::Ge(pt, p_);
    }

    virtual scalar sigmaB(const point& pt, const label zone) const override
    {
        return transmissionMMS::sigmaB(pt, zone, p_);
    }

    virtual scalar divGiGradPhiI(const point& pt, const scalar t) const override
    {
        return transmissionMMS::divGiGradPhiI(pt, t, p_);
    }

    virtual scalar phiESource(const point& pt, const scalar t, const label zone) const override
    {
        return transmissionMMS::phiESource(pt, t, zone, p_);
    }

    virtual scalar interfaceFlux
    (
        const point& pt,
        const vector& n,
        const scalar t
    ) const override
    {
        return transmissionMMS::interfaceFlux(pt, n, t, p_);
    }

    virtual vector interfaceNormal(const point&) const override
    {
        return normal();
    }

    virtual label nMisplacedCells
    (
        const fvMesh& mesh,
        const labelUList& zone
    ) const override;

    virtual void validateInterface
    (
        const fvMesh& baseMesh,
        const labelUList& zone,
        const fvMesh& heartMesh
    ) const override;

    virtual labelList decomposition
    (
        const fvMesh& mesh,
        const labelUList& zone
    ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
