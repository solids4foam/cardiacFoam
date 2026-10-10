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
    sphericalShellTransmission

Description
    Transmission manufactured solution on the spherical shell: a myocardium
    shell r1 < |x| < r2 around the bath ball bloodPool, inside the outer bath
    zone bath, with fibres rotating with radius about a fixed frame. Fields
    from manufacturedShellTransmissionBathBidomainReference.H. Mesh vertices
    must lie in the closed range of |x| of their cell's zone, so interface
    vertices lie on the spheres.

SourceFiles
    sphericalShellTransmission.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalShellTransmission_H
#define sphericalShellTransmission_H

#include "transmissionManufacturedSolution.H"
#include "manufacturedShellTransmissionBathBidomainReference.H"

namespace Foam
{

class sphericalShellTransmission
:
    public transmissionManufacturedSolution
{
    shellTransmissionMMS::Params p_;

    //- Absolute tolerance on |x|^2 for vertices on a sphere
    scalar rhoTolerance() const;

public:

    //- Exact-solution parameters from a verification dictionary.
    static shellTransmissionMMS::Params readParams(const dictionary& dict);

    //- Check the reader against the reference defaults and samples
    static void selfTest();

    explicit sphericalShellTransmission(const dictionary& dict);

    virtual word geometry() const override
    {
        return "sphericalShell";
    }

    virtual word zoneName(const label zone) const override;

    virtual label zoneOf(const point& pt) const override
    {
        return shellTransmissionMMS::zoneOf(pt, p_);
    }

    virtual bool acceptsSolutionD(const label nD) const override
    {
        return nD == 3;
    }

    virtual scalar Vm(const point& pt, const scalar t) const override
    {
        return shellTransmissionMMS::Vm(pt, t, p_);
    }

    virtual scalar phiE(const point& pt, const scalar t, const label zone) const override
    {
        return shellTransmissionMMS::phiE(pt, t, zone, p_);
    }

    virtual scalar phiI(const point& pt, const scalar t) const override
    {
        return shellTransmissionMMS::phiI(pt, t, p_);
    }

    virtual vector gradPhiE(const point& pt, const scalar t, const label zone) const override
    {
        return shellTransmissionMMS::gradPhiE(pt, t, zone, p_);
    }

    virtual vector gradVm(const point& pt, const scalar t) const override
    {
        return shellTransmissionMMS::gradVm(pt, t, p_);
    }

    virtual scalar u1(const point& pt, const scalar t) const override
    {
        return shellTransmissionMMS::u1(pt, t, p_);
    }

    virtual scalar u2(const point& pt, const scalar t) const override
    {
        return shellTransmissionMMS::u2(pt, t, p_);
    }

    virtual symmTensor Gi(const point& pt) const override
    {
        return shellTransmissionMMS::Gi(pt, p_);
    }

    virtual symmTensor Ge(const point& pt) const override
    {
        return shellTransmissionMMS::Ge(pt, p_);
    }

    virtual scalar sigmaB(const point& pt, const label zone) const override
    {
        return shellTransmissionMMS::sigmaB(pt, zone, p_);
    }

    virtual scalar divGiGradPhiI(const point& pt, const scalar t) const override
    {
        return shellTransmissionMMS::divGiGradPhiI(pt, t, p_);
    }

    virtual scalar phiESource(const point& pt, const scalar t, const label zone) const override
    {
        return shellTransmissionMMS::phiESource(pt, t, zone, p_);
    }

    virtual scalar interfaceFlux
    (
        const point& pt,
        const vector& n,
        const scalar t
    ) const override
    {
        return shellTransmissionMMS::interfaceFlux(pt, n, t, p_);
    }

    virtual vector interfaceNormal(const point& pt) const override;

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

    //- The eight parts of each sphere cut by three planes parallel to the coordinate planes
    virtual label nInterfaceParts() const override
    {
        return 8;
    }

    //- The planes are offset from the coordinate planes so that no face centre of a symmetric mesh lies on one
    virtual label interfacePart(const point& Cf) const override
    {
        return
            label(Cf.x() > 0.0123)
          + 2*label(Cf.y() > -0.0217)
          + 4*label(Cf.z() > 0.0311);
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
