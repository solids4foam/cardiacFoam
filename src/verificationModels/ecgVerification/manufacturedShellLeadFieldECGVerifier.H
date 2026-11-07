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
    manufacturedShellLeadFieldECGVerifier

Description
    Lead-field ECG verifier against the spherical-shell transmission
    solution of shell_mms.py (heart shell between the bath zones bloodPool
    and bath), with the channels of leadFieldECGVerifierBase. Its
    verificationModel takes no keys beyond type.

    The conductor must be the bidomain one, pinned, with no ground or
    surface-current patches, and its zones, Gi, Ge and torso conductivity
    must be those of the solution at the cell centres.
    writeManufacturedFields (setManufacturedFields) writes them:
    0/ConductivityIntracellular, 0/ConductivityExtracellular and
    0/bodyAndOrgansConductivity (the outer-bath formula in the myocardium,
    where the conductor does not read it), each cell with the formula of its
    cellZone. The phiE source of the solution and its exact outward current
    on the outer boundary enter both static channels as one extra source,
    without its volume mean: the direct solve adds it to b, the lead field
    adds sum z_e V s.

SourceFiles
    manufacturedShellLeadFieldECGVerifier.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedShellLeadFieldECGVerifier_H
#define manufacturedShellLeadFieldECGVerifier_H

#include "leadFieldECGVerifierBase.H"
#include "ecgVerification/manufacturedShellLeadFieldECGReference.H"
#include "PtrList.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedShellLeadFieldECGVerifier
:
    public leadFieldECGVerifierBase
{
    //- The parameters of the solution
    shellTransmissionMMS::Params shellParams_;

    //- The lead fields z_e
    PtrList<volScalarField> leadFields_;

    //- The solution zone of every base cell (its cellZone)
    labelList shellZones_;

    //- phiE source and outer-boundary current at time t as one source with
    //  the dimensions dims, volume mean removed
    tmp<volScalarField> shellSource
    (
        const scalar t,
        const dimensionSet& dims
    ) const;


    // Reference hooks

        virtual word referenceName() const override;

        //- Check the conductor against the solution and read the lead
        //  fields
        virtual void initialiseReference
        (
            const leadFieldConductor& lfc
        ) override;

        virtual scalar exactVm
        (
            const point& p,
            const scalar t
        ) const override;

        //- The formula of zone, or of the zone containing p when zone < 0
        virtual scalar referencePhi
        (
            const point& p,
            const scalar t,
            const label zone
        ) const override;

        virtual label cellZone(const label cellI) const override
        {
            return shellZones_[cellI];
        }

        virtual tmp<volScalarField> staticSource
        (
            const volScalarField& VmExact,
            const scalar t,
            scalarField& leadFieldStatic
        ) const override;


public:

    OverrideTypeName("manufacturedShellLeadFieldECGVerifier");

    manufacturedShellLeadFieldECGVerifier
    (
        const electroStateProvider& stateProvider,
        const word& domainName,
        const dictionary& dict,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual ~manufacturedShellLeadFieldECGVerifier() = default;

    virtual bool read(const dictionary& dict) override;

    //- Write Gi, Ge and the torso conductivity of the solution on the
    //  undecomposed mesh
    virtual void writeManufacturedFields(const fvMesh& mesh) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
