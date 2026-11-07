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
    manufacturedLeadFieldECGVerifier

Description
    Lead-field ECG verifier against the slab and layered references
    (manufacturedLeadFieldECGReference.H), with the channels of
    leadFieldECGVerifierBase.

    Keys of verificationModel:
        reference   slab;    // required: slab or layered
        layerInterfaces (-0.5 1.5);  // layered only, required: xL < 0, xR > 1

    slab needs an insulated conductor: no ground or surface-current
    patches. layered is pinned (no ground patches); it takes alpha from the
    surface-current patches of the conductor (faces normal to x, current
    times n_x equal to alpha on every face, no net current; none for
    alpha = 0) and the layer conductivities from its torso cells, which must
    each lie in one layer with one isotropic value per layer. Both prescribe
    no fields: their conductivities are those of the case.

SourceFiles
    manufacturedLeadFieldECGVerifier.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedLeadFieldECGVerifier_H
#define manufacturedLeadFieldECGVerifier_H

#include "leadFieldECGVerifierBase.H"
#include "ecgVerification/manufacturedLeadFieldECGReference.H"
#include "Enum.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedLeadFieldECGVerifier
:
    public leadFieldECGVerifierBase
{
public:

    //- Exact solutions
    enum class referenceType
    {
        slab,
        layered
    };

    static const Enum<referenceType> referenceTypeNames;


private:

    referenceType reference_;

    //- Gi_xx/(Gi_xx + Ge_xx)
    scalar slabK_;

    //- Layered: the layers and the current density -alpha along x
    manufacturedLeadFieldLayers layers_;
    scalar alpha_;

    //- Layered: alpha from the surface currents, the layer conductivities
    //  from the torso cells
    void initialiseLayered(const leadFieldConductor& lfc);


    // Reference hooks

        virtual word referenceName() const override;

        virtual void initialiseReference
        (
            const leadFieldConductor& lfc
        ) override;

        virtual scalar exactVm
        (
            const point& p,
            const scalar t
        ) const override;

        virtual scalar referencePhi
        (
            const point& p,
            const scalar t,
            const label zone
        ) const override;

        virtual void reportParameters(Ostream& os) const override;

        virtual void writeParameters(Ostream& os) const override;


public:

    OverrideTypeName("manufacturedLeadFieldECGVerifier");

    manufacturedLeadFieldECGVerifier
    (
        const electroStateProvider& stateProvider,
        const word& domainName,
        const dictionary& dict,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual ~manufacturedLeadFieldECGVerifier() = default;

    //- Read the verifier entry (called by the constructor)
    bool read(const dictionary& dict);
};

} // End namespace Foam

#endif

// ************************************************************************* //
