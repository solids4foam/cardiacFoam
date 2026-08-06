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
    manufacturedFDABathBidomainVerifier

Description
    Manufactured FDA bath bidomain verifier.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedFDABathBidomainVerifier_H
#define manufacturedFDABathBidomainVerifier_H

#include "electroVerificationModel.H"
#include "overrideTypeName.H"
#include "Enum.H"

namespace Foam
{

class manufacturedFDABathBidomainVerifier
:
    public electroVerificationModel
{
    volScalarField* phiEPtr_;
    const labelUList* phiEHeartCellMapPtr_;
    Switch enabled_;
    bool useExplicitAlgorithm_;
    mutable bool errorsReported_;
    scalar k_;
    scalar alpha_;

public:

    //- Which of the two FDA bidomain-with-bath boundary variants this case
    //  runs. The choice drives both the boundary conditions and the error
    //  metric, which is why it lives here rather than being set by hand in
    //  0/phiE: with an electrode pair the exact phiE is only defined up to an
    //  arbitrary C(t), so its error must be taken after a gauge shift, and a
    //  case configured one way but measured the other would report silently
    //  wrong norms.
    enum class bathVariant
    {
        groundElectrode,   //!< Dirichlet phiE = 0 at x = -1, I_E = +alpha at x = 2
        electrodePair      //!< I_E = -alpha at x = -1, I_E = +alpha at x = 2
    };

    static const Enum<bathVariant> bathVariantNames;

private:

    bathVariant variant_;
    fileName outputFileName_;

    //- Mean of a field over the mesh, volume weighted. Used to remove the
    //  arbitrary constant before comparing phiE in the electrodePair variant.
    static scalar volumeMean(const fvMesh& m, const scalarField& f);

    const dictionary& verificationDict() const;

public:

    OverrideTypeName("manufacturedFDABathBidomainVerifier");

    explicit manufacturedFDABathBidomainVerifier(const dictionary& dict);

    virtual void bindBidomainField(volScalarField& phiE) override;

    virtual void bindBidomainField
    (
        volScalarField& phiE,
        const labelUList& heartCellMap
    ) override;

    virtual void unbindBidomainField() override;

    virtual wordList preProcessFieldNames(const ionicModel& model) const override;

    virtual wordList requiredPostProcessFieldNames
    (const ionicModel& model) const override;

    virtual bool shouldPostProcess
    (const ionicModel& model, const volScalarField& Vm) const override;

    virtual void preProcess
    (
        ionicModel& model,
        volScalarField& Vm,
        PtrList<volScalarField>& fields
    ) override;

    virtual void postProcess
    (
        const ionicModel& model,
        const volScalarField& Vm,
        const PtrList<volScalarField>& fields
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
