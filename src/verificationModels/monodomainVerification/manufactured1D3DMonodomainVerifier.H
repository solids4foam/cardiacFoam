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
    manufactured1D3DMonodomainVerifier

Description
    Diagnostic verifier for the staged 1D-3D manufactured monodomain workflow.

    This verifier is intentionally non-intrusive:
    - it does not overwrite Vm during preProcess,
    - it binds the live myocardium sourceField so it can inspect the summed
      source seen by the PDE solve,
    - it writes coupled-source diagnostics at final time.

Author
    OpenAI. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufactured1D3DMonodomainVerifier_H
#define manufactured1D3DMonodomainVerifier_H

#include "electroVerificationModel.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufactured1D3DMonodomainVerifier : public electroVerificationModel
{
    const volScalarField* sourceFieldPtr_;
    bool reportBaseline3DReference_;
    mutable bool diagnosticsWritten_;

public:

    OverrideTypeName("manufactured1D3DMonodomainVerifier");

    explicit manufactured1D3DMonodomainVerifier(const dictionary& dict);

    virtual void bindSourceField(const volScalarField& sourceField) override;

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
