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

#ifndef myocardiumDomainInterface_H
#define myocardiumDomainInterface_H

#include "electroDomainInterface.H"
#include "electroStateProvider.H"
#include "electromechanicalSignalProvider.H"
#include "electroDomainCouplingEndpoints.H"
#include "autoPtr.H"
#include "pimpleControl.H"

namespace Foam
{

class fvMesh;
class dictionary;
class ionicModel;
class electroVerificationModel;

class myocardiumDomainInterface
:
    public electroDomainInterface,
    public tissueCouplingEndpoint,
    public electroStateProvider
{
public:

    using electroDomainInterface::advance;

    static autoPtr<myocardiumDomainInterface> New
    (
        const fvMesh& mesh,
        const dictionary& electroProperties,
        PtrList<volScalarField>& outFields,
        const wordList& postProcessFieldNames,
        PtrList<volScalarField>& postProcessFields,
        autoPtr<ionicModel>& ionicModelPtr,
        autoPtr<electroVerificationModel>& verificationModelPtr,
        scalar initialDeltaT
    );

    virtual ~myocardiumDomainInterface() = default;

    virtual const ElectromechanicalSignalProvider* provider() const
    {
        return nullptr;
    }

    virtual const volScalarField& activationTime() const = 0;

    virtual void advance
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr
    ) = 0;

    virtual scalar suggestExplicitDeltaT(scalar maxCo) const = 0;

    virtual bool shouldPostProcess() const = 0;

    virtual void exportStates() = 0;

    virtual void postProcess() = 0;

    virtual bool applyModelTimeControls(Time& runTime) const
    {
        (void)runTime;
        return false;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
