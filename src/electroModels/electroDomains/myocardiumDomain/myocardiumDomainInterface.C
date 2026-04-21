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

#include "myocardiumDomainInterface.H"
#include "myocardiumDomain.H"
#include "eikonalMyocardiumDomain.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
#include "error.H"

namespace Foam
{

namespace
{

word myocardiumSolverType(const dictionary& electroProperties)
{
    const word coeffDictName = electroProperties.dictName();

    return coeffDictName.endsWith("Coeffs")
        ? word(coeffDictName.substr(0, coeffDictName.size() - 6))
        : coeffDictName;
}

} // End anonymous namespace


autoPtr<myocardiumDomainInterface> myocardiumDomainInterface::New
(
    const fvMesh& mesh,
    const dictionary& electroProperties,
    PtrList<volScalarField>& outFields,
    const wordList& postProcessFieldNames,
    PtrList<volScalarField>& postProcessFields,
    autoPtr<ionicModel>& ionicModelPtr,
    autoPtr<electroVerificationModel>& verificationModelPtr,
    scalar initialDeltaT
)
{
    const word solverType = myocardiumSolverType(electroProperties);

    if (solverType == "eikonalSolver")
    {
        ionicModelPtr.clear();
        verificationModelPtr.clear();

        return autoPtr<myocardiumDomainInterface>
        (
            new EikonalMyocardiumDomain(mesh, electroProperties)
        );
    }

    ionicModelPtr =
        ionicModel::New
        (
            electroProperties,
            MyocardiumDomain::configuredCellCount(mesh, electroProperties),
            initialDeltaT
        );

    verificationModelPtr =
        electroVerificationModel::New(electroProperties);

    return autoPtr<myocardiumDomainInterface>
    (
        MyocardiumDomain::New
        (
            mesh,
            electroProperties,
            outFields,
            postProcessFieldNames,
            postProcessFields,
            ionicModelPtr(),
            verificationModelPtr.get()
        ).ptr()
    );
}

} // End namespace Foam

// ************************************************************************* //
