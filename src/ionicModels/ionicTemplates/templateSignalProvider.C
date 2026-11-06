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

#include "templateSignalProvider.H"
#include "IOdictionary.H"

namespace Foam
{

namespace
{

const volScalarField& lookupActivationTime(const fvMesh& mesh)
{
    if (!mesh.foundObject<volScalarField>("activationTime"))
    {
        FatalErrorInFunction
            << "Template signals need an activationTime field on mesh "
            << mesh.name() << "." << exit(FatalError);
    }
    return mesh.lookupObject<volScalarField>("activationTime");
}

} // End unnamed namespace


templateSignalProvider::templateSignalProvider
(
    const fvMesh& mesh,
    const dictionary& templatesDict
)
:
    activationTime_(lookupActivationTime(mesh)),
    templates_(),
    weights_()
{
    const IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary* hetDictPtr =
        ionicTemplates::findHeterogeneityDict(electroProperties);

    if (!hetDictPtr)
    {
        FatalIOErrorInFunction(templatesDict)
            << "Template signals need an ionicHeterogeneity block in "
            << "electroProperties or its myocardium solver Coeffs."
            << exit(FatalIOError);
    }

    templates_ = ionicTemplates::generate(templatesDict, *hetDictPtr, mesh.time());
    weights_ = ionicTemplates::regionWeights(mesh, *hetDictPtr);

    if (weights_.size() != templates_.size())
    {
        FatalErrorInFunction
            << weights_.size() << " region weights for " << templates_.size()
            << " templates." << exit(FatalError);
    }
}


scalar templateSignalProvider::signal
(
    const label i,
    const CouplingSignal s
) const
{
    const scalar tau = activationTime_[i];
    const scalar localTime =
        tau < 0 ? -GREAT : activationTime_.time().value() - tau;

    scalar value = 0;
    forAll(templates_, r)
    {
        const scalar w = weights_[r][i];
        if (w == 0)
        {
            continue;
        }
        const ionicTemplates::DynamicTemplate& tpl = templates_[r];
        const scalarField& values =
            s == CouplingSignal::CAI ? tpl.valuesCai : tpl.valuesMv;
        value +=
            w
           *ionicTemplates::evaluateTemplate
            (
                localTime,
                tpl.times.cdata(),
                values.cdata(),
                values.size()
            );
    }
    return value;
}

} // End namespace Foam

// ************************************************************************* //
