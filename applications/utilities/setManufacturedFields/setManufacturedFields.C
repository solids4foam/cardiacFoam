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

Application
    setManufacturedFields

Description
    Write the fields a manufactured solution prescribes (conductivities and
    the cell decomposition) for the verification model selected in
    constant/electroProperties: the myocardium's, or, when the myocardium
    has none, the one of the only ECG domain with a verification model.
    Serial only.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "electroVerificationModel.H"
#include "eikonalVerificationModel.H"
#include "ecgVerificationModel.H"
#include "dimVoltage.H"

using namespace Foam;

namespace
{

// The case mesh as the state of an ECG domain whose verification model only
// writes its manufactured fields: no state channels
class meshStateProvider
:
    public electroStateProvider
{
    const fvMesh& mesh_;

public:

    explicit meshStateProvider(const fvMesh& mesh)
    :
        mesh_(mesh)
    {}

    virtual const fvMesh& mesh() const
    {
        return mesh_;
    }
};

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word solverName(electroProperties.get<word>("myocardiumSolver"));
    const dictionary& coeffs =
        electroProperties.subDict(solverName + "Coeffs");

    const word modelType(electroVerificationModel::selectedType(coeffs));

    // The ECG domains with a verification model
    wordList verifiedEcgDomains;
    if (const dictionary* ecgPtr = electroProperties.findDict("ecgDomains"))
    {
        for (const entry& domain : *ecgPtr)
        {
            if (!domain.isDict())
            {
                continue;
            }
            const word ecgType
            (
                ecgVerificationModel::selectedType(domain.dict())
            );
            if (!ecgType.empty() && ecgType != "none")
            {
                verifiedEcgDomains.append(domain.keyword());
            }
        }
    }

    if
    (
        (modelType.empty() || modelType == "none")
     && !verifiedEcgDomains.empty()
    )
    {
        // The myocardium has no manufactured solution: the ECG domain's
        // verification model writes the fields
        if (verifiedEcgDomains.size() > 1)
        {
            FatalIOErrorInFunction(electroProperties)
                << "The ECG domains " << verifiedEcgDomains << " each have "
                << "a verification model and the myocardium has none: the "
                << "fields of one would overwrite those of the others."
                << exit(FatalIOError);
        }

        const meshStateProvider provider(mesh);
        const word& domainName = verifiedEcgDomains.first();

        // The fields do not depend on the electrodes
        ecgVerificationModel::New
        (
            provider,
            domainName,
            electroProperties.subDict("ecgDomains").subDict(domainName),
            wordList(),
            List<vector>()
        )->writeManufacturedFields(mesh);
    }
    else if (eikonalVerificationModel::dictionaryConstructorTable(modelType))
    {
        // The verifier reads the conductivity only once the solver uses it
        const volTensorField conductivity
        (
            IOobject
            (
                "conductivity",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor(dimless, Zero)
        );

        eikonalVerificationModel::New
        (
            coeffs,
            mesh,
            conductivity,
            dimensionedScalar("chi", dimArea/dimVolume, coeffs),
            dimensionedScalar
            (
                "cm",
                dimCurrent*dimTime/(dimVoltage*dimArea),
                coeffs
            ),
            dimensionedScalar("c0", coeffs),
            Switch(coeffs.lookup("eikonalAdvectionDiffusionApproach"))
        )->writeManufacturedFields(mesh);
    }
    else
    {
        autoPtr<electroVerificationModel> verifier
        (
            electroVerificationModel::New(coeffs)
        );

        if (!verifier)
        {
            FatalIOErrorInFunction(electroProperties)
                << "No verification model is selected: "
                << solverName << "Coeffs needs a verificationModel dictionary "
                << "whose type names a manufactured-solution verifier."
                << exit(FatalIOError);
        }

        verifier->writeManufacturedFields(mesh);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
