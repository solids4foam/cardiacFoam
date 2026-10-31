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
    constant/electroProperties. Serial only.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IOdictionary.H"
#include "electroVerificationModel.H"
#include "eikonalVerificationModel.H"
#include "dimVoltage.H"
#include "ecgVerification/manufacturedLeadFieldShellFields.H"

using namespace Foam;

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

    // An ECG domain verified against the sphericalShell lead-field reference
    // takes the conductivities of that solution
    bool shellLeadField = false;
    if (const dictionary* ecgPtr = coeffs.findDict("ecgDomains"))
    {
        for (const entry& domain : *ecgPtr)
        {
            const dictionary* verifierPtr =
                domain.isDict()
              ? domain.dict().findDict("verificationModel")
              : nullptr;
            if
            (
                verifierPtr
             && verifierPtr->getOrDefault<word>("type", "")
             == "manufacturedLeadFieldECGVerifier"
             && verifierPtr->getOrDefault<word>("reference", "")
             == "sphericalShell"
            )
            {
                shellLeadField = true;
            }
        }
    }
    if (shellLeadField)
    {
        writeManufacturedLeadFieldShellFields(mesh);
    }

    if (shellLeadField && (modelType.empty() || modelType == "none"))
    {
        // The myocardium has no manufactured solution of its own
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
