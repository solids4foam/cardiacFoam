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
    buildLeadFields

Description
    Lead fields of a leadFieldECG domain: for each electrode, the reciprocal
    solve A z = delta with the forward heart+torso operator A of
    passiveConductor and a unit current delta in the electrode cell.

    The unknown is named phiE on the base and heart meshes, so the forward's
    fvSchemes entries and wall-gradient states apply. Ground patches are
    fixedValue 0, other non-coupled patches zeroGradient; without them the
    reference cell is pinned to 0. Assembly, reference and solve repeat
    until the relative change of z between passes is below 1e-10: with the
    fvSolution entry leadField until the change first falls below it, then
    with leadFieldFinal (relTol 0, tight tolerance) until it falls below it
    again.

    Every electrode lies in a torso cell; an electrode in a myocardium cell
    is fatal (leadFieldConductor).

    Writes constant/[region/]leadField_<electrode> and leadFieldSignature.
    With -parallel each processor writes its own processorN/constant, so a
    parallel ECG run needs lead fields built on its decomposition (run after
    decomposePar) or decomposed from a serial build.

Usage
    buildLeadFields [-region <name>] [-ecgDomain <name>]

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "leadFieldConductor.H"
#include "sealedWallPolicy.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include <limits>

using namespace Foam;

namespace
{

// IOobject that is neither read nor written automatically.
IOobject fieldIO
(
    const word& name, const word& instance, const fvMesh& mesh, bool reg
)
{
    return IOobject
    (
        name, instance, mesh, IOobject::NO_READ, IOobject::NO_WRITE, reg
    );
}


// The requested leadFieldECG entry of ecgDomains, or the only one.
word selectEcgDomain(const dictionary& ecgDomains, const word& requested)
{
    const word lfType("leadFieldECG");
    wordList names;
    for (const entry& e : ecgDomains)
    {
        const dictionary* d = e.dictPtr();
        if (d && d->getOrDefault<word>("ecgSolver", "pseudoECG") == lfType)
        {
            names.append(e.keyword());
        }
    }

    if (requested.empty() ? names.size() != 1 : !names.found(requested))
    {
        FatalIOErrorInFunction(ecgDomains)
            << "Need exactly one ecgDomains entry with ecgSolver "
            << "leadFieldECG, or -ecgDomain naming one of them. Entries: "
            << names << ", -ecgDomain: '" << requested << "'."
            << exit(FatalIOError);
    }

    return requested.empty() ? names[0] : requested;
}


} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote("Write the lead fields of a leadFieldECG domain.");
    #include "addRegionOption.H"
    argList::addOption
    (
        "ecgDomain", "name",
        "leadFieldECG entry of ecgDomains (default: the only one)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Lead fields and signature keep every digit
    IOstream::defaultPrecision(std::numeric_limits<scalar>::max_digits10);

    // constant/[region/]electroProperties of the undecomposed case
    const IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties", runTime.caseConstant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE, false
        )
    );
    const dictionary& electroCoeffs = electroProperties.subDict
    (
        electroProperties.get<word>("myocardiumSolver") + "Coeffs"
    );
    const word ecgDomainName = selectEcgDomain
    (
        electroCoeffs.subDict("ecgDomains"),
        args.getOrDefault<word>("ecgDomain", word::null)
    );

    const leadFieldConductor lfc(mesh, electroCoeffs, ecgDomainName);
    const passiveConductor& conductor = lfc.conductor();
    const fvMesh& heart = lfc.heartMesh();
    const labelList& heartToBase = conductor.heartCellToBaseCell();
    const dictionary& solverDict = mesh.solverDict("leadField");
    const dictionary& finalSolverDict = mesh.solverDict("leadFieldFinal");
    const scalar convergenceTol = 1e-10;
    const label maxPasses = 100;

    // z is in 1/(conductivity length), so that A z = delta in 1/volume.
    const dimensionedScalar zZero
    (
        "zero", dimless/(conductor.sigmaTotal().dimensions()*dimLength), 0.0
    );

    wordList patchTypes
    (
        mesh.boundary().size(), zeroGradientFvPatchScalarField::typeName
    );
    for (const word& patchName : lfc.groundPatches())
    {
        patchTypes[mesh.boundaryMesh().findPatchID(patchName)] =
            fixedValueFvPatchScalarField::typeName;
    }

    volScalarField z
    (
        fieldIO("phiE", runTime.timeName(), mesh, true), mesh, zZero,
        patchTypes
    );
    sealedWallPolicy::publish(z, false);

    volScalarField heartZ
    (
        fieldIO("phiE", runTime.timeName(), heart, true), heart, zZero,
        "zeroGradient"
    );
    sealedWallPolicy::publish(heartZ, conductor.sealedHeartBoundary());

    volScalarField delta
    (
        fieldIO("leadFieldSource", runTime.timeName(), mesh, false), mesh,
        dimensionedScalar("zero", dimless/dimVolume, 0.0)
    );

    const wordList& names = lfc.electrodeNames();
    const labelList& cells = lfc.electrodeCells();

    forAll(names, electrodeI)
    {
        Info<< "Lead field " << names[electrodeI] << nl;

        z == zZero;
        delta.primitiveFieldRef() = 0.0;
        if (cells[electrodeI] >= 0)
        {
            delta.primitiveFieldRef()[cells[electrodeI]] =
                1.0/mesh.V()[cells[electrodeI]];
        }

        const scalarField& zI = z.primitiveField();
        scalarField zPrev(zI.size());
        scalar change = GREAT;
        label pass = 0;
        label finalPasses = 0;

        // Fixed point of the explicit corrections, started from z = 0.
        // Once the change is below the tolerance, the passes solve with
        // leadFieldFinal until it is below the tolerance again.
        while (finalPasses == 0 || change >= convergenceTol)
        {
            if (pass == maxPasses)
            {
                FatalErrorInFunction
                    << "Lead field of electrode " << names[electrodeI]
                    << " not converged after " << pass
                    << " passes; last relative change " << change
                    << " (tolerance " << convergenceTol << ")."
                    << exit(FatalError);
            }

            zPrev = zI;

            // The heart companion takes the base values before each
            // assembly, as the forward's heart phiE does.
            scalarField& heartZI = heartZ.primitiveFieldRef();
            forAll(heartToBase, heartCellI)
            {
                heartZI[heartCellI] = z[heartToBase[heartCellI]];
            }
            heartZ.correctBoundaryConditions();

            tmp<fvScalarMatrix> tEqn(conductor.assemble(z, heartZ, &delta));
            conductor.applyReference(tEqn.ref(), z, lfc.referenceCell(), 0.0);

            const bool finalPass =
                finalPasses > 0 || (pass >= 2 && change < convergenceTol);
            if (finalPass)
            {
                ++finalPasses;
            }
            solve(tEqn.ref(), finalPass ? finalSolverDict : solverDict);

            ++pass;
            change = gMax(mag(zI - zPrev))/max(gMax(mag(zI)), VSMALL);
        }

        Info<< "  " << pass << " passes (" << finalPasses
            << " with leadFieldFinal), relative change " << change << nl;

        const word fieldName("leadField_" + names[electrodeI]);
        volScalarField
        (
            fieldIO(fieldName, runTime.constant(), mesh, false), z
        ).write();
    }

    IOdictionary
    (
        fieldIO("leadFieldSignature", runTime.constant(), mesh, false),
        lfc.signature()
    ).regIOobject::write();

    Info<< "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
