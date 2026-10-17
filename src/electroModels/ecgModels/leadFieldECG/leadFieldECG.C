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

#include "leadFieldECG.H"
#include "ecgDomain.H"
#include "addToRunTimeSelectionTable.H"
#include <limits>

namespace Foam
{

defineTypeNameAndDebug(leadFieldECG, 0);
addToRunTimeSelectionTable(ecgSolver, leadFieldECG, dictionary);


leadFieldECG::leadFieldECG(const dictionary& dict)
:
    leadFieldECG(dict.parent().parent(), dict.dictName())
{}


leadFieldECG::leadFieldECG
(
    const dictionary& electroCoeffs,
    const word& ecgDomainName
)
:
    electroCoeffs_(electroCoeffs),
    ecgDomainName_(ecgDomainName),
    conductorPtr_(),
    heartWeights_(),
    constants_()
{}


void leadFieldECG::load
(
    const electroStateProvider& provider,
    const wordList& electrodeNames
)
{
    if (electroCoeffs_.dictName() == "eikonalSolverCoeffs")
    {
        FatalErrorInFunction
            << "leadFieldECG with eikonalSolver requires Task 9"
            << exit(FatalError);
    }

    const fvMeshSubset* subsetPtr = provider.meshSubsetPtr();
    if (!subsetPtr)
    {
        FatalIOErrorInFunction(electroCoeffs_)
            << "ECG domain '" << ecgDomainName_ << "' needs the myocardium on "
            << "a cellZone of the heart+torso mesh; set cellZone in "
            << electroCoeffs_.name() << "."
            << exit(FatalIOError);
    }

    const volScalarField* VmPtr = provider.VmPtr();
    if (!VmPtr)
    {
        FatalErrorInFunction
            << "ECG domain '" << ecgDomainName_ << "' needs the myocardium Vm."
            << exit(FatalError);
    }

    conductorPtr_.reset
    (
        new leadFieldConductor
        (
            *subsetPtr, *VmPtr, electroCoeffs_, ecgDomainName_
        )
    );
    const leadFieldConductor& lfc = conductorPtr_();
    const wordList& names = lfc.electrodeNames();
    const fvMesh& mesh = provider.baseMesh();

    if (names != electrodeNames)
    {
        FatalErrorInFunction
            << "Electrodes " << names << " of the lead-field conductor differ "
            << "from the ECG domain's " << electrodeNames << "."
            << exit(FatalError);
    }

    // Files in constant/[region/] of this (processor) case
    const auto io = [&mesh](const word& name)
    {
        return IOobject
        (
            name, mesh.time().constant(), mesh,
            IOobject::MUST_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER
        );
    };

    wordList files(1, word("leadFieldSignature"));
    for (const word& name : names)
    {
        files.append("leadField_" + name);
    }
    // The lead fields are volScalarFields, looked up in this processor's
    // constant; the signature is a global dictionary
    const auto hasFile = [&io](const word& file)
    {
        return
            file == "leadFieldSignature"
          ? io(file).typeHeaderOk<IOdictionary>(false)
          : io(file).typeHeaderOk<volScalarField>(false);
    };
    for (const word& file : files)
    {
        if (!hasFile(file))
        {
            FatalErrorInFunction
                << "ECG domain '" << ecgDomainName_ << "' has no lead-field "
                << "file " << io(file).objectRelPath() << "." << nl
                << "run buildLeadFields"
                << exit(FatalError);
        }
    }

    OStringStream diffs;
    diffs.precision(std::numeric_limits<scalar>::digits10);
    if (!lfc.matches(IOdictionary(io("leadFieldSignature")), diffs))
    {
        FatalErrorInFunction
            << "The lead fields of ECG domain '" << ecgDomainName_ << "' do "
            << "not match this case:" << nl << diffs.str().c_str()
            << "run buildLeadFields"
            << exit(FatalError);
    }

    const bool hasConstant =
        max(mag(lfc.groundValues())) > 0
     || max(mag(lfc.surfaceCurrentValues())) > 0
     || (lfc.groundPatches().empty() && lfc.referenceValue() != 0);
    const tmp<scalarField> tb
    (
        hasConstant ? lfc.constantSource() : tmp<scalarField>()
    );

    const labelList& heartToBase = lfc.conductor().heartCellToBaseCell();
    heartWeights_.setSize(names.size());
    constants_.setSize(names.size(), 0.0);

    forAll(names, e)
    {
        const volScalarField z(io("leadField_" + names[e]), mesh);

        heartWeights_[e].setSize(heartToBase.size());
        forAll(heartToBase, c)
        {
            heartWeights_[e][c] = z[heartToBase[c]]*mesh.V()[heartToBase[c]];
        }

        if (hasConstant)
        {
            constants_[e] = gSum(z.primitiveField()*tb());
        }
    }

    Info<< "leadFieldECG: " << names.size() << " lead fields, constant terms "
        << constants_ << endl;
}


void leadFieldECG::evaluate
(
    const volScalarField& Vm,
    scalarField& values
) const
{
    const passiveConductor& conductor = conductorPtr_->conductor();
    const labelList& heartToBase = conductor.heartCellToBaseCell();
    const tmp<volScalarField> tb(conductor.heartSource(Vm));

    scalarField bHeart(heartToBase.size());
    forAll(heartToBase, c)
    {
        bHeart[c] = tb()[heartToBase[c]];
    }

    values.setSize(heartWeights_.size());
    forAll(values, e)
    {
        values[e] = gSum(heartWeights_[e]*bHeart) + constants_[e];
    }
}


void leadFieldECG::solve
(
    ecgDomain& domain,
    scalar t0,
    scalar dt,
    scalarField& values
)
{
    (void)t0;
    (void)dt;

    if (!loaded())
    {
        load(domain, domain.electrodeNames());
    }

    evaluate(domain.Vm(), values);
}

} // End namespace Foam

// ************************************************************************* //
