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
#include <limits>

namespace Foam
{

leadFieldECG::leadFieldECG
(
    const dictionary& electroProperties,
    const word& ecgDomainName
)
:
    electroProperties_(electroProperties),
    ecgDomainName_(ecgDomainName),
    conductorPtr_(),
    heartWeights_(),
    constants_()
{}


void leadFieldECG::load
(
    const electroStateProvider& provider,
    const wordList& electrodeNames,
    const List<vector>& electrodePositions
)
{
    conductorPtr_.reset
    (
        new leadFieldConductor
        (
            provider,
            electroProperties_,
            ecgDomainName_,
            electrodeNames,
            electrodePositions
        )
    );
    const leadFieldConductor& lfc = conductorPtr_();
    const passiveConductor::torsoSettings& torso = lfc.conductor().torso();
    const wordList& names = lfc.electrodeNames();
    const fvMesh& mesh = provider.baseMesh();

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
        max(mag(torso.groundValues)) > 0
     || max(mag(torso.surfaceCurrentValues)) > 0
     || (torso.pinsReference() && torso.phiEReferenceValue != 0);
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

    Info<< "torsoECG lead fields: " << names.size()
        << " lead fields, constant terms " << constants_ << endl;
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

} // End namespace Foam

// ************************************************************************* //
