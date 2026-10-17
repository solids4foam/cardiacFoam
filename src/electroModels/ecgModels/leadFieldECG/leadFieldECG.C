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
#include "IOdictionary.H"
#include "addToRunTimeSelectionTable.H"
namespace Foam
{
defineTypeNameAndDebug(leadFieldECG, 0);
addToRunTimeSelectionTable(ecgSolver, leadFieldECG, dictionary);
leadFieldECG::leadFieldECG(const dictionary& dict) : data_(), heartOnly_(false)
{
    const word mode(dict.lookupOrDefault<word>("importMode","originalConductor"));
    if (mode!="originalConductor" && mode!="heartOnly")
        FatalIOErrorInFunction(dict) << "importMode must be originalConductor or heartOnly." << exit(FatalIOError);
    heartOnly_=mode=="heartOnly";
}
void leadFieldECG::initialise(const ecgDomain& domain)
{
    const IOdictionary props(IOobject("leadFieldProperties", domain.mesh().time().constant(),
        domain.baseMesh(), IOobject::MUST_READ, IOobject::NO_WRITE));
    data_.reset(new leadFieldData(domain, props, heartOnly_));
}
wordList leadFieldECG::channelNames(const ecgDomain&) const { return data_->names(); }
List<vector> leadFieldECG::channelPositions(const ecgDomain&) const { return List<vector>(); }
void leadFieldECG::solve(ecgDomain& domain, scalar, scalar, scalarField& values)
{
    data_->project(domain, values);
}
}
