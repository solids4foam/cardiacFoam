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

#include "eikonalLeadFieldECG.H"
#include "ecgDomain.H"
#include "ecgModelIO.H"
#include "IOdictionary.H"
#include "sealedWallPolicy.H"
#include "gradScheme.H"
#include "Pstream.H"
#include "dimVoltage.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
defineTypeNameAndDebug(eikonalLeadFieldECG, 0);
addToRunTimeSelectionTable(ecgSolver, eikonalLeadFieldECG, dictionary);

eikonalLeadFieldECG::eikonalLeadFieldECG(const dictionary& dict)
:
    templates_(dict),
    start_(dict.subDict("sampling").get<scalar>("start")),
    end_(dict.subDict("sampling").get<scalar>("end")),
    delta_(dict.subDict("sampling").get<scalar>("deltaT")),
    heartOnly_(false), sampled_(false), data_(), output_(), lastValues_()
{
    if (!templates_.homogeneous())
    {
        FatalIOErrorInFunction(dict)
            << "eikonalLeadFieldECG requires a spatially uniform template: "
            << "set templateType endo, mid, epi or manufactured. Regional and "
            << "personalized templates require a separate spatial-gradient model."
            << exit(FatalIOError);
    }
    if (!std::isfinite(start_) || !std::isfinite(end_)
     || !std::isfinite(delta_) || delta_ <= 0 || end_ < start_)
    {
        FatalIOErrorInFunction(dict) << "Invalid eikonalLeadFieldECG sampling interval."
            << exit(FatalIOError);
    }
    const word mode(dict.lookupOrDefault<word>("importMode","originalConductor"));
    if (mode != "originalConductor" && mode != "heartOnly")
    {
        FatalIOErrorInFunction(dict) << "importMode must be originalConductor or heartOnly."
            << exit(FatalIOError);
    }
    heartOnly_=mode=="heartOnly";
}

void eikonalLeadFieldECG::initialise(const ecgDomain& domain)
{
    const IOdictionary props(IOobject("leadFieldProperties",domain.mesh().time().constant(),
        domain.baseMesh(),IOobject::MUST_READ,IOobject::NO_WRITE));
    data_.reset(new leadFieldData(domain,props,heartOnly_));
}

wordList eikonalLeadFieldECG::channelNames(const ecgDomain&) const
{
    return data_->names();
}

List<vector> eikonalLeadFieldECG::channelPositions(const ecgDomain&) const
{
    return List<vector>();
}

void eikonalLeadFieldECG::solve(ecgDomain& domain,scalar,scalar,scalarField& values)
{
    if (sampled_) { values=lastValues_; return; }
    templates_.initialise(domain);
    const fvMesh& mesh=domain.mesh();
    const volScalarField& source=domain.activationTime();
    // Use fresh private gradient storage and copy the source owner's wall policy.
    // The activation field and its cached solver gradient remain untouched.
    volScalarField activation(IOobject("eikonalLeadFieldObservationTau",
        mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),source);
    const labelIOList* states=mesh.findObject<labelIOList>
        (sealedWallPolicy::statesName(source.name()));
    if (!states)
    {
        FatalErrorInFunction << "Activation-time owner must publish its wall-gradient policy."
            << exit(FatalError);
    }
    sealedWallPolicy::publish(mesh,activation.name(),sealedWallPolicy::states(mesh,source.name()));
    activation.correctBoundaryConditions();
    const tmp<volVectorField> gradient=fv::gradScheme<scalar>::New
    (
        mesh,mesh.gradScheme("eikonalLeadFieldGradTau")
    )().grad(activation,"grad("+activation.name()+')');
    volVectorField gradVm(IOobject("eikonalLeadFieldGradVm",mesh.time().timeName(),mesh,
        IOobject::NO_READ,IOobject::NO_WRITE),mesh,
        dimensionedVector("zero",dimVoltage/dimLength,Zero));
    if (Pstream::master())
        output_=ecgModelIO::openTimeSeries(mesh.time().globalPath()/"postProcessing",
            "eikonalLeadFieldECG.dat",data_->names());
    for (label i=0; ; ++i)
    {
        const scalar time=start_+i*delta_;
        if (time>end_+SMALL) break;
        templates_.reconstructGradVm(time,activation,gradient(),gradVm);
        data_->projectGradient(gradVm,lastValues_);
        if (Pstream::master()) ecgModelIO::writeRow(output_.ref(),time,lastValues_);
        domain.recordVerification(time,lastValues_);
    }
    sampled_=true;
    values=lastValues_;
    Info<< "eikonalLeadFieldECG: wrote uniform " << templates_.templateType()
        << " template recordings from " << start_ << " to " << end_
        << " s to postProcessing/eikonalLeadFieldECG.dat" << nl;
}
}
