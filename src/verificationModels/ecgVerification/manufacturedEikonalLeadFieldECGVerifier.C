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

#include "manufacturedEikonalLeadFieldECGVerifier.H"
#include "eikonalVerification/manufacturedEikonalReference.H"
#include "IOdictionary.H"
#include "ecgModelIO.H"
#include "Pstream.H"
#include "addToRunTimeSelectionTable.H"
#include <cmath>

namespace Foam
{
defineTypeNameAndDebug(manufacturedEikonalLeadFieldECGVerifier,0);
addToRunTimeSelectionTable(ecgVerificationModel,manufacturedEikonalLeadFieldECGVerifier,dictionary);

manufacturedEikonalLeadFieldECGVerifier::manufacturedEikonalLeadFieldECGVerifier
(
    const electroStateProvider& state,const word& name,const dictionary& dict,
    const wordList& channels,const List<vector>& positions
)
:
    ecgVerificationModel(state,name,channels,positions),
    scales_(),nodes_(),weights_(),checkNodes_(),checkWeights_(),
    maxError_(0),maxDelta_(0),samples_(0),output_()
{
    read(dict);
    if (Pstream::master())
    {
        wordList columns(5*channels.size());
        forAll(channels,i)
        {
            columns[5*i]=channels[i]+"_numeric";
            columns[5*i+1]=channels[i]+"_exact";
            columns[5*i+2]=channels[i]+"_error";
            columns[5*i+3]=channels[i]+"_check";
            columns[5*i+4]=channels[i]+"_quadratureDelta";
        }
        output_=ecgModelIO::openTimeSeries(mesh_.time().globalPath()/"postProcessing",
            "manufacturedEikonalLeadFieldECG_"+name+".dat",columns);
    }
}

ecgVerificationModel::Requirements manufacturedEikonalLeadFieldECGVerifier::requirements() const
{
    Requirements r; r.needActivationTime=true; return r;
}

bool manufacturedEikonalLeadFieldECGVerifier::read(const dictionary& dict)
{
    const dictionary& cfg=dict.subDict("verificationModel");
    if (mesh_.nGeometricD()!=3)
        FatalIOErrorInFunction(cfg) << "Eikonal lead-field MMS requires the 3D unit cube."
            << exit(FatalIOError);
    validateManufacturedEikonalUnitDomain(mesh_,3);
    scales_=scalarField(cfg.lookup("scales"));
    if (scales_.size()!=electrodeNames_.size())
        FatalIOErrorInFunction(cfg) << "scales must match lead-field channels." << exit(FatalIOError);
    const label q=cfg.lookupOrDefault<label>("referenceQuadratureOrder",96);
    const label check=cfg.lookupOrDefault<label>("checkQuadratureOrder",64);
    if (check<=0 || q<=check)
        FatalIOErrorInFunction(cfg) << "Require 0 < checkQuadratureOrder < referenceQuadratureOrder."
            << exit(FatalIOError);
    pseudoECGManufacturedQuadratureRule(q,nodes_,weights_);
    pseudoECGManufacturedQuadratureRule(check,checkNodes_,checkWeights_);
    // Check the prescribed sensitivity contract independently of the observer.
    const IOdictionary props(IOobject("leadFieldProperties",mesh_.time().constant(),mesh_,
        IOobject::MUST_READ,IOobject::NO_WRITE,false));
    forAll(electrodeNames_,i)
    {
        if (!std::isfinite(scales_[i]))
            FatalIOErrorInFunction(cfg) << "scales must be finite." << exit(FatalIOError);
        const volVectorField field(IOobject(props.subDict("leadFields").get<word>(electrodeNames_[i]),
            props.lookupOrDefault<word>("fieldInstance","0"),mesh_,
            IOobject::MUST_READ,IOobject::NO_WRITE,false),mesh_);
        scalar mismatch=0;
        forAll(field,c) mismatch=max(mismatch,mag(field[c]-scales_[i]*vector(1,2,3)));
        reduce(mismatch,maxOp<scalar>());
        if (field.dimensions()!=dimensionSet(0,-2,0,0,0,0,0) || mismatch>1e-12)
            FatalIOErrorInFunction(cfg) << "Reference requires L_i=scale_i*(1,2,3) m^-2."
                << exit(FatalIOError);
    }
    return true;
}

scalar manufacturedEikonalLeadFieldECGVerifier::reference
(
    scalar time,const UList<scalar>& nodes,const UList<scalar>& weights
) const
{
    // Divergence theorem: constant L turns the volume integral L.grad(Vm)
    // into face voltage differences. The six faces are symmetric in x,y,z.
    // This oracle does not use a numerical or template-reconstructed gradient.
    const scalar pi2=2*constant::mathematical::pi;
    scalar integral=0;
    forAll(nodes,i) forAll(nodes,j)
    {
        const scalar a=0.5*(nodes[i]+1)+0.5*(nodes[j]+1);
        integral+=0.25*weights[i]*weights[j]
            *(Foam::sin(pi2*(time-Foam::exp(1+a)))
             -Foam::sin(pi2*(time-Foam::exp(a))));
    }
    return 6*integral;
}

void manufacturedEikonalLeadFieldECGVerifier::record(const List<scalar>& values)
{
    record(mesh_.time().value(),values);
}

void manufacturedEikonalLeadFieldECGVerifier::record(scalar time,const List<scalar>& values)
{
    if (values.size()!=scales_.size())
        FatalErrorInFunction << "Lead-field MMS channel count changed." << exit(FatalError);
    const scalar exact=reference(time,nodes_,weights_);
    const scalar check=reference(time,checkNodes_,checkWeights_);
    scalarField row(5*values.size());
    forAll(values,i)
    {
        const scalar error=mag(values[i]-scales_[i]*exact);
        const scalar delta=mag(scales_[i]*(exact-check));
        maxError_=max(maxError_,error); maxDelta_=max(maxDelta_,delta);
        row[5*i]=values[i]; row[5*i+1]=scales_[i]*exact; row[5*i+2]=error;
        row[5*i+3]=scales_[i]*check; row[5*i+4]=delta;
    }
    ++samples_;
    if (Pstream::master()) ecgModelIO::writeRow(output_.ref(),time,row);
}

void manufacturedEikonalLeadFieldECGVerifier::end()
{
    if (Pstream::master())
    {
        OFstream out(mesh_.time().globalPath()/"postProcessing"/
            ("manufacturedEikonalLeadFieldECGSummary_"+domainName_+".dat"));
        out.precision(12);
        out << "# Manufactured eikonal lead-field ECG summary" << nl
            << "samples " << samples_ << nl << "maxError " << maxError_ << nl
            << "maxQuadratureDelta " << maxDelta_ << nl;
    }
}
}
