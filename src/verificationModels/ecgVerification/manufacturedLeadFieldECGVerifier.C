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

#include "manufacturedLeadFieldECGVerifier.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "ecgModelIO.H"
#include "Pstream.H"
#include "addToRunTimeSelectionTable.H"
namespace Foam
{
defineTypeNameAndDebug(manufacturedLeadFieldECGVerifier, 0);
addToRunTimeSelectionTable(ecgVerificationModel, manufacturedLeadFieldECGVerifier, dictionary);
manufacturedLeadFieldECGVerifier::manufacturedLeadFieldECGVerifier
(
    const electroStateProvider& state, const word& name, const dictionary& dict,
    const wordList& channels, const List<vector>& positions
)
:
    ecgVerificationModel(state, name, channels, positions), scales_(),
    maxError_(0), samples_(0), output_()
{
    read(dict);
    if (Pstream::master())
    {
        wordList columns(3*channels.size());
        forAll(channels,i)
        {
            columns[3*i] = channels[i]+"_numeric";
            columns[3*i+1] = channels[i]+"_exact";
            columns[3*i+2] = channels[i]+"_error";
        }
        output_ = ecgModelIO::openTimeSeries(mesh_.time().globalPath()/"postProcessing",
            "manufacturedLeadFieldECG_"+name+".dat", columns);
    }
}
ecgVerificationModel::Requirements manufacturedLeadFieldECGVerifier::requirements() const
{
    Requirements r; r.needVm = true; return r;
}
bool manufacturedLeadFieldECGVerifier::read(const dictionary& dict)
{
    const dictionary& config = dict.subDict("verificationModel");
    scales_ = scalarField(config.lookup("scales"));
    if (scales_.size() != electrodeNames_.size())
    {
        FatalIOErrorInFunction(config) << "scales must contain one value per lead-field channel."
            << exit(FatalIOError);
    }
    return true;
}
void manufacturedLeadFieldECGVerifier::record(const List<scalar>& values)
{
    record(mesh_.time().value(), values);
}
void manufacturedLeadFieldECGVerifier::record(scalar time, const List<scalar>& values)
{
    // Independent continuum integral for FDA's (1,2,3) cosine voltage on
    // the unit cube and L=kappa*(sin(pi*x)cos(2pi*y)cos(3pi*z),
    // 2cos(pi*x)sin(2pi*y)cos(3pi*z),3cos(pi*x)cos(2pi*y)sin(3pi*z)).
    // Use the existing source oracle for its exact temporal normalization.
    scalarField x(1,0), y(1,0), z(1,0), v;
    computeManufacturedV(v,x,y,z,time,3);
    scalarField row(3*values.size());
    forAll(values,i)
    {
        const scalar exact = -14*constant::mathematical::pi*scales_[i]*v[0]/8;
        const scalar err = mag(values[i]-exact);
        maxError_ = max(maxError_,err);
        row[3*i] = values[i]; row[3*i+1] = exact; row[3*i+2] = err;
    }
    ++samples_;
    if (Pstream::master()) ecgModelIO::writeRow(output_.ref(),time,row);
}
void manufacturedLeadFieldECGVerifier::end()
{
    if (Pstream::master())
    {
        OFstream out(mesh_.time().globalPath()/"postProcessing"/
            ("manufacturedLeadFieldECGSummary_"+domainName_+".dat"));
        out.precision(12);
        out << "# Manufactured lead-field ECG summary" << nl
            << "samples " << samples_ << nl << "maxError " << maxError_ << nl;
    }
}
}
