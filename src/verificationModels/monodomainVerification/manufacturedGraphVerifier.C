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

#include "manufacturedGraphVerifier.H"

#include "OFstream.H"
#include "OSspecific.H"
#include "ionicModel.H"
#include "monodomainVerification/manufacturedFDAReference.H"
#include "verificationUtils.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

using namespace verificationUtils;

defineTypeNameAndDebug(manufacturedGraphVerifier, 0);
addToRunTimeSelectionTable
(
    graphVerificationModel,
    manufacturedGraphVerifier,
    dictionary
);

manufacturedGraphVerifier::manufacturedGraphVerifier
(
    const dictionary& dict
)
:
    graphVerificationModel(dict),
    errorsReported_(false)
{}


void manufacturedGraphVerifier::preProcess
(
    const Time& runTime,
    ionicModel* model,
    scalarField& Vm1D,
    const pointField& nodeLocations
)
{
    if (!model)
    {
        return;
    }

    const label dimension = model->geometricDimension();
    const scalar t = runTime.value();

    scalarField X(nodeLocations.component(vector::X));
    scalarField Y(nodeLocations.component(vector::Y));
    scalarField Z(nodeLocations.component(vector::Z));

    computeManufacturedV(Vm1D, X, Y, Z, t, dimension);

    PtrList<scalarField> stateFields;
    wordList stateNames;
    PtrList<scalarField> algFields;
    wordList algNames;

    model->exportFields(stateNames, stateFields, algNames, algFields);

    label u1Field = -1, u2Field = -1, u3Field = -1;
    forAll(stateNames, i)
    {
        if (stateNames[i] == "u1") u1Field = i;
        else if (stateNames[i] == "u2") u2Field = i;
        else if (stateNames[i] == "u3") u3Field = i;
    }

    if (u1Field != -1 && u2Field != -1 && u3Field != -1)
    {
        computeManufacturedU(stateFields[u1Field], stateFields[u2Field], stateFields[u3Field], X, Y, Z, t, dimension);
        model->importFields(Vm1D, stateNames, stateFields, algNames, algFields);
    }
}


void manufacturedGraphVerifier::postProcess
(
    const Time& runTime,
    const ionicModel* model,
    const scalarField& Vm1D,
    const pointField& nodeLocations
)
{
    if (errorsReported_ || !model)
    {
        return;
    }

    const scalar t = runTime.value();
    const label dimension = model->geometricDimension();

    scalarField X(nodeLocations.component(vector::X));
    scalarField Y(nodeLocations.component(vector::Y));
    scalarField Z(nodeLocations.component(vector::Z));

    scalarField VmExact, u1Exact, u2Exact, u3Exact;
    computeManufacturedV(VmExact, X, Y, Z, t, dimension);
    computeManufacturedU(u1Exact, u2Exact, u3Exact, X, Y, Z, t, dimension);

    PtrList<scalarField> stateFields;
    wordList stateNames;
    PtrList<scalarField> algFields;
    wordList algNames;
    
    const_cast<ionicModel*>(model)->exportFields(stateNames, stateFields, algNames, algFields);

    label u1Field = -1, u2Field = -1, u3Field = -1;
    forAll(stateNames, i)
    {
        if (stateNames[i] == "u1") u1Field = i;
        else if (stateNames[i] == "u2") u2Field = i;
        else if (stateNames[i] == "u3") u3Field = i;
    }

    auto VmNorms = computeNorms(Vm1D, VmExact);
    Tuple2<Tuple2<scalar, scalar>, scalar> u1Norms(Tuple2<scalar, scalar>(0,0),0);
    Tuple2<Tuple2<scalar, scalar>, scalar> u2Norms(Tuple2<scalar, scalar>(0,0),0);

    if (u1Field != -1) u1Norms = computeNorms(stateFields[u1Field], u1Exact);
    if (u2Field != -1) u2Norms = computeNorms(stateFields[u2Field], u2Exact);

    const fileName outputDir(runTime.path()/"postProcessing");
    mkDir(outputDir);
    const fileName outputFile
    (
        outputDir
      / (
            "graph_"
          + dimensionName(dimension)
          + "_"
          + Foam::name(nodeLocations.size())
          + "_nodes.dat"
        )
    );

    if (Pstream::master())
    {
        Info << "\nManufactured-solution graph error summary (t = " << t << "):" << nl
             << "-------------------------------------------------" << nl
             << "Field     L1-error       L2-error       Linf-error" << nl
             << "Vm1D   " << VmNorms.first().first() << "   "
             << VmNorms.first().second() << "   " << VmNorms.second() << nl
             << "u1     " << u1Norms.first().first() << "   "
             << u1Norms.first().second() << "   " << u1Norms.second() << nl
             << "u2     " << u2Norms.first().first() << "   "
             << u2Norms.first().second() << "   " << u2Norms.second() << nl
             << "-------------------------------------------------" << endl;

        OFstream out(outputFile);
        out << "Manufactured-solution graph error summary (t = " << t << "):\n";
        out << "Field     L1-error       L2-error       Linf-error\n";
        out << "Vm1D   " << VmNorms.first().first() << "   "
            << VmNorms.first().second() << "   " << VmNorms.second() << "\n";
        out << "u1     " << u1Norms.first().first() << "   "
            << u1Norms.first().second() << "   " << u1Norms.second() << "\n";
        out << "u2     " << u2Norms.first().first() << "   "
            << u2Norms.first().second() << "   " << u2Norms.second() << "\n";
        out << "-------------------------------------------------\n\n";
    }

    errorsReported_ = true;
}

} // End namespace Foam

// ************************************************************************* //
