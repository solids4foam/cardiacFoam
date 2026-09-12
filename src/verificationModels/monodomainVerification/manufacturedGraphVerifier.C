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

namespace
{

label stateIndex
(
    const wordList& stateNames,
    const word& stateName
)
{
    const label idx = stateNames.find(stateName);

    if (idx < 0)
    {
        FatalErrorInFunction
            << "manufacturedGraphVerifier requires ionic state '"
            << stateName << "' but available states are "
            << stateNames
            << exit(FatalError);
    }

    return idx;
}


wordList stateNames(const ionicModel& model)
{
    if (!model.ioStateNames() || model.ioNumStates() <= 0)
    {
        FatalErrorInFunction
            << "manufacturedGraphVerifier requires ionic model "
            << model.type()
            << " to expose state metadata."
            << exit(FatalError);
    }

    wordList names(model.ioNumStates());
    const char* const* rawNames = model.ioStateNames();
    forAll(names, i)
    {
        names[i] = rawNames[i];
    }

    return names;
}

} // End anonymous namespace


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
    const pointField& nodeLocations,
    label localStartNode
)
{
    if (!model)
    {
        FatalErrorInFunction
            << "manufacturedGraphVerifier requires an ionic model."
            << exit(FatalError);
    }

    const label dimension = model->geometricDimension();
    const scalar t = runTime.value();

    scalarField X(nodeLocations.component(vector::X));
    scalarField Y(nodeLocations.component(vector::Y));
    scalarField Z(nodeLocations.component(vector::Z));

    computeManufacturedV(Vm1D, X, Y, Z, t, dimension);

    const PtrList<scalarField>* statesPtr = model->ioStatesPtr();
    if (!statesPtr)
    {
        FatalErrorInFunction
            << "manufacturedGraphVerifier requires ionic model "
            << model->type()
            << " to expose mutable state storage."
            << exit(FatalError);
    }

    const wordList names = stateNames(*model);
    const label u1Idx = stateIndex(names, "u1");
    const label u2Idx = stateIndex(names, "u2");
    const label u3Idx = stateIndex(names, "u3");

    PtrList<scalarField>& states =
        const_cast<PtrList<scalarField>&>(*statesPtr);

    scalarField localX(states.size(), 0.0);
    scalarField localY(states.size(), 0.0);
    scalarField localZ(states.size(), 0.0);

    forAll(states, localNodeI)
    {
        const label nodeI = localStartNode + localNodeI;

        if (nodeI < 0 || nodeI >= nodeLocations.size())
        {
            FatalErrorInFunction
                << "Local graph node " << nodeI
                << " is outside graph node range [0, "
                << nodeLocations.size() - 1 << "]."
                << exit(FatalError);
        }

        localX[localNodeI] = X[nodeI];
        localY[localNodeI] = Y[nodeI];
        localZ[localNodeI] = Z[nodeI];
    }

    scalarField u1Exact, u2Exact, u3Exact;
    computeManufacturedU
    (
        u1Exact,
        u2Exact,
        u3Exact,
        localX,
        localY,
        localZ,
        t,
        dimension
    );

    forAll(states, localNodeI)
    {
        states[localNodeI][u1Idx] = u1Exact[localNodeI];
        states[localNodeI][u2Idx] = u2Exact[localNodeI];
        states[localNodeI][u3Idx] = u3Exact[localNodeI];
    }
}


void manufacturedGraphVerifier::postProcess
(
    const Time& runTime,
    const ionicModel* model,
    const scalarField& Vm1D,
    const pointField& nodeLocations,
    label localStartNode
)
{
    if (errorsReported_)
    {
        return;
    }

    if (!model)
    {
        FatalErrorInFunction
            << "manufacturedGraphVerifier requires an ionic model."
            << exit(FatalError);
    }

    const scalar t = runTime.value();
    const label dimension = model->geometricDimension();

    scalarField X(nodeLocations.component(vector::X));
    scalarField Y(nodeLocations.component(vector::Y));
    scalarField Z(nodeLocations.component(vector::Z));

    scalarField VmExact, u1Exact, u2Exact, u3Exact;
    computeManufacturedV(VmExact, X, Y, Z, t, dimension);
    computeManufacturedU(u1Exact, u2Exact, u3Exact, X, Y, Z, t, dimension);

    const PtrList<scalarField>* statesPtr = model->ioStatesPtr();
    if (!statesPtr)
    {
        FatalErrorInFunction
            << "manufacturedGraphVerifier requires ionic model "
            << model->type()
            << " to expose state storage."
            << exit(FatalError);
    }

    const wordList names = stateNames(*model);
    const label u1Idx = stateIndex(names, "u1");
    const label u2Idx = stateIndex(names, "u2");

    scalarField u1(statesPtr->size(), 0.0);
    scalarField u2(statesPtr->size(), 0.0);
    scalarField u1LocalExact(statesPtr->size(), 0.0);
    scalarField u2LocalExact(statesPtr->size(), 0.0);

    forAll(*statesPtr, localNodeI)
    {
        const label nodeI = localStartNode + localNodeI;

        if (nodeI < 0 || nodeI >= nodeLocations.size())
        {
            FatalErrorInFunction
                << "Local graph node " << nodeI
                << " is outside graph node range [0, "
                << nodeLocations.size() - 1 << "]."
                << exit(FatalError);
        }

        u1[localNodeI] = (*statesPtr)[localNodeI][u1Idx];
        u2[localNodeI] = (*statesPtr)[localNodeI][u2Idx];
        u1LocalExact[localNodeI] = u1Exact[nodeI];
        u2LocalExact[localNodeI] = u2Exact[nodeI];
    }

    auto VmNorms = computeNorms(Vm1D, VmExact);
    auto u1Norms = computeNorms(u1, u1LocalExact);
    auto u2Norms = computeNorms(u2, u2LocalExact);

    const fileName outputDir(runTime.globalPath()/"postProcessing");
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
