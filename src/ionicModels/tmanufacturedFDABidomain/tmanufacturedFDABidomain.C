/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include <math.h>
#include "tmanufacturedFDABidomain.H"
#include "tmanufacturedFDABidomain_2014.H"
#include "tmanufacturedFDABidomain_2014Names.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "ionicSelector.H"
#include "volFields.H"

namespace Foam
{
    defineTypeNameAndDebug(tmanufacturedFDABidomain, 0);
    addToRunTimeSelectionTable(ionicModel, tmanufacturedFDABidomain, dictionary);
}

Foam::tmanufacturedFDABidomain::tmanufacturedFDABidomain
(
    const dictionary& dict,
    const label num,
    const scalar initialDeltaT,
    const Switch solveVmWithinODESolver
)
:
    ionicModel(dict, num, initialDeltaT, solveVmWithinODESolver),
    STATES_(num),
    CONSTANTS_(BIDOMAIN_NUM_CONSTANTS, 0.0),
    ALGEBRAIC_(num),
    RATES_(num),
    k_
    (
        dict.found("manufacturedBidomain")
      ? dict.subDict("manufacturedBidomain").lookupOrDefault<scalar>
        (
            "k",
            1.0/Foam::sqrt(2.0)
        )
      : dict.lookupOrDefault<scalar>("k", 1.0/Foam::sqrt(2.0))
    )
{
    setTissue(ionicSelector::selectDimension(dict, supportedDimensions()));

    forAll(STATES_, integrationPtI)
    {
        STATES_.set(integrationPtI, new scalarField(BIDOMAIN_NUM_STATES, 0.0));
        ALGEBRAIC_.set
        (
            integrationPtI,
            new scalarField(BIDOMAIN_NUM_ALGEBRAIC, 0.0)
        );
        RATES_.set
        (
            integrationPtI,
            new scalarField(BIDOMAIN_NUM_STATES, 0.0)
        );

        tmanufacturedFDABidomainInitConsts
        (
            CONSTANTS_.data(),
            RATES_[integrationPtI].data(),
            STATES_[integrationPtI].data(),
            tissue(),
            k_
        );
    }
}


Foam::tmanufacturedFDABidomain::~tmanufacturedFDABidomain()
{}


Foam::List<Foam::word> Foam::tmanufacturedFDABidomain::supportedDimensions() const
{
    return {"1D", "2D", "3D"};
}


void Foam::tmanufacturedFDABidomain::solveODE
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& Im
)
{
    const scalar tStart = stepStartTime;
    const scalar tEnd = tStart + deltaT;
    const label monitorCell = 0;

    forAll(STATES_, integrationPtI)
    {
        scalarField& S = STATES_[integrationPtI];
        scalarField& A = ALGEBRAIC_[integrationPtI];
        scalarField& R = RATES_[integrationPtI];

        scalar& h = ionicModel::step()[integrationPtI];

        S[BidomainV] = Vm[integrationPtI];
        h = min(h, deltaT);

        if (integrationPtI == monitorCell)
        {
            debugPrintFields(integrationPtI, tStart, tEnd, h);
        }

        odeSolver().solve(tStart, tEnd, S, h);

        ::tmanufacturedFDABidomainComputeVariables
        (
            tEnd,
            CONSTANTS_.data(),
            R.data(),
            S.data(),
            A.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        if (integrationPtI == monitorCell)
        {
            debugPrintFields(integrationPtI, tStart, tEnd, h);
        }

        Im[integrationPtI] = A[BidomainIion] / CONSTANTS_[BidomainCm];
    }
}


void Foam::tmanufacturedFDABidomain::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALG(BIDOMAIN_NUM_ALGEBRAIC, 0.0);

    ::tmanufacturedFDABidomainComputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        ALG.data(),
        tissue(),
        solveVmWithinODESolver()
    );
}

// ************************************************************************* //
