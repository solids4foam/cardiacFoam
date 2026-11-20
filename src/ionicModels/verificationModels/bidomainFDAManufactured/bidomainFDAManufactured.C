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

#include "bidomainFDAManufactured.H"
#include "bidomainFDAManufactured_2014.H"
#include "bidomainFDAManufactured_2014Names.H"
#include "addToRunTimeSelectionTable.H"
#include "ionicModel.H"
#include "ionicModelIO.H"
#include "ionicSelector.H"
#include "volFields.H"

#include <math.h>

namespace Foam
{
    defineTypeNameAndDebug(bidomainFDAManufactured, 0);
    addToRunTimeSelectionTable(ionicModel, bidomainFDAManufactured, dictionary);
}


Foam::bidomainFDAConstants Foam::bidomainFDABenchmarkConstants
(
    const dictionary& coeffs,
    const label dimension
)
{
    const word source
    (
        coeffs.getOrDefault<word>("conductivitySource", "uniform")
    );

    if (source != "uniform")
    {
        FatalIOErrorInFunction(coeffs)
            << "The FDA bidomain benchmark needs uniform conductivities, but "
            << "conductivitySource is '" << source << "'."
            << exit(FatalIOError);
    }

    const symmTensor Gi
    (
        dimensionedSymmTensor("conductivityIntracellular", coeffs).value()
    );
    const symmTensor Ge
    (
        dimensionedSymmTensor("conductivityExtracellular", coeffs).value()
    );

    const tensor GiT(Gi & tensor(I));
    const tensor GeT(Ge & tensor(I));

    // The separable manufactured fields admit no mixed conductivity between active directions
    for (label i = 0; i < dimension; ++i)
    {
        for (label j = 0; j < dimension; ++j)
        {
            if
            (
                i != j
             && (
                    mag(GiT(i, j)) > SMALL*mag(Gi)
                 || mag(GeT(i, j)) > SMALL*mag(Ge)
                )
            )
            {
                FatalIOErrorInFunction(coeffs)
                    << "The FDA bidomain benchmark needs conductivity tensors "
                    << "without mixed entries between active directions; got "
                    << "Gi = " << Gi << " and Ge = " << Ge << "."
                    << exit(FatalIOError);
            }
        }
    }

    // phiE = -k Vm + s(t) needs k (Gi + Ge) = Gi in every active direction, to input rounding
    const scalar k = GiT.xx()/(GiT.xx() + GeT.xx());
    for (label d = 0; d < dimension; ++d)
    {
        const scalar kd = GiT(d, d)/(GiT(d, d) + GeT(d, d));
        if (mag(kd - k) > 1e-6*k)
        {
            FatalIOErrorInFunction(coeffs)
                << "The FDA bidomain benchmark needs Ge proportional to Gi: "
                << "Gi/(Gi + Ge) is " << k << " in x but " << kd
                << " in active direction " << d << "."
                << exit(FatalIOError);
        }
    }

    // div(Gm grad Vm) = -sum_d Gm_dd (d pi)^2 Vm for Vm ~ cos(pi x) cos(2 pi y) cos(3 pi z)
    scalar beta = 0;
    for (label d = 0; d < dimension; ++d)
    {
        const scalar Gm = GiT(d, d)*GeT(d, d)/(GiT(d, d) + GeT(d, d));
        beta -= Gm*sqr((d + 1)*constant::mathematical::pi);
    }

    // The solver checks the units of chi and cm
    return bidomainFDAConstants
    {
        dimensionedScalar("chi", coeffs).value(),
        dimensionedScalar("cm", coeffs).value(),
        k,
        beta
    };
}


Foam::bidomainFDAManufactured::bidomainFDAManufactured
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
    RATES_(num)
{
    setTissue(ionicSelector::selectDimension(dict, supportedDimensions()));

    const bidomainFDAConstants benchmark
    (
        bidomainFDABenchmarkConstants(dict, tissue())
    );

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

        bidomainFDAManufacturedInitConsts
        (
            CONSTANTS_.data(),
            RATES_[integrationPtI].data(),
            STATES_[integrationPtI].data(),
            benchmark.Cm,
            benchmark.chi,
            benchmark.beta
        );
    }

    applyIonicConstantOverrides();
}


Foam::bidomainFDAManufactured::~bidomainFDAManufactured()
{}


Foam::List<Foam::word> Foam::bidomainFDAManufactured::supportedDimensions() const
{
    return {"1D", "2D", "3D"};
}


void Foam::bidomainFDAManufactured::solveODE
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

        setActiveVmRate(integrationPtI);

        if (integrationPtI == monitorCell)
        {
            debugPrintFields(integrationPtI, tStart, tEnd, h);
        }

        odeSolver().solve(tStart, tEnd, S, h);

        ::bidomainFDAManufacturedComputeVariables
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

    clearVmRate();
}


void Foam::bidomainFDAManufactured::evaluateIonicCurrent
(
    const scalar t,
    const scalarField& Vm,
    scalarField& Im
)
{
    scalarField S(BIDOMAIN_NUM_STATES, 0.0);
    scalarField A(BIDOMAIN_NUM_ALGEBRAIC, 0.0);
    scalarField R(BIDOMAIN_NUM_STATES, 0.0);

    forAll(STATES_, integrationPtI)
    {
        S = STATES_[integrationPtI];
        S[BidomainV] = Vm[integrationPtI];
        A = 0.0;
        R = 0.0;

        ::bidomainFDAManufacturedComputeVariables
        (
            t,
            CONSTANTS_.data(),
            R.data(),
            S.data(),
            A.data(),
            tissue(),
            solveVmWithinODESolver()
        );

        Im[integrationPtI] = A[BidomainIion] / CONSTANTS_[BidomainCm];
    }
}


void Foam::bidomainFDAManufactured::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    scalarField ALG(BIDOMAIN_NUM_ALGEBRAIC, 0.0);

    ::bidomainFDAManufacturedComputeVariables
    (
        t,
        CONSTANTS_.data(),
        dydt.data(),
        const_cast<scalarField&>(y).data(),
        ALG.data(),
        tissue(),
        solveVmWithinODESolver()
    );

    if (!solveVmWithinODESolver())
    {
        dydt[BidomainV] = activeVmRate();
    }
}

// ************************************************************************* //
