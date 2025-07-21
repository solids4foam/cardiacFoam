/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "BuenoOrovioIonicModel.H"
#include "heaviside.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(BuenoOrovioIonicModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::BuenoOrovioIonicModel::BuenoOrovioIonicModel
(
    const dictionary& dict,
    const label num
)
:
    ODESystem(),
    dict_(dict),
    nIntegrationPoints_(num),
    monitorID_(readInt(dict.lookup("monitorID"))),
    monitorFile_("monitorFile"),
    odeSolver_(ODESolver::New(*this, dict_)),
    gatingVariables_(num),
    gatingDerivatives_(num),
    u_(0.0),
    uO_(readScalar(dict.lookup("uO"))),
    uU_(readScalar(dict.lookup("uU"))),
    thetaV_(readScalar(dict.lookup("thetaV"))),
    thetaW_(readScalar(dict.lookup("thetaW"))),
    thetaVMinus_(readScalar(dict.lookup("thetaVMinus"))),
    thetaO_(readScalar(dict.lookup("thetaO"))),
    tauV1Minus_(readScalar(dict.lookup("tauV1Minus"))),
    tauV2Minus_(readScalar(dict.lookup("tauV2Minus"))),
    tauVPlus_(readScalar(dict.lookup("tauVPlus"))),
    tauW1Minus_(readScalar(dict.lookup("tauW1Minus"))),
    tauW2Minus_(readScalar(dict.lookup("tauW2Minus"))),
    kWMinus_(readScalar(dict.lookup("kWMinus"))),
    uWMinus_(readScalar(dict.lookup("uWMinus"))),
    tauWPlus_(readScalar(dict.lookup("tauWPlus"))),
    tauFi_(readScalar(dict.lookup("tauFi"))),
    tauO1_(readScalar(dict.lookup("tauO1"))),
    tauO2_(readScalar(dict.lookup("tauO2"))),
    tauSo1_(readScalar(dict.lookup("tauSo1"))),
    tauSo2_(readScalar(dict.lookup("tauSo2"))),
    kSo_(readScalar(dict.lookup("kSo"))),
    uSo_(readScalar(dict.lookup("uSo"))),
    tauS1_(readScalar(dict.lookup("tauS1"))),
    tauS2_(readScalar(dict.lookup("tauS2"))),
    kS_(readScalar(dict.lookup("kS"))),
    uS_(readScalar(dict.lookup("uS"))),
    tauSi_(readScalar(dict.lookup("tauSi"))),
    tauWInfty_(readScalar(dict.lookup("tauWInfty"))),
    wInftyStar_(readScalar(dict.lookup("wInftyStar")))
{
    // Lookup the initial values for the gating variables
    const scalar initialV(readScalar(dict.lookup("initialV")));
    const scalar initialW(readScalar(dict.lookup("initialW")));
    const scalar initialS(readScalar(dict.lookup("initialS")));

    // Initialise the gating variables and derivatives
    // Each integration point has three gating variables
    forAll(gatingVariables_, gatingVarI)
    {
        gatingVariables_.set(gatingVarI, new scalarField(nEqns(), 0.0));
        gatingDerivatives_.set(gatingVarI, new scalarField(nEqns(), 0.0));

        gatingVariables_[gatingVarI][0] = initialV;
        gatingVariables_[gatingVarI][1] = initialW;
        gatingVariables_[gatingVarI][2] = initialS;
    }

    // Check monitor is within range
    if (monitorID_ < 0 && monitorID_ > (nIntegrationPoints_ - 1))
    {
        FatalError
            << "monitorID is out of range!" << abort(FatalError);
    }

    // Write header to monitor file
    monitorFile_
        << "# t step u v w s Jfi Jso Jsi" << endl;
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::BuenoOrovioIonicModel::calculateCurrent
(
    const scalar stepStartTime,
    const scalar deltaT,
    const scalarField& Vm,
    scalarField& totalJ
)
{
    if (totalJ.size() != nIntegrationPoints_)
    {
        FatalErrorIn("void Foam::BuenoOrovioIonicModel::calculateCurrent(...)")
            << "totalJ.size() != nIntegrationPoints" << abort(FatalError);
    }

    if (Vm.size() != nIntegrationPoints_)
    {
        FatalErrorIn("void Foam::BuenoOrovioIonicModel::calculateCurrent(...)")
            << "Vm.size() != nIntegrationPoints" << abort(FatalError);
    }

    // Update the ODE system for each integration point
    // TODO: this only makes sense if the time-stpe is solved once: otherwise I
    // need to store the old values and only update them for new time-steps
    // One solution would be to use volScalarFields; an added benefit is
    // visualisation, but a downside is that we cannot easily use pointFields or
    // surfaceFields
    forAll(gatingVariables_, integrationPtI)
    {
        // Take a reference to the gating variables for this integration point
        scalarField& yStart = gatingVariables_[integrationPtI];
        // scalarField& dyStart = gatingDerivatives_[integrationPtI];

        // Set t to the initial time and step size
        // Note: ODE solves updates these so we re-create them for each point
        // scalar t = stepStartTime;
        const scalar tStart = stepStartTime;
        const scalar tEnd = stepStartTime + deltaT;

        // Set step to deltaT
        scalar step = deltaT;

        // Take references to improve readability
        scalar& u = u_;
        const scalar& v = yStart[0];
        const scalar& w = yStart[1];
        const scalar& s = yStart[2];

        // Store the voltage for this integration point as a way to pass it to
        // the derivatives function
        u = (1000*Vm[integrationPtI] + 84)/85.7;
        //u = Vm[integrationPtI];

        // Calculate ODE derivatives: do I need to do this? Is it for
        // initialising yStart?
        //derivatives(t, yStart, dyStart);

        // Update ODE system
        // Note: this updates the value of t to the end time
        // TODO: We should check if it suceeded
        odeSolver_->solve(tStart, tEnd, yStart, step);
        //odeSolver_->solve(t, yStart, step);

        // Time constants
        const scalar tauSo =
            tauSo1_
          + (tauSo2_ - tauSo1_)*(1.0 + tanh(kSo_*(u - uSo_)))/2.0;
        const scalar tauO =
            (1.0 - heaviside(u - thetaO_))*tauO1_ + heaviside(u - thetaO_)*tauO2_;

        // Calculate the three currents
        const scalar Jfi =
            -v*heaviside(u - thetaV_)*(u - thetaV_)*(uU_ - u)/tauFi_;
        const scalar Jso =
            (u - uO_)*(1.0 - heaviside(u - thetaW_))/tauO
          + heaviside(u - thetaW_)/tauSo;
        const scalar Jsi = -heaviside(u - thetaW_)*w*s/tauSi_;

        if (integrationPtI == monitorID_)
        {
            monitorFile_
                << tEnd << " "
                << step << " "
                << u << " "
                << yStart[0] << " "
                << yStart[1] << " "
                << yStart[2] << " "
                << Jfi << " "
                << Jso << " "
                << Jsi << endl;
        }

        // Update the total current
        totalJ[integrationPtI] = Jfi + Jso + Jsi;
    }
}


void Foam::BuenoOrovioIonicModel::derivatives
(
    const scalar t,
    const scalarField& y,
    scalarField& dydt
) const
{
    // Note: it is assumed that u_ holds the voltage for the current integration point

    // Take references to improve readability
    const scalar& u = u_;
    const scalar& v = y[0];
    const scalar& w = y[1];
    const scalar& s = y[2];
    scalar& dvdt = dydt[0];
    scalar& dwdt = dydt[1];
    scalar& dsdt = dydt[2];

    // Calculate the time constants
    const scalar tauVMinus =
        (1.0 - heaviside(u - thetaVMinus_))*tauV1Minus_
      + heaviside(u - thetaVMinus_)*tauV2Minus_;
    const scalar tauWMinus =
        tauW1Minus_
      + (tauW2Minus_ - tauW1Minus_)*(1.0 + tanh(kWMinus_*(u - uWMinus_)))/2.0;
    const scalar tauS =
        (1.0 - heaviside(u - thetaW_))*tauS1_ + heaviside(u - thetaW_)*tauS2_;
    const scalar vInfty = heaviside(thetaVMinus_ - u);
    const scalar wInfty =
        (1.0 - heaviside(u - thetaO_))*(1.0 - u/tauWInfty_)
      + heaviside(u - thetaO_)*wInftyStar_;

    // Calculate current derivatives
    dvdt =
        (1.0 - heaviside(u - thetaV_))*(vInfty - v)/tauVMinus
      - heaviside(u - thetaV_)*v/tauVPlus_;
    dwdt =
        (1.0 - heaviside(u - thetaW_))*(wInfty - w)/tauWMinus
      - heaviside(u - thetaW_)*w/tauWPlus_;
    dsdt = ((1 + tanh(kS_*(u - uS_)))/2.0 - s)/tauS;
}


void Foam::BuenoOrovioIonicModel::jacobian
(
    const scalar t,
    const scalarField& y,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    notImplemented("void Foam::BuenoOrovioIonicModel::jacobian(...)");

    // const scalar T = c[nSpecie_];
    // const scalar p = c[nSpecie_ + 1];

    // forAll(c_, i)
    // {
    //     c_[i] = max(c[i], 0.0);
    // }

    // dfdc = Zero;

    // // Length of the first argument must be nSpecie_
    // omega(c_, T, p, dcdt);

    // forAll(reactions_, ri)
    // {
    //     const Reaction<ThermoType>& R = reactions_[ri];

    //     const scalar kf0 = R.kf(p, T, c_);
    //     const scalar kr0 = R.kr(kf0, p, T, c_);

    //     forAll(R.lhs(), j)
    //     {
    //         const label sj = R.lhs()[j].index;
    //         scalar kf = kf0;
    //         forAll(R.lhs(), i)
    //         {
    //             const label si = R.lhs()[i].index;
    //             const scalar el = R.lhs()[i].exponent;
    //             if (i == j)
    //             {
    //                 if (el < 1.0)
    //                 {
    //                     if (c_[si] > SMALL)
    //                     {
    //                         kf *= el*pow(c_[si], el - 1.0);
    //                     }
    //                     else
    //                     {
    //                         kf = 0.0;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     kf *= el*pow(c_[si], el - 1.0);
    //                 }
    //             }
    //             else
    //             {
    //                 kf *= pow(c_[si], el);
    //             }
    //         }

    //         forAll(R.lhs(), i)
    //         {
    //             const label si = R.lhs()[i].index;
    //             const scalar sl = R.lhs()[i].stoichCoeff;
    //             dfdc(si, sj) -= sl*kf;
    //         }
    //         forAll(R.rhs(), i)
    //         {
    //             const label si = R.rhs()[i].index;
    //             const scalar sr = R.rhs()[i].stoichCoeff;
    //             dfdc(si, sj) += sr*kf;
    //         }
    //     }

    //     forAll(R.rhs(), j)
    //     {
    //         const label sj = R.rhs()[j].index;
    //         scalar kr = kr0;
    //         forAll(R.rhs(), i)
    //         {
    //             const label si = R.rhs()[i].index;
    //             const scalar er = R.rhs()[i].exponent;
    //             if (i == j)
    //             {
    //                 if (er < 1.0)
    //                 {
    //                     if (c_[si] > SMALL)
    //                     {
    //                         kr *= er*pow(c_[si], er - 1.0);
    //                     }
    //                     else
    //                     {
    //                         kr = 0.0;
    //                     }
    //                 }
    //                 else
    //                 {
    //                     kr *= er*pow(c_[si], er - 1.0);
    //                 }
    //             }
    //             else
    //             {
    //                 kr *= pow(c_[si], er);
    //             }
    //         }

    //         forAll(R.lhs(), i)
    //         {
    //             const label si = R.lhs()[i].index;
    //             const scalar sl = R.lhs()[i].stoichCoeff;
    //             dfdc(si, sj) += sl*kr;
    //         }
    //         forAll(R.rhs(), i)
    //         {
    //             const label si = R.rhs()[i].index;
    //             const scalar sr = R.rhs()[i].stoichCoeff;
    //             dfdc(si, sj) -= sr*kr;
    //         }
    //     }
    // }

    // // Calculate the dcdT elements numerically
    // const scalar delta = 1.0e-3;

    // omega(c_, T + delta, p, dcdt_);
    // for (label i=0; i<nSpecie_; i++)
    // {
    //     dfdc(i, nSpecie_) = dcdt_[i];
    // }

    // omega(c_, T - delta, p, dcdt_);
    // for (label i=0; i<nSpecie_; i++)
    // {
    //     dfdc(i, nSpecie_) = 0.5*(dfdc(i, nSpecie_) - dcdt_[i])/delta;
    // }

    // dfdc(nSpecie_, nSpecie_) = 0;
    // dfdc(nSpecie_ + 1, nSpecie_) = 0;
}


// ************************************************************************* //
