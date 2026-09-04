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

Class
    Foam::LandNiedererTWorld

Description
    OpenFOAM wrapper for TWorld's six-state Land-style contraction subsystem
    implemented as a standalone activeTensionModel.

    The model is driven by [Ca2+]_i (in mM) from an upstream ionicModel
    through the ElectromechanicalSignalProvider (CouplingSignal::CAI).
    Sarcomere stretch lambda is received via the solveAtPoint interface and
    used to compute the Frank-Starling length-dependence factor (Lfac) and
    the length-sensitive Ca50 dynamically at each timestep.

    Six ODE state variables are maintained per integration point:
      Ca_TRPN   -- normalised Ca-troponin occupancy
      TmBlocked -- tropomyosin-blocked regulatory unit fraction
      XW        -- weak-binding crossbridge fraction
      XS        -- strong-binding crossbridge fraction
      ZETAS     -- velocity distortion, strong crossbridges
      ZETAW     -- velocity distortion, weak crossbridges

    Active tension (Ta, kPa) is an algebraic output computed from the
    crossbridge state at each step.

    Equations are extracted from the TWorld 2025 ionic model contraction
    subsystem. This is not the original seven-state Land et al. (2017)
    intact-human model; use LandNiederer for that model.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    LandNiedererTWorld.C
\*---------------------------------------------------------------------------*/

#ifndef LandNiedererTWorld_H
#define LandNiedererTWorld_H

#include "activeTensionModel.H"
#include "overrideTypeName.H"

#include "ODESystem.H"
#include "ODESolver.H"
#include "autoPtr.H"
#include "PtrList.H"
#include "scalarField.H"
#include "scalarMatrices.H"

#include "LandNiedererTWorld_2025Names.H"

namespace Foam
{

class LandNiedererTWorld
:
    public activeTensionModel,
    public ODESystem
{
    // Private data

        //- ODE solver
        autoPtr<ODESolver> odeSolver_;

        //- Per-integration-point ODE state storage
        PtrList<scalarField> STATES_;

        //- Per-integration-point algebraic variable storage
        PtrList<scalarField> ALGEBRAIC_;

        //- Per-integration-point ODE rate storage
        PtrList<scalarField> RATES_;

        //- Shared (broadcast) constant storage
        mutable scalarField CONSTANTS_;

        //- Per-integration-point previous lambda (for lambda_rate computation)
        scalarField prevLambda_;

        //- Lambda at current integration point (for derivatives() sub-steps)
        mutable scalar currentLambda_;

        //- lambda_rate at current integration point (for derivatives() sub-steps)
        mutable scalar currentLambdaRate_;


protected:

    // I/O metadata hooks (override base class nullptr defaults)

        virtual const PtrList<scalarField>* ioStatesPtr() const override
        {
            return &STATES_;
        }

        virtual const PtrList<scalarField>* ioAlgebraicPtr() const override
        {
            return &ALGEBRAIC_;
        }

        virtual const PtrList<scalarField>* ioRatesPtr() const override
        {
            return &RATES_;
        }

        virtual const scalarField* ioConstantsPtr() const override
        {
            return &CONSTANTS_;
        }

        virtual label ioNumConstants() const override
        {
            return NUM_CONSTANTS;
        }

        virtual const char* const* ioConstantNames() const override;

        virtual const char* const* ioStateNames() const override;

        virtual label ioNumStates() const override
        {
            return NUM_STATES;
        }

        virtual const char* const* ioAlgebraicNames() const override;

        virtual label ioNumAlgebraic() const override
        {
            return NUM_ALGEBRAIC;
        }

        virtual bool supportsRestartState() const override
        {
            return true;
        }

        virtual bool readRestartState(const fvMesh& mesh) override;

        virtual void writeRestartState(const fvMesh& mesh) const override;

        virtual bool restartTension(scalarField& Ta) const override;

        //- This model is driven by intracellular calcium.
        virtual CouplingSignal driveSignal() const override
        {
            return CouplingSignal::CAI;
        }

        //- Advance ODE and compute active tension for integration point i.
        //  driveVal : [Ca2+]_i in mM
        //  lambda   : sarcomere stretch ratio
        //  Ta       : output active tension [kPa]
        virtual void solveAtPoint
        (
            label i,
            scalar driveVal,
            scalar lambda,
            scalar& Ta
        ) override;


public:

    OverrideTypeName("LandNiedererTWorld");

    // Constructors

        LandNiedererTWorld
        (
            const dictionary& dict,
            const label nIntegrationPoints
        );

    // Destructor

        virtual ~LandNiedererTWorld() = default;


    // Member Functions

        //- Declare that this model needs [Ca2+]_i from the provider
        virtual Requirements requirements() const override;

        //- Run the ODE to resting steady state at the given Ca_i [mM].
        //  Called automatically by sequentialElectroMechanical after the
        //  ionic signal provider is wired, so that all integration points
        //  start from the correct resting equilibrium.
        virtual void preconditionToRestingState(const scalar restingCai) override;

        //- Number of ODEs (= NUM_STATES = 6)
        virtual label nEqns() const override;

        //- ODE RHS — called by the ODE solver during sub-stepping.
        //  Reads currentDriveSignal_ (Cai), currentLambda_, currentLambdaRate_.
        virtual void derivatives
        (
            const scalar t,
            const scalarField& y,
            scalarField& dydt
        ) const override;

        //- Jacobian not implemented; use explicit ODE solvers (e.g. RK45).
        virtual void jacobian
        (
            const scalar t,
            const scalarField& y,
            scalarField& dfdx,
            scalarSquareMatrix& dfdy
        ) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
