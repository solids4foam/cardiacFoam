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

#ifndef hemoMetricsFunctionObject_H
#define hemoMetricsFunctionObject_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "wallFvPatch.H"
#include "OFstream.H"

namespace Foam
{
namespace functionObjects
{

class hemoMetricsFunctionObject : public fvMeshFunctionObject
{
    hemoMetricsFunctionObject(const hemoMetricsFunctionObject&) = delete;
    void operator=(const hemoMetricsFunctionObject&) = delete;

protected:
    // User options
    wordList patchNames_;
    scalar rho_;                 // density [kg/m3]
    scalar mu_;                  // dynamic viscosity [Pa·s]
    word UName_;                 // velocity field name

    // Cycle handling
    scalar cyclePeriod_;         // [s]
    scalar startPhase_;          // [s]

    // Fields
    autoPtr<volVectorField> wss_;    // wall shear stress vector
    autoPtr<volScalarField> magWSS_; // |WSS|

    autoPtr<volScalarField> intMagWSS_; // ∫|WSS| dt
    autoPtr<volVectorField> intWSS_;    // ∫WSS dt

    autoPtr<volScalarField> TAWSS_;
    autoPtr<volScalarField> OSI_;
    autoPtr<volScalarField> RRT_;

    // Timing
    scalar cycleStartTime_{-GREAT};

    // Helpers
    void ensureFields();
    void computeWSS();               // uses constant mu_
    void accumulate();
    void finalizeCycleAndWrite();

    bool isWallPatch(const fvPatch& pp) const;

public:
    TypeName("hemoMetrics");

    hemoMetricsFunctionObject(const word& name, const Time&, const dictionary&);

    // Read dictionary changes
    virtual bool read(const dictionary&);

    // Hooks
    virtual void start();
    virtual bool execute();
    virtual bool end();
    virtual bool write();
};

} // namespace functionObjects
} // namespace Foam

#endif
