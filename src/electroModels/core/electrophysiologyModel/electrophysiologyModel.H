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
    Foam::electrophysiologyModel

Description
    Concrete electroModel for myocardium-centred electrophysiology runs.

    This class bridges the electroModel selection machinery with the domain
    system builder.  It is registered in the electroModel run-time selection
    table under \e "monodomainSolver", \e "bidomainSolver", and
    \e "eikonalSolver", keeping the base electroModel fully agnostic to
    domain-level concerns.

    Responsibilities:
      - own the IonicModel and optional verification model when required
      - own the output / post-process field lists
      - delegate spatial assembly to electrophysicsSystemBuilder
      - forward evolve() / setDeltaT() to the domain system

    The concrete myocardium domain is selected from the active Coeffs
    sub-dictionary. Reaction-diffusion runs construct MyocardiumDomain,
    while activation-time runs construct EikonalMyocardiumDomain.

Author
    Simao Nieto de Castro, UCD.

SourceFiles
    electrophysiologyModel.C
\*---------------------------------------------------------------------------*/

#ifndef electrophysiologyModel_H
#define electrophysiologyModel_H

#include "electroModel.H"
#include "ionicModel.H"
#include "electroVerificationModel.H"
#include "overrideTypeName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class electrophysiologyModel Declaration
\*---------------------------------------------------------------------------*/

class electrophysiologyModel
:
    public electroModel
{
    //- Run-time selectable ionic model
    autoPtr<ionicModel> ionicModelPtr_;

    //- Optional verification / post-processing model
    autoPtr<electroVerificationModel> verificationModelPtr_;

    //- Ionic state export fields (written to disk)
    PtrList<volScalarField> outFields_;

    //- Scratch fields for verification post-processing
    PtrList<volScalarField> postProcessFields_;

    //- Disallow copy
    electrophysiologyModel(const electrophysiologyModel&) = delete;
    void operator=(const electrophysiologyModel&) = delete;

    //- Pre-read the myocardiumSolver keyword so the correct Coeffs sub-dict
    //  is passed to the base class constructor.
    static word readSolverType(Time& runTime, const word& region);

public:

    //- Runtime type information
    //  The class registers under myocardiumSolver dict keys via
    //  addNamedToRunTimeSelectionTable, so the typeName here is for debug
    //  messages only.
    OverrideTypeName("electrophysiologyModel");

    // Constructors

        //- Construct from time and region.
        //  Reads myocardiumSolver from electroProperties to determine which
        //  Coeffs sub-dict to use, then builds the full domain system.
        electrophysiologyModel
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

        virtual void setDeltaT(Time& runTime) override;

    //- Destructor
    virtual ~electrophysiologyModel() = default;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
