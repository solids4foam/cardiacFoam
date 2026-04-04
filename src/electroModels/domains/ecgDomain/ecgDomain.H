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
    Foam::ECGDomain

Description
    Abstract base class for ECG electro domains.

    Follows the same runtime-selection pattern as ionicModel /
    activeTensionModel. Concrete subclasses implement evolve() to evaluate
    body-surface potentials from the electrophysiology state exposed by an
    electroStateProvider.

    ECGDomain is a post-primary electrophysiology domain fed one-way by
    the electrophysiology state exposed through an electroStateProvider.
    Current implementations consume Vm and conductivity directly from the
    provider. Future torso-domain ECG models may additionally use a dedicated
    one-way coupling model.

    Common members available to all subclasses:
      mesh_             — const ref to the heart fvMesh
      Vm()              — transmembrane voltage from the state provider
      conductivity()    — conductivity tensor from the state provider
      electrodeNames_   — electrode labels  (empty if not configured)
      electrodePositions_ — electrode positions (empty if not configured)

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    ECGDomain.C
\*---------------------------------------------------------------------------*/

#ifndef ECGDomain_H
#define ECGDomain_H

#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "electroModel/electrophysicsDomain.H"
#include "electroModel/electroStateProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ECGDomain Declaration
\*---------------------------------------------------------------------------*/

class ECGDomain
:
    public electrophysicsDomain
{
protected:

    const electroStateProvider& stateProvider_;
    const fvMesh& mesh_;

    //- Electrode labels (empty when electrodes not configured)
    wordList electrodeNames_;

    //- Electrode positions (empty when electrodes not configured)
    List<vector> electrodePositions_;

    //- Read electrode positions from the 'electrodes' subdict (if present)
    void readElectrodes(const dictionary& dict);

public:

    TypeName("ECGDomain");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ECGDomain,
        dictionary,
        (
            const electroStateProvider& stateProvider,
            const dictionary& dict,
            const word& domainName
        ),
        (stateProvider, dict, domainName)
    );

    // Constructors

        //- Construct from an electrophysiology state provider and ECG
        //  domain coefficients dictionary.
        ECGDomain
        (
            const electroStateProvider& stateProvider,
            const dictionary& dict,
            const word& domainName
        );

    // Selectors

        //- Return a new ECGDomain selected from dict entry
        //  'ECGDomain'.
        static autoPtr<ECGDomain> New
        (
            const electroStateProvider& stateProvider,
            const dictionary& dict,
            const word& domainName = "ECG"
        );

    // Destructor
    virtual ~ECGDomain() = default;

    // Access

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        const volScalarField& Vm() const;

        const volTensorField& conductivity() const;

        const wordList& electrodeNames() const { return electrodeNames_; }
        const List<vector>& electrodePositions() const
        {
            return electrodePositions_;
        }

    // Member Functions

        virtual const Time& time() const override
        {
            return mesh_.time();
        }

        //- Advance the ECG domain for the current timestep.
        virtual void evolve
        (
            scalar t0,
            scalar dt
        ) = 0;

        virtual void advance
        (
            scalar t0,
            scalar dt
        ) override
        {
            evolve(t0, dt);
        }

        //- Refresh model coefficients from the ECG sub-dictionary.
        //  Subclasses may extend this, but should preserve output layout.
        virtual bool read(const dictionary& dict);

        //- Finalize any model-specific reporting at the end of the run
        virtual void end() override {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
