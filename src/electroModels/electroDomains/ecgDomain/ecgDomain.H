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
    Concrete ECG domain class.

    Owns the electrode configuration, output I/O, and an ECGSolver kernel
    selected at runtime from the 'ecgSolver' dictionary key.  This mirrors
    the pattern used by MyocardiumDomain: the domain is a single concrete
    class that holds state and delegates the numerical algorithm to a
    runtime-selectable ECGSolver.

    To add a new ECG algorithm (e.g. TorsoECGSolver), register a new
    ECGSolver subclass — no changes to ECGDomain are needed.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    ecgDomain.C
\*---------------------------------------------------------------------------*/

#ifndef ECGDomain_H
#define ECGDomain_H

#include "volFields.H"
#include "autoPtr.H"
#include "OFstream.H"
#include "ecgSolver.H"
#include "electroDomainInterface.H"
#include "electroStateProvider.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ecgVerificationModel;

/*---------------------------------------------------------------------------*\
                       Class ECGDomain Declaration
\*---------------------------------------------------------------------------*/

class ECGDomain
:
    public electroDomainInterface
{
    //- Output stream for pseudoECG.dat
    autoPtr<OFstream> outputPtr_;

    //- Runtime-selected ECG evaluation kernel
    autoPtr<ECGSolver> solverPtr_;

    //- Optional ECG verification extension
    autoPtr<ecgVerificationModel> verificationModelPtr_;

    //- Reused output values buffer
    scalarField numericValues_;

    //- Disallow copy construct
    ECGDomain(const ECGDomain&) = delete;

    //- Disallow assignment
    void operator=(const ECGDomain&) = delete;

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

    // Constructors

        ECGDomain
        (
            const electroStateProvider& stateProvider,
            const dictionary& dict,
            const word& domainName = "ECG"
        );

    // Destructor
    virtual ~ECGDomain();

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
        );

        virtual void advance
        (
            scalar t0,
            scalar dt
        ) override
        {
            evolve(t0, dt);
        }

        //- Refresh model coefficients from the ECG sub-dictionary.
        virtual bool read(const dictionary& dict);

        //- Finalize any model-specific reporting at the end of the run
        virtual void end() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
