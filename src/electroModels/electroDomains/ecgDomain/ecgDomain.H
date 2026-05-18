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
    Foam::ecgDomain

Description
    Concrete ECG domain class.

    Owns the electrode configuration, output I/O, and an ecgSolver kernel
    selected at runtime from the 'ecgSolver' dictionary key.  This mirrors
    the pattern used by myocardiumDomain: the domain is a single concrete
    class that holds state and delegates the numerical algorithm to a
    runtime-selectable ecgSolver.

    To add a new ECG algorithm (e.g. torsoECGSolver), register a new
    ecgSolver subclass — no changes to ecgDomain are needed.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    ecgDomain.C
\*---------------------------------------------------------------------------*/

#ifndef ecgDomain_H
#define ecgDomain_H

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
                       Class ecgDomain Declaration
\*---------------------------------------------------------------------------*/

class ecgDomain
:
    public electroDomainInterface,
    public electroStateProvider
{
    const electroStateProvider& stateProvider_;
    const fvMesh& mesh_;

    word solverType_;

    //- Output stream for pseudoECG.dat or bathECG.dat
    autoPtr<OFstream> outputPtr_;
    //- Runtime-selected ECG evaluation kernel
    autoPtr<ecgSolver> solverPtr_;
    //- Optional ECG verification extension
    autoPtr<ecgVerificationModel> verificationModelPtr_;
    //- Reused output values buffer
    scalarField numericValues_;

    const dictionary* inheritedElectrodePositionsPtr_;
    const dictionary* inheritedManufacturedBidomainPtr_;
    const dictionary* inheritedPotentialDomainPtr_;

    //- Disallow copy construct
    ecgDomain(const ecgDomain&) = delete;

    //- Disallow assignment
    void operator=(const ecgDomain&) = delete;

    //- Electrode labels (empty when electrodes not configured)
    wordList electrodeNames_;

    //- Electrode positions (empty when electrodes not configured)
    List<vector> electrodePositions_;

    //- Read electrode positions from the 'electrodes' subdict (if present)
    void readElectrodes(const dictionary& dict);

public:

    TypeName("ecgDomain");

    // Constructors

        ecgDomain
        (
            const electroStateProvider& stateProvider,
            const dictionary& dict,
            const word& domainName = "ECG",
            const dictionary* inheritedElectrodePositionsPtr = nullptr,
            const dictionary* inheritedManufacturedBidomainPtr = nullptr,
            const dictionary* inheritedPotentialDomainPtr = nullptr
        );

    // Destructor
    virtual ~ecgDomain();

    // Access

        virtual const fvMesh& mesh() const override
        {
            return mesh_;
        }

        virtual const fvMesh& baseMesh() const override
        {
            return stateProvider_.baseMesh();
        }

        virtual const labelUList* subsetFaceMapPtr() const override
        {
            return stateProvider_.subsetFaceMapPtr();
        }

        virtual const labelUList* subsetCellMapPtr() const override
        {
            return stateProvider_.subsetCellMapPtr();
        }

        const volScalarField& Vm() const;

        const volTensorField& conductivity() const;

        const volScalarField& phiE() const;

        const wordList& electrodeNames() const { return electrodeNames_; }
        const List<vector>& electrodePositions() const
        {
            return electrodePositions_;
        }

        virtual const volScalarField* VmPtr() const override
        {
            return stateProvider_.VmPtr();
        }

        virtual const volScalarField* phiEPtr() const override
        {
            return stateProvider_.phiEPtr();
        }

        virtual const volTensorField* conductivityPtr() const override
        {
            return stateProvider_.conductivityPtr();
        }

        virtual const volTensorField* intracellularConductivityPtr()
            const override
        {
            return stateProvider_.intracellularConductivityPtr();
        }

        virtual const volTensorField* extracellularConductivityPtr()
            const override
        {
            return stateProvider_.extracellularConductivityPtr();
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

        virtual void write() override {}

        //- Finalize any model-specific reporting at the end of the run
        virtual void end() override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
