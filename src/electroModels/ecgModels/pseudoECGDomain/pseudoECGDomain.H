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
    Foam::PseudoECGDomain

Description
    ECG model using the Gima-Rudy dipole approximation (pseudo-ECG).

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    PseudoECGDomain.C
\*---------------------------------------------------------------------------*/

#ifndef PseudoECGDomain_H
#define PseudoECGDomain_H

#include "ecgDomain.H"
#include "ecgSolver.H"
#include "overrideTypeName.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ecgVerificationModel;

/*---------------------------------------------------------------------------*\
                     Class PseudoECGDomain Declaration
\*---------------------------------------------------------------------------*/

class PseudoECGDomain
:
    public ECGDomain
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
    PseudoECGDomain(const PseudoECGDomain&);

    //- Disallow assignment
    void operator=(const PseudoECGDomain&);


public:

    //- Runtime type information
    OverrideTypeName("pseudoECGDomain");

    // Constructors

        //- Construct from components
        PseudoECGDomain
        (
            const electroStateProvider& stateProvider,
            const dictionary& dict,
            const word& domainName
        );

    //- Destructor
    virtual ~PseudoECGDomain();

    // Member Functions

        //- Advance the pseudo-ECG domain for the current time-step
        virtual void evolve
        (
            scalar t0,
            scalar dt
        ) override;

        //- Refresh coefficients while preserving the existing output layout
        virtual bool read(const dictionary& dict) override;

        //- Flush manufactured diagnostics at the end of the run
        virtual void end() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
