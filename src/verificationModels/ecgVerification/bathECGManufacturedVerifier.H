/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    bathECGManufacturedVerifier

Description
    Bath ECG manufactured verifier.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef bathECGManufacturedVerifier_H
#define bathECGManufacturedVerifier_H

#include "ecgVerificationModel.H"
#include "OFstream.H"
#include "Switch.H"
#include "overrideTypeName.H"

namespace Foam
{

class bathECGManufacturedVerifier
:
    public ecgVerificationModel
{
    autoPtr<OFstream> outputPtr_;

    Switch enabled_;
    scalar k_;
    scalar alpha_;
    Switch groundElectrode_;

    label sampleCount_;
    scalar fieldErrorL1Sum_;
    scalar fieldErrorL2Sum_;
    scalar fieldErrorLinf_;
    List<scalar> electrodeErrorL1Sum_;
    List<scalar> electrodeErrorL2Sum_;
    List<scalar> electrodeErrorLinf_;
    bool summaryWritten_;

    void initialiseOutput();
    void resizeStatistics();
    void writeSummary();

    bathECGManufacturedVerifier(const bathECGManufacturedVerifier&) = delete;
    void operator=(const bathECGManufacturedVerifier&) = delete;

public:

    OverrideTypeName("bathECGManufacturedVerifier");

    bathECGManufacturedVerifier
    (
        const electroStateProvider& stateProvider,
        const dictionary& dict,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual ~bathECGManufacturedVerifier() = default;

    virtual Requirements requirements() const override;

    virtual bool read(const dictionary& dict) override;

    virtual void record(const List<scalar>& numericValues) override;

    virtual void end() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
