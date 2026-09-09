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
    manufacturedEikonalECGVerifier

Description
    Manufactured verifier for eikonal template-surrogate ECG.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedEikonalECGVerifier_H
#define manufacturedEikonalECGVerifier_H

#include "ecgVerificationModel.H"
#include "OFstream.H"
#include "Switch.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedEikonalECGVerifier
:
    public ecgVerificationModel
{
    autoPtr<OFstream> outputPtr_;

    Switch enabled_;
    label dimension_;
    label referenceQuadratureOrder_;
    List<label> checkQuadratureOrders_;
    vector k_;

    label sampleCount_;
    List<scalar> referenceErrorL1Sum_;
    List<scalar> referenceErrorL2Sum_;
    List<scalar> referenceErrorLinf_;
    List<List<scalar>> checkErrorL1Sum_;
    List<List<scalar>> checkErrorL2Sum_;
    List<List<scalar>> checkErrorLinf_;
    List<List<scalar>> referenceDeltaL1Sum_;
    List<List<scalar>> referenceDeltaL2Sum_;
    List<List<scalar>> referenceDeltaLinf_;

    bool summaryWritten_;

    void resizeStorage();
    void initialiseOutput();
    void updateStatistics
    (
        const List<scalar>& numericValues,
        const List<List<scalar>>& checkReferenceValues,
        const List<scalar>& referenceValues
    );
    void writeSummary();

    const dictionary& manufacturedDict(const dictionary& dict) const;

    manufacturedEikonalECGVerifier
    (
        const manufacturedEikonalECGVerifier&
    ) = delete;
    void operator=(const manufacturedEikonalECGVerifier&) = delete;

public:

    OverrideTypeName("manufacturedEikonalECGVerifier");

    manufacturedEikonalECGVerifier
    (
        const electroStateProvider& stateProvider,
        const word& domainName,
        const dictionary& dict,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual ~manufacturedEikonalECGVerifier() = default;

    virtual Requirements requirements() const override;

    virtual bool read(const dictionary& dict) override;

    // Base class override for implicitly timed workflows
    virtual void record(const List<scalar>& numericValues) override;

    // Overload for eikonal ECG, which is driven by a post-processing pass
    // over activation times rather than a time-stepping loop, and thus
    // requires explicit time stamps.

    virtual void record
    (
        scalar sampleTime,
        const List<scalar>& numericValues
    ) override;

    virtual void end() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
