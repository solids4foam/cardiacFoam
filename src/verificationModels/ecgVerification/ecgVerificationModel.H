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

#ifndef ecgVerificationModel_H
#define ecgVerificationModel_H

#include "List.H"
#include "autoPtr.H"
#include "electroStateProvider.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class ecgVerificationModel
{
public:

    struct Requirements
    {
        bool needVm = false;
        bool needPhiE = false;
        bool needConductivity = false;
    };

protected:

    const electroStateProvider& stateProvider_;
    const fvMesh& mesh_;
    wordList electrodeNames_;
    List<vector> electrodePositions_;

    const volScalarField& requireVm() const;
    const volScalarField& requirePhiE() const;
    const volTensorField& requireConductivity() const;

public:

    TypeName("ecgVerificationModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ecgVerificationModel,
        dictionary,
        (
            const electroStateProvider& stateProvider,
            const dictionary& dict,
            const wordList& electrodeNames,
            const List<vector>& electrodePositions
        ),
        (stateProvider, dict, electrodeNames, electrodePositions)
    );

    ecgVerificationModel
    (
        const electroStateProvider& stateProvider,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual ~ecgVerificationModel() = default;

    static autoPtr<ecgVerificationModel> New
    (
        const electroStateProvider& stateProvider,
        const dictionary& dict,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    static word selectedType(const dictionary& dict);

    virtual Requirements requirements() const = 0;

    void validateProvider() const;

    void updateElectrodes
    (
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual bool read(const dictionary& dict) = 0;

    virtual void record(const List<scalar>& numericValues) = 0;

    virtual void end() {}
};

} // End namespace Foam

#endif

// ************************************************************************* //
