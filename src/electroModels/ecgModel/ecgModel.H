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
    Foam::ecgModel

Description
    Abstract base class for ECG computation models.

    Follows the same runtime-selection pattern as ionicModel /
    activeTensionModel. Concrete subclasses implement compute() to evaluate
    body-surface potentials from the transmembrane voltage field Vm.

Author
    Simao Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    ecgModel.C
\*---------------------------------------------------------------------------*/

#ifndef ecgModel_H
#define ecgModel_H

#include "volFields.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ecgModel Declaration
\*---------------------------------------------------------------------------*/

class ecgModel
{
protected:

    //- Transmembrane voltage owned by the parent electroModel
    const volScalarField& Vm_;

    //- Heart mesh
    const fvMesh& mesh_;

    //- Monodomain conductivity tensor provided by the parent electroModel
    const volTensorField& conductivity_;

    //- Electrode labels (empty when electrodes are not configured)
    wordList electrodeNames_;

    //- Electrode positions (empty when electrodes are not configured)
    List<vector> electrodePositions_;

    //- Read electrode positions from the electrodes sub-dictionary
    void readElectrodes(const dictionary& dict);


public:

    TypeName("ecgModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ecgModel,
        dictionary,
        (
            const volScalarField& Vm,
            const dictionary& dict,
            const volTensorField& conductivity
        ),
        (Vm, dict, conductivity)
    );

    // Constructors

        //- Construct from Vm field, coefficients dictionary, and conductivity
        ecgModel
        (
            const volScalarField& Vm,
            const dictionary& dict,
            const volTensorField& conductivity
        );

    // Selectors

        //- Return a new ecgModel selected from dict entry ecgModel
        static autoPtr<ecgModel> New
        (
            const volScalarField& Vm,
            const dictionary& dict,
            const volTensorField& conductivity
        );

    //- Destructor
    virtual ~ecgModel() = default;

    // Member Functions

        //- Return the conductivity tensor used by the ECG model
        const volTensorField& conductivity() const
        {
            return conductivity_;
        }

        //- Compute ECG output for the current time-step
        virtual void compute() = 0;

        //- Refresh model coefficients from the ECG sub-dictionary
        virtual bool read(const dictionary& dict);

        //- Finalise any model-specific reporting at the end of the run
        virtual void end() {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
