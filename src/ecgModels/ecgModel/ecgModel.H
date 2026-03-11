/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ecgModel

Description
    Abstract base class for ECG computation models.

    Follows the same runtime-selection pattern as ionicModel /
    activeTensionModel.  Concrete subclasses implement compute() to evaluate
    body-surface potentials from the transmembrane voltage field Vm.

    ecgModel is a standalone library (src/ecgModels) with no dependency on
    electroModels.  It receives Vm as a const reference from the owning
    ecgElectro electroModel.

    Common members available to all subclasses:
      Vm_               — const ref to monodomain Vm field
      mesh_             — const ref to the heart fvMesh
      Gi()              — conductivity tensor used by the ECG model
                           (prefers monodomain conductivity field when provided,
                           else falls back to local Gi initialisation)
      sigmaT_           — torso conductivity
      electrodeNames_   — electrode labels  (empty if not configured)
      electrodePositions_ — electrode positions (empty if not configured)

    If monodomain conductivity is provided by the owning electroModel, it is
    reused directly. Otherwise Gi is initialised by trying READ_IF_PRESENT from
    disk, then falling back to a uniform tensor built from the 'Gi' symmTensor
    entry in the dict.

Author
    Simão Nieto de Castro, UCD.
    Philip Cardiff, UCD.

SourceFiles
    ecgModel.C

\*---------------------------------------------------------------------------*/

#ifndef ecgModel_H
#define ecgModel_H

#include "volFields.H"
#include "dimensionedScalar.H"
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
    // Private Member Functions

        //- Initialise Gi: READ_IF_PRESENT from disk, else uniform dict tensor
        static tmp<volTensorField> initialiseGi
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


protected:

    //- Transmembrane voltage (owned by ecgElectro via monoDomainElectro)
    const volScalarField& Vm_;

    //- Heart mesh
    const fvMesh& mesh_;

    //- Owned fallback conductivity tensor when monodomain conductivity is not
    //  provided by the electroModel
    autoPtr<volTensorField> GiOwnedPtr_;

    //- Conductivity tensor used by the ECG model (either monodomain
    //  conductivity or owned fallback Gi)
    const volTensorField* GiPtr_;

    //- Torso / infinite-medium conductivity [S/m]
    dimensionedScalar sigmaT_;

    //- Electrode labels (empty when electrodes not configured)
    wordList electrodeNames_;

    //- Electrode positions (empty when electrodes not configured)
    List<vector> electrodePositions_;

    //- Read electrode positions from the 'electrodes' subdict (if present)
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
            const volTensorField* conductivityPtr
        ),
        (Vm, dict, conductivityPtr)
    );

    // Constructors

        //- Construct from Vm field and model coefficients dictionary.
        //  If conductivityPtr is non-null, it is reused as Gi.
        ecgModel
        (
            const volScalarField& Vm,
            const dictionary& dict,
            const volTensorField* conductivityPtr
        );

    // Selectors

        //- Return a new ecgModel selected from dict entry 'ecgModel'
        static autoPtr<ecgModel> New
        (
            const volScalarField& Vm,
            const dictionary& dict,
            const volTensorField* conductivityPtr = nullptr
        );

    // Destructor
    virtual ~ecgModel() = default;

    // Access

        const volTensorField& Gi() const { return *GiPtr_; }
        const dimensionedScalar& sigmaT() const { return sigmaT_; }

    // Member Functions

        //- Compute ECG output for the current timestep
        virtual void compute() = 0;

        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
