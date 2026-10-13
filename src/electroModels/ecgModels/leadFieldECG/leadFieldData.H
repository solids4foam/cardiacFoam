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

#ifndef leadFieldData_H
#define leadFieldData_H

#include "electroStateProvider.H"
#include "dictionary.H"
#include "PtrList.H"
#include "volFields.H"

namespace Foam
{

// Fixed conductivity-weighted sensitivity, including reciprocal-current
// normalization but excluding cell volumes. No propagation state is modified.
class leadFieldData
{
    wordList names_;
    PtrList<volVectorField> fields_;
    autoPtr<volTensorField> Gi_;
    autoPtr<volTensorField> Ge_;
    vectorField signatures_;

public:
    leadFieldData(const electroStateProvider&, const dictionary&, bool heartOnly = false);
    const wordList& names() const { return names_; }
    void project(const electroStateProvider&, scalarField&) const;

    static tmp<volTensorField> readPhysical
    (
        const fvMesh& mesh,
        const dictionary& dict,
        bool intracellular
    );
    // Geometry is retained separately because heartOnly deliberately cuts bath faces.
    static vectorField meshSignatures(const fvMesh&);
    static vectorField topologySignatures
    (
        const fvMesh&, const scalarField& heartMask, bool fullConductor
    );
    static void validateFieldNames(const dictionary&);
    static void writeProvenance
    (
        const fvMesh&, const volTensorField&, const volTensorField&,
        const dictionary&,
        const word& gradientConstruction = "myocardialNeighbourLeastSquares"
    );
};

}
#endif
