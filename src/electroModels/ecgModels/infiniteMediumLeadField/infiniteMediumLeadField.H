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
    infiniteMediumLeadField

Description
    Lead field of an infinite homogeneous extracellular medium, shared by
    pseudoECG and eikonalECG.

    The unipolar potential at electrode e is

        phi_e = 1/(4 pi sigma_e) sum_cells (grad Vm & z_e) ,
        z_e   = (G & r) V/|r|^3 ,  r = x - x_e ,

    with G the domain conductivity, x and V the cell centre and volume, and
    sigma_e the scalar extracellular conductivity. A cell centre within VSMALL
    of an electrode contributes nothing.

    Reads from the ECG domain entry:
        sigmaExtracellular   required, positive

SourceFiles
    infiniteMediumLeadField.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef infiniteMediumLeadField_H
#define infiniteMediumLeadField_H

#include "dictionary.H"
#include "vectorField.H"
#include "scalarField.H"

namespace Foam
{

class ecgDomain;

class infiniteMediumLeadField
{
    scalar sigmaExtracellular_;
    bool leadVectorsCalculated_;
    List<List<vector>> leadVectors_;

    void calculateLeadVectors(const ecgDomain& domain);

public:

    //- Read sigmaExtracellular from the ECG domain entry dict
    explicit infiniteMediumLeadField(const dictionary& dict);

    //- 1/(4 pi sigma_e)
    scalar scale() const;

    //- Potential at each electrode of domain from the cell gradient of Vm,
    //  summed over all processors. The lead vectors are built on the first
    //  call.
    void evaluate
    (
        const ecgDomain& domain,
        const vectorField& gradVm,
        scalarField& values
    );
};

} // End namespace Foam

#endif

// ************************************************************************* //
