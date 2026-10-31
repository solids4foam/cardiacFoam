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
    Foam::templateSignalProvider

Description
    Vm and Cai of each cell from its region templates at t - activationTime.
    A cell with activationTime < 0 stays at the template's first sample.

SourceFiles
    templateSignalProvider.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef templateSignalProvider_H
#define templateSignalProvider_H

#include "electromechanicalSignalProvider.H"
#include "ionicTemplates.H"
#include "volFields.H"

namespace Foam
{

class templateSignalProvider
:
    public ElectromechanicalSignalProvider
{
    const volScalarField& activationTime_;
    List<ionicTemplates::DynamicTemplate> templates_;
    List<scalarField> weights_;

public:

    //- Construct from the mesh holding activationTime and a
    //  personalizedTemplates dictionary
    templateSignalProvider(const fvMesh& mesh, const dictionary& templatesDict);

    virtual ~templateSignalProvider() = default;

    virtual bool hasSignal(const CouplingSignal s) const override
    {
        return
            s == CouplingSignal::VM
         || (s == CouplingSignal::CAI && !templates_[0].valuesCai.empty());
    }

    virtual scalar signal(const label i, const CouplingSignal s) const override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
