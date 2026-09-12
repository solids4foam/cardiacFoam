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
    monodomain1DSolver

Description
    Implicit monodomain cable-equation solver for a tree conduction-system graph.

Author
    Simao Nieto de Castro. All rights reserved.

SourceFiles
    monodomain1DSolver.C
\*---------------------------------------------------------------------------*/

#ifndef monodomain1DSolver_H
#define monodomain1DSolver_H

#include "conductionSystemSolver.H"
#include "fvMesh.H"
#include "overrideTypeName.H"
#include "scalarField.H"
#include "scalarMatrices.H"

namespace Foam
{

class monodomain1DSolver
:
    public conductionSystemSolver
{
    mutable scalarField appliedCurrentBuffer_;

    // Vm threshold above which a node is considered activated.
    // Defaults to 0.0 (physical-unit models: rest≈-0.084 V, peak≈+0.05 V).
    scalar activationVmThreshold_;

public:

    OverrideTypeName("monodomain1DSolver");

    monodomain1DSolver(const fvMesh&, const dictionary& dict)
    :
        activationVmThreshold_(dict.lookupOrDefault<scalar>("activationVmThreshold", 0.0))
    {}

    virtual void advance
    (
        conductionSystemDomain& domain,
        scalar t0,
        scalar dt
    ) override;

    virtual bool writesGraphVm() const override
    {
        return true;
    }

    virtual bool writesGraphIion() const override
    {
        return true;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
