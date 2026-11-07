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
    ecgSolver

Description
    Abstract base class for ECG solvers.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef ecgSolver_H
#define ecgSolver_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"
#include "word.H"

namespace Foam
{

class dictionary;
class ecgDomain;
class myocardiumDomainInterface;

class ecgSolver
{
public:

    //- The state an ECG solver reads
    enum class stateSource
    {
        //- The myocardium domain (Vm, activation time, conductivities)
        myocardium,

        //- The bath domain's global phiE (torsoECG with torsoCoupling bath)
        bath
    };

    TypeName("ecgSolver");

    declareRunTimeSelectionTable
    (
        autoPtr,
        ecgSolver,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    ecgSolver() = default;

    //- The ecgSolver of an ECG domain entry (required): the single reader
    //  of the key
    static word selectedType(const dictionary& dict);

    static autoPtr<ecgSolver> New(const dictionary& dict);

    virtual ~ecgSolver() = default;

    virtual bool writesOwnTimeSeries() const
    {
        return false;
    }

    //- The state this solver reads
    virtual stateSource source() const
    {
        return stateSource::myocardium;
    }

    //- Fatal when this myocardium cannot drive the solver (default: any can)
    virtual void checkMyocardium
    (
        const myocardiumDomainInterface& myocardium,
        const word& domainName
    ) const
    {
        (void)myocardium;
        (void)domainName;
    }

    //- Fatal when the domain's state lacks a channel this solver reads
    //  (default: none required at construction)
    virtual void checkState(const ecgDomain& domain) const
    {
        (void)domain;
    }

    virtual void solve
    (
        ecgDomain& domain,
        scalar t0,
        scalar dt,
        scalarField& values
    ) = 0;
};

} // End namespace Foam

#endif

// ************************************************************************* //
