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

#include "heartPhiETrace.H"
#include "torsoCoupling.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::Enum<Foam::heartPhiETrace::type> Foam::heartPhiETrace::names
({
    { type::zeroGradient, "zeroGradient" },
    { type::global, "global" }
});


Foam::heartPhiETrace::type Foam::heartPhiETrace::readWithTorso
(
    const dictionary& myocardiumCoeffs
)
{
    const type trace = names.get("bathHeartPhiETrace", myocardiumCoeffs);

    if
    (
        trace == type::global
     && !myocardiumCoeffs.get<Switch>("sealedHeartBoundary")
    )
    {
        FatalIOErrorInFunction(myocardiumCoeffs)
            << "bathHeartPhiETrace global requires sealedHeartBoundary true"
            << exit(FatalIOError);
    }

    return trace;
}


Foam::heartPhiETrace::type Foam::heartPhiETrace::read
(
    const dictionary& myocardiumCoeffs
)
{
    if
    (
        torsoCoupling::readFromCoeffs(myocardiumCoeffs)
     != torsoCoupling::type::none
    )
    {
        return readWithTorso(myocardiumCoeffs);
    }

    if (myocardiumCoeffs.found("bathHeartPhiETrace"))
    {
        FatalIOErrorInFunction(myocardiumCoeffs)
            << "bathHeartPhiETrace sets the heart phiE at the heart-torso "
            << "wall, but electroProperties has no torso"
            << exit(FatalIOError);
    }

    // Without a torso the heart wall is the outer boundary, where the heart
    // phiE is zeroGradient
    return type::zeroGradient;
}


// ************************************************************************* //
