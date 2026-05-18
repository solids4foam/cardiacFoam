/*---------------------------------------------------------------------------*\
License
    This file is part of cardiacFoam.

    cardiacFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    cardiacFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    State / algebraic / constant indexing for the bath-bidomain variant of
    the FDA manufactured solution. Mirrors the heart-only variant; kept in
    a separate header so the two ionic models never share linkage-level
    symbols and can diverge independently.
\*---------------------------------------------------------------------------*/

#ifndef bathBidomainFDAManufactured_2014Names_H
#define bathBidomainFDAManufactured_2014Names_H

#define BATH_BIDOMAIN_NUM_STATES 4
#define BATH_BIDOMAIN_NUM_ALGEBRAIC 1
#define BATH_BIDOMAIN_NUM_CONSTANTS 3

enum BathBidomainStates
{
    BathBidomainV = 0,
    BathBidomainU1,
    BathBidomainU2,
    BathBidomainU3
};

enum BathBidomainAlgebraic
{
    BathBidomainIion = 0
};

enum BathBidomainConstants
{
    BathBidomainCm = 0,
    BathBidomainBeta,
    BathBidomainChi
};

#endif

// ************************************************************************* //
