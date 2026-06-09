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

Description
    Hardcoded restitution templates for 1D Purkinje models.
\*---------------------------------------------------------------------------*/

#ifndef restitutionTemplates_H
#define restitutionTemplates_H

#include "scalar.H"

namespace Foam
{
namespace restitutionTemplates
{

static const int purkinjeSize = 6;

static const scalar purkinjeDI[] = {
    0.01, 0.05, 0.10, 0.20, 0.50, 1.00
};

static const scalar purkinjeAPD[] = {
    0.100, 0.150, 0.200, 0.250, 0.300, 0.320
};

static const scalar purkinjeCV[] = {
    0.50, 1.00, 1.50, 2.00, 2.50, 3.00
};

} // namespace restitutionTemplates
} // namespace Foam

#endif
