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

#include "vtkUnstructuredReader.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkUnstructuredReader::printFieldStats(const objectRegistry& obj)
{
    const UPtrList<const Type> fields(obj.csorted<Type>());

    if (!fields.empty())
    {
        Info<< "Read " << fields.size() << ' ' << Type::typeName
            << " fields:" << nl
            << "Size\tName" << nl
            << "----\t----" << nl;

        for (const Type& field : fields)
        {
            Info<< field.size() << '\t' << field.name() << nl;
        }
        Info<< endl;
    }
}


// ************************************************************************* //
