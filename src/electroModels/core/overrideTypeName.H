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
    Drop-in replacement for TypeName() in derived classes.

    The standard OpenFOAM TypeName() macro declares
        virtual const word& type() const { return typeName; }
    without the C++11 override specifier, producing
    -Winconsistent-missing-override warnings in derived classes.

    OverrideTypeName() is identical to TypeName() but adds override,
    silencing that warning without touching the OpenFOAM macro infrastructure.

Usage
    In any class that (transitively) overrides a base that also uses TypeName:

        OverrideTypeName("MyDerivedClass");

    The corresponding .C file still uses the normal defineTypeNameAndDebug().

\*---------------------------------------------------------------------------*/

#ifndef overrideTypeName_H
#define overrideTypeName_H

#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Declare a ClassName() with a virtual type() marked override.
//  Use in place of TypeName() in derived classes.
#define OverrideTypeName(TypeNameString)                                       \
    ClassName(TypeNameString);                                                 \
    virtual const ::Foam::word& type() const override { return typeName; }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
