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

Namespace
    Foam::dictionaryKeyCheck

Description
    Rejects dictionary keywords and word values that no reader accepts, so a
    misspelt or invented entry stops the run instead of being ignored.

    The error names the entry, the closest valid name and every valid name.
    Keywords that start with '_' are case variables for $-substitution and
    #ifeq; no reader uses them and they are always accepted.

SourceFiles
    dictionaryKeyCheck.C
\*---------------------------------------------------------------------------*/

#ifndef dictionaryKeyCheck_H
#define dictionaryKeyCheck_H

#include "dictionary.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dictionaryKeyCheck
{

//- Stop with a FatalIOError if dict holds a keyword that allowed does not
//  match.
void checkKeys(const dictionary& dict, const wordRes& allowed);

//- checkKeys on every sub-dictionary of dict, for blocks of user-named
//  entries such as ecgDomains. A non-dictionary entry is rejected unless
//  allowedEntries matches it; a sub-dictionary named in skipped is not
//  checked.
void checkNamedBlocks
(
    const dictionary& dict,
    const wordRes& allowed,
    const wordRes& allowedEntries = wordRes(),
    const wordList& skipped = wordList()
);

//- Stop with a FatalIOError if value, read from keyword of dict, is not one
//  of valid.
void checkValue
(
    const dictionary& dict,
    const word& keyword,
    const word& value,
    const wordList& valid
);

//- The word read from keyword (or deflt when keyword is absent), after
//  checkValue.
word getWord
(
    const dictionary& dict,
    const word& keyword,
    const word& deflt,
    const wordList& valid
);

//- The entry of candidates closest to name, or word::null when none is
//  close enough to be a likely misspelling.
word closestMatch(const word& name, const wordList& candidates);

} // End namespace dictionaryKeyCheck
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
