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

#include "dictionaryKeyCheck.H"
#include "DynamicList.H"
#include "FlatOutput.H"

#include <algorithm>
#include <cctype>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

// Edit distance between a and b, ignoring case
std::size_t editDistance(const std::string& a, const std::string& b)
{
    std::vector<std::size_t> prev(b.size() + 1), curr(b.size() + 1);

    for (std::size_t j = 0; j <= b.size(); ++j)
    {
        prev[j] = j;
    }

    for (std::size_t i = 1; i <= a.size(); ++i)
    {
        curr[0] = i;

        for (std::size_t j = 1; j <= b.size(); ++j)
        {
            const bool same =
                std::tolower(static_cast<unsigned char>(a[i - 1]))
             == std::tolower(static_cast<unsigned char>(b[j - 1]));

            curr[j] = std::min
            ({
                prev[j] + 1,
                curr[j - 1] + 1,
                prev[j - 1] + (same ? 0 : 1)
            });
        }

        std::swap(prev, curr);
    }

    return prev[b.size()];
}


// The literal (non-regex) entries of allowed, sorted
Foam::wordList literalNames(const Foam::wordRes& allowed)
{
    Foam::DynamicList<Foam::word> names(allowed.size());

    for (const Foam::wordRe& w : allowed)
    {
        if (w.isLiteral())
        {
            names.append(Foam::word(w));
        }
    }

    Foam::wordList result(std::move(names));
    Foam::sort(result);
    return result;
}


void reportUnknown
(
    const Foam::dictionary& dict,
    const Foam::wordList& unknown,
    const Foam::wordRes& allowed
)
{
    using namespace Foam;

    const wordList valid(literalNames(allowed));

    FatalIOErrorIn("dictionaryKeyCheck::checkKeys", dict)
        << "Unknown keyword" << (unknown.size() > 1 ? "s" : "")
        << " in dictionary " << dict.relativeName() << ':' << nl;

    for (const word& key : unknown)
    {
        const word suggestion(dictionaryKeyCheck::closestMatch(key, valid));

        FatalIOError << "    '" << key << "'";

        if (!suggestion.empty())
        {
            FatalIOError << ", did you mean '" << suggestion << "'?";
        }

        FatalIOError << nl;
    }

    FatalIOError << nl << "Valid keywords are " << flatOutput(valid);

    if (valid.size() < allowed.size())
    {
        FatalIOError << ", plus keywords matching:";

        for (const wordRe& w : allowed)
        {
            if (!w.isLiteral())
            {
                FatalIOError << ' ' << w;
            }
        }
    }

    FatalIOError
        << nl << "Keywords starting with '_' are case variables and are "
        << "always accepted." << exit(FatalIOError);
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::dictionaryKeyCheck::closestMatch
(
    const word& name,
    const wordList& candidates
)
{
    // A likely misspelling: a few edits, at most a third of the name
    const std::size_t maxDistance =
        std::max<std::size_t>(1, std::min<std::size_t>(3, name.size()/3));

    word best;
    std::size_t bestDistance = maxDistance + 1;

    for (const word& candidate : candidates)
    {
        const std::size_t d = editDistance(name, candidate);

        if (d < bestDistance)
        {
            best = candidate;
            bestDistance = d;
        }
    }

    // Otherwise a shortened name, e.g. Icoupling for IcouplingSource
    if (best.empty() && name.size() >= 3)
    {
        for (const word& candidate : candidates)
        {
            if (candidate.starts_with(name))
            {
                return candidate;
            }
        }
    }

    return best;
}


void Foam::dictionaryKeyCheck::checkKeys
(
    const dictionary& dict,
    const wordRes& allowed
)
{
    DynamicList<word> unknown;

    for (const entry& e : dict)
    {
        const word& key = e.keyword();

        if (key.starts_with('_') || allowed.match(key))
        {
            continue;
        }

        unknown.append(key);
    }

    if (!unknown.empty())
    {
        reportUnknown(dict, unknown, allowed);
    }
}


void Foam::dictionaryKeyCheck::checkNamedBlocks
(
    const dictionary& dict,
    const wordRes& allowed,
    const wordRes& allowedEntries,
    const wordList& skipped
)
{
    DynamicList<word> unknown;

    for (const entry& e : dict)
    {
        const word& key = e.keyword();

        if (key.starts_with('_') || skipped.found(key))
        {
            continue;
        }

        if (e.isDict())
        {
            checkKeys(e.dict(), allowed);
        }
        else if (!allowedEntries.match(key))
        {
            unknown.append(key);
        }
    }

    if (!unknown.empty())
    {
        FatalIOErrorInFunction(dict)
            << "Dictionary " << dict.relativeName()
            << " holds named sub-dictionaries, but " << unknown
            << (unknown.size() > 1 ? " are not dictionaries" : " is not a dictionary")
            << '.' << exit(FatalIOError);
    }
}


void Foam::dictionaryKeyCheck::checkValue
(
    const dictionary& dict,
    const word& keyword,
    const word& value,
    const wordList& valid
)
{
    if (valid.found(value))
    {
        return;
    }

    const word suggestion(closestMatch(value, valid));

    FatalIOErrorInFunction(dict)
        << "Unknown " << keyword << " '" << value << "' in dictionary "
        << dict.relativeName();

    if (!suggestion.empty())
    {
        FatalIOError << ", did you mean '" << suggestion << "'?";
    }

    FatalIOError
        << nl << "Valid " << keyword << " values are " << flatOutput(valid)
        << exit(FatalIOError);
}


Foam::word Foam::dictionaryKeyCheck::getWord
(
    const dictionary& dict,
    const word& keyword,
    const word& deflt,
    const wordList& valid
)
{
    const word value(dict.getOrDefault<word>(keyword, deflt));
    checkValue(dict, keyword, value, valid);
    return value;
}


// ************************************************************************* //
