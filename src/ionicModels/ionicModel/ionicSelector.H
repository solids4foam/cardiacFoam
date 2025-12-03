#ifndef ionicSelector_H
#define ionicSelector_H

#include "dictionary.H"
#include "List.H"
#include "word.H"

namespace Foam {

class ionicSelector
{
public:
    // Main interface
    static label selectTissueOrDimension
    (
        const dictionary& dict,
        const bool hasManufacturedSolution,
        const List<word>& supportedTissues,
        const List<word>& supportedDimensions
    );

private:
    static label tissueFlag(const word& t);
    static label dimensionFlag(const word& d);
};

} // End namespace Foam

#endif

