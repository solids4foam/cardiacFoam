#ifndef ionicSelector_H
#define ionicSelector_H

#include "List.H"
#include "dictionary.H"
#include "word.H"

namespace Foam {

class ionicSelector {
public:
  // Main interface
  static label selectTissue(const dictionary &dict,
                            const List<word> &supportedTissues);

  static label tissueFlag(const word &t);
  static label dimensionFlag(const word &d);
};

} // End namespace Foam

#endif
