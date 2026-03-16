#include "IOstreams.H"
#include "error.H"
#include "ionicSelector.H"

namespace Foam {

label ionicSelector::selectTissue(const dictionary &dict,
                                  const List<word> &supportedTissues) {
  if (!dict.found("tissue")) {
    FatalErrorInFunction << "Normal model requires 'tissue' entry"
                         << exit(FatalError);
  }

  word t;
  dict.lookup("tissue") >> t;

  if (!supportedTissues.contains(t)) {
    FatalErrorInFunction << "Unsupported tissue '" << t
                         << "' allowed: " << supportedTissues
                         << exit(FatalError);
  }

  return tissueFlag(t);
}

// Map tissues
label ionicSelector::tissueFlag(const word &t) {
  return (t == "epicardialCells")    ? 1
         : (t == "mCells")           ? 2
         : (t == "endocardialCells") ? 3
         : (t == "myocyte")          ? 4
                                     : -1;
}

// Map dimensions
label ionicSelector::dimensionFlag(const word &d) {
  return (d == "1D") ? 1 : (d == "2D") ? 2 : (d == "3D") ? 3 : -1;
}

} // namespace Foam
