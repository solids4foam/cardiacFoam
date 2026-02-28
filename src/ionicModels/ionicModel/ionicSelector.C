#include "ionicSelector.H"
#include "error.H"
#include "IOstreams.H"

namespace Foam {


    label ionicSelector::selectTissueOrDimension
    (
        const dictionary& dict,
        const bool hasManufacturedSolution,
        const List<word>& supportedTissues,
        const List<word>& supportedDimensions
    )
    {
        if (hasManufacturedSolution)
        {
            if (!dict.found("dimension"))
            {
                FatalErrorInFunction
                    << "Manufactured model requires 'dimension' entry"
                    << exit(FatalError);
            }

            word d; dict.lookup("dimension") >> d;

            if (!supportedDimensions.contains(d))
            {
                FatalErrorInFunction
                    << "Unsupported dimension '" << d << "' allowed: "
                    << supportedDimensions
                    << exit(FatalError);
            }

            Info<< "MS: using dimension '" << d << "'" << nl;
            return dimensionFlag(d);
        }
        else
        {
            if (!dict.found("tissue"))
            {
                FatalErrorInFunction
                    << "Normal model requires 'tissue' entry"
                    << exit(FatalError);
            }

            word t; dict.lookup("tissue") >> t;

            if (!supportedTissues.contains(t))
            {
                FatalErrorInFunction
                    << "Unsupported tissue '" << t << "' allowed: "
                    << supportedTissues
                    << exit(FatalError);
            }

            //Info<< "Using tissue '" << t << "'" << nl;

            return tissueFlag(t);
        }
    }

    // Map tissues
    label ionicSelector::tissueFlag(const word& t)
    {
        return (t == "epicardialCells")  ? 1 :
            (t == "mCells")           ? 2 :
            (t == "endocardialCells") ? 3 :
            (t == "myocyte")          ? 4 : -1;
    }

    // Map dimensions
    label ionicSelector::dimensionFlag(const word& d)
    {
        return (d == "1D") ? 1 :
            (d == "2D") ? 2 :
            (d == "3D") ? 3 : -1;
    }

} // namespace Foam
