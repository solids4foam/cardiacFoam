/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/

#ifndef solverStartupLog_H
#define solverStartupLog_H

#include "IOstreams.H"
#include "word.H"

namespace Foam
{
namespace solverStartupLog
{

inline void readingDictionary(const word& dictName)
{
    Info<< "Reading " << dictName << nl << endl;
}

inline void sectionDivider()
{
    Info<< "---------------------------------------------" << nl;
}

inline void logSolverMode(const bool solveExplicit)
{
    Info<< nl
        << "Using " << (solveExplicit ? "explicit" : "implicit")
        << " solver" << nl << endl;
}

} // End namespace solverStartupLog
} // End namespace Foam

#endif

