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

Class
    prePacingIO

Description
    Static utility class for reading the optional, opt-in
    constant/prePacingProperties dictionary. Absence of the file disables
    prePacing everywhere.
\*---------------------------------------------------------------------------*/

#ifndef prePacingIO_H
#define prePacingIO_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "word.H"

namespace Foam
{
namespace prePacingIO
{

struct PrePacingConfig
{
    bool enabled = false;
    scalar tolerance = 1e-4;
    label minBeats = 10;
    label maxBeats = 3000;
    scalar autorhythmicCheckInterval = 200.0;   // ms
    bool autorhythmicCheckIntervalSet = false;

    //- singleCellStimulus-format S1 protocol (ms) used while pre-pacing
    bool hasStimulus = false;
    dictionary stimulus;
};

inline PrePacingConfig configFor(const fvMesh& mesh, const word& regionName)
{
    PrePacingConfig cfg;

    IOobject dictIO
    (
        "prePacingProperties",
        mesh.time().constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (!dictIO.typeHeaderOk<IOdictionary>(true))
    {
        return cfg;
    }

    const IOdictionary dict(dictIO);

    cfg.enabled = true;
    cfg.tolerance = dict.lookupOrDefault<scalar>("tolerance", cfg.tolerance);
    cfg.minBeats = dict.lookupOrDefault<label>("minBeats", cfg.minBeats);
    cfg.maxBeats = dict.lookupOrDefault<label>("maxBeats", cfg.maxBeats);
    cfg.autorhythmicCheckIntervalSet = dict.found("autorhythmicCheckInterval");
    cfg.autorhythmicCheckInterval =
        dict.lookupOrDefault<scalar>
        (
            "autorhythmicCheckInterval", cfg.autorhythmicCheckInterval
        );
    if (dict.isDict("stimulus"))
    {
        cfg.hasStimulus = true;
        cfg.stimulus = dict.subDict("stimulus");
    }

    if (dict.found("regions"))
    {
        const dictionary& regionsDict = dict.subDict("regions");
        if (regionsDict.found(regionName))
        {
            const dictionary& regionDict = regionsDict.subDict(regionName);
            cfg.enabled =
                regionDict.lookupOrDefault<bool>("enabled", cfg.enabled);
            cfg.tolerance =
                regionDict.lookupOrDefault<scalar>
                (
                    "tolerance", cfg.tolerance
                );
            cfg.minBeats =
                regionDict.lookupOrDefault<label>("minBeats", cfg.minBeats);
            cfg.maxBeats =
                regionDict.lookupOrDefault<label>("maxBeats", cfg.maxBeats);
            if (regionDict.found("autorhythmicCheckInterval"))
            {
                cfg.autorhythmicCheckIntervalSet = true;
                cfg.autorhythmicCheckInterval =
                    regionDict.get<scalar>("autorhythmicCheckInterval");
            }
            if (regionDict.isDict("stimulus"))
            {
                cfg.hasStimulus = true;
                cfg.stimulus = regionDict.subDict("stimulus");
            }
        }
    }

    return cfg;
}

} // End namespace prePacingIO
} // End namespace Foam

#endif

// ************************************************************************* //
