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
    Reader for the optional, opt-in constant/prePacingProperties
    dictionary. Absence of the file disables prePacing everywhere.

    Top-level entries apply to every tissue region; regions.<name>
    overrides them for one region, where <name> is an ionicHeterogeneity
    region name, or the 'tissue' entry when there is no heterogeneity.

        tolerance               1e-4;
        minBeats                10;
        maxBeats                3000;
        singleCellStimulus      { stim_start ...; stim_period_S1 ...; ... }
        beatComparisonInterval  800;    // ms, only without a stimulus
        regions { endocardialCells { maxBeats 5000; } }
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

    //- ms between state comparisons when no singleCellStimulus is given
    scalar beatComparisonInterval = 200.0;
    bool beatComparisonIntervalSet = false;

    //- Pacing protocol (ms), same keys as electroProperties
    bool hasSingleCellStimulus = false;
    dictionary singleCellStimulus;
};

//- Apply the entries of dict on top of cfg.
inline void readEntries(const dictionary& dict, PrePacingConfig& cfg)
{
    cfg.enabled = dict.lookupOrDefault<bool>("enabled", cfg.enabled);
    cfg.tolerance = dict.lookupOrDefault<scalar>("tolerance", cfg.tolerance);
    cfg.minBeats = dict.lookupOrDefault<label>("minBeats", cfg.minBeats);
    cfg.maxBeats = dict.lookupOrDefault<label>("maxBeats", cfg.maxBeats);

    if (dict.found("beatComparisonInterval"))
    {
        cfg.beatComparisonIntervalSet = true;
        cfg.beatComparisonInterval =
            dict.get<scalar>("beatComparisonInterval");
    }

    if (dict.isDict("singleCellStimulus"))
    {
        cfg.hasSingleCellStimulus = true;
        cfg.singleCellStimulus = dict.subDict("singleCellStimulus");
    }
}

//- Settings for one tissue region; regionName may be word::null for the
//  top-level settings only.
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
    readEntries(dict, cfg);

    if (dict.isDict("regions") && dict.subDict("regions").isDict(regionName))
    {
        readEntries(dict.subDict("regions").subDict(regionName), cfg);
    }

    return cfg;
}

} // End namespace prePacingIO
} // End namespace Foam

#endif

// ************************************************************************* //
