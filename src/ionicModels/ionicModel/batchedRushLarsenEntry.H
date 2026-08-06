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
    Shared dispatch-table entry shape used by every batched ionic model's
    Rush-Larsen extraction. One entry per state. The entry carries fields
    for the support API (`rushLarsenParametersFromHotPathSupport`, reads
    supportValues). All models use the compact/hot-path support path.

    Each model's dispatch table is:
        std::array<batchedRushLarsenEntry, NUM_STATES>

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef batchedRushLarsenEntry_H
#define batchedRushLarsenEntry_H

#include <cmath>
#include <cstdint>

namespace Foam
{

enum class RLScalarSource : std::uint8_t
{
    none = 0,
    algebraic,   // algebraicValues[index]
    constant,    // CONSTANTS_[index]
    literal      // literal field
};

enum class RLSupportSource : std::uint8_t
{
    none = 0,
    support      // supportValues[supportIndex]
};

struct batchedRushLarsenEntry
{
    RLScalarSource scalarTauSource    = RLScalarSource::none;
    int            scalarTauIndex     = -1;
    double         scalarTauLiteral   = 0.0;

    RLScalarSource scalarGInfSource   = RLScalarSource::none;
    int            scalarGInfIndex    = -1;
    double         scalarGInfLiteral  = 0.0;

    RLSupportSource supportTauSource    = RLSupportSource::none;
    int             supportTauIndex     = -1;

    RLSupportSource supportGInfSource   = RLSupportSource::none;
    int             supportGInfIndex    = -1;
};


inline constexpr batchedRushLarsenEntry rlNone() { return {}; }

inline constexpr batchedRushLarsenEntry rlSupport
(
    int tauSupportIdx, int gInfSupportIdx
)
{
    batchedRushLarsenEntry e;
    e.supportTauSource  = RLSupportSource::support;
    e.supportTauIndex   = tauSupportIdx;
    e.supportGInfSource = RLSupportSource::support;
    e.supportGInfIndex  = gInfSupportIdx;
    return e;
}

inline constexpr batchedRushLarsenEntry rlScalarAlgAndSupport
(
    int scalarTauAlg, int scalarGInfAlg,
    int supportTau,   int supportGInf
)
{
    batchedRushLarsenEntry e;
    e.scalarTauSource   = RLScalarSource::algebraic;
    e.scalarTauIndex    = scalarTauAlg;
    e.scalarGInfSource  = RLScalarSource::algebraic;
    e.scalarGInfIndex   = scalarGInfAlg;
    e.supportTauSource  = RLSupportSource::support;
    e.supportTauIndex   = supportTau;
    e.supportGInfSource = RLSupportSource::support;
    e.supportGInfIndex  = supportGInf;
    return e;
}

inline constexpr batchedRushLarsenEntry rlScalarConstTauAlgInfAndSupport
(
    int scalarTauConst, int scalarGInfAlg,
    int supportTau, int supportGInf
)
{
    batchedRushLarsenEntry e;
    e.scalarTauSource   = RLScalarSource::constant;
    e.scalarTauIndex    = scalarTauConst;
    e.scalarGInfSource  = RLScalarSource::algebraic;
    e.scalarGInfIndex   = scalarGInfAlg;
    e.supportTauSource  = RLSupportSource::support;
    e.supportTauIndex   = supportTau;
    e.supportGInfSource = RLSupportSource::support;
    e.supportGInfIndex  = supportGInf;
    return e;
}

inline constexpr batchedRushLarsenEntry rlScalarLitTauAlgInfAndSupport
(
    double scalarTauLit, int scalarGInfAlg,
    int supportTau, int supportGInf
)
{
    batchedRushLarsenEntry e;
    e.scalarTauSource   = RLScalarSource::literal;
    e.scalarTauLiteral  = scalarTauLit;
    e.scalarGInfSource  = RLScalarSource::algebraic;
    e.scalarGInfIndex   = scalarGInfAlg;
    e.supportTauSource  = RLSupportSource::support;
    e.supportTauIndex   = supportTau;
    e.supportGInfSource = RLSupportSource::support;
    e.supportGInfIndex  = supportGInf;
    return e;
}

template<class Constants, class SupportValues, class Scalar>
inline bool resolveSupportRushLarsenEntry
(
    const batchedRushLarsenEntry& entry,
    const Constants&,
    const SupportValues& supportValues,
    const Scalar minTau,
    Scalar& steadyState,
    Scalar& tau
)
{
    if (entry.supportTauSource == RLSupportSource::none) return false;

    switch (entry.supportTauSource)
    {
        case RLSupportSource::support:
            tau = supportValues[entry.supportTauIndex];
            break;
        case RLSupportSource::none:
            return false;
    }

    if (tau <= minTau || !std::isfinite(tau)) return false;

    switch (entry.supportGInfSource)
    {
        case RLSupportSource::support:
            steadyState = supportValues[entry.supportGInfIndex];
            break;
        case RLSupportSource::none:
            return false;
    }

    return std::isfinite(steadyState);
}

template<class Constants, class Algebraics, class SupportValues>
inline void projectScalarRushLarsenEntryToSupport
(
    const batchedRushLarsenEntry& entry,
    const Constants& constants,
    const Algebraics& algebraicValues,
    SupportValues& supportValues
)
{
    if (entry.supportTauSource == RLSupportSource::support)
    {
        switch (entry.scalarTauSource)
        {
            case RLScalarSource::algebraic:
                supportValues[entry.supportTauIndex] =
                    algebraicValues[entry.scalarTauIndex];
                break;
            case RLScalarSource::constant:
                supportValues[entry.supportTauIndex] =
                    constants[entry.scalarTauIndex];
                break;
            case RLScalarSource::literal:
                supportValues[entry.supportTauIndex] =
                    entry.scalarTauLiteral;
                break;
            case RLScalarSource::none:
                break;
        }
    }

    if (entry.supportGInfSource == RLSupportSource::support)
    {
        switch (entry.scalarGInfSource)
        {
            case RLScalarSource::algebraic:
                supportValues[entry.supportGInfIndex] =
                    algebraicValues[entry.scalarGInfIndex];
                break;
            case RLScalarSource::constant:
                supportValues[entry.supportGInfIndex] =
                    constants[entry.scalarGInfIndex];
                break;
            case RLScalarSource::literal:
                supportValues[entry.supportGInfIndex] =
                    entry.scalarGInfLiteral;
                break;
            case RLScalarSource::none:
                break;
        }
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //
