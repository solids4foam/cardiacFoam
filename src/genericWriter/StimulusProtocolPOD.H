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
    Plain-Old-Data (POD) projection of `Foam::StimulusProtocol`, designed to
    be safely passed to and used inside batched / device-side ionic-model
    evaluators.

    Why this exists
    ---------------
    `Foam::StimulusProtocol` (in stimulusIO.H) pulls in `scalarField.H`,
    `dictionary.H`, `OFstream.H`, `boundBox.H`, etc. - none of which are
    device-callable and many of which are `<iostream>`/<fstream>-heavy.
    A `__device__` evaluator cannot include any of them.

    `StimulusProtocolPOD` is a flat, trivially-copyable, zero-init-able
    struct containing only `double` and `int` members. The companion
    `computeStimulusPOD()` is a header-only inline (annotated
    `__host__ __device__` when compiling under CUDA/HIP) that mirrors the
    behavior of `Foam::stimulusIO::computeStimulus()`.

    To pass a stimulus from host to device, build a `StimulusProtocolPOD`
    on the host (use `Foam::stimulusIO::toPOD()` from stimulusIO.H), then
    pass-by-value into the kernel - no pointer chasing, no virtual calls,
    no Foam types crossing the device boundary.
\*---------------------------------------------------------------------------*/

#ifndef StimulusProtocolPOD_H
#define StimulusProtocolPOD_H

#include <cmath>

// Host/device qualifier - empty on host-only builds, expanded under CUDA
// or HIP. Centralised here so the macro can be redefined or undefined by
// callers if needed.
#if defined(__CUDACC__) || defined(__HIPCC__)
    #define CARDIAC_HOST_DEVICE __host__ __device__
#else
    #define CARDIAC_HOST_DEVICE
#endif

namespace Foam
{

//- Flat, trivially-copyable mirror of Foam::StimulusProtocol. POD-only:
//  safe to memcpy, safe to pass by value into a CUDA/HIP kernel, safe to
//  store in unified memory.
struct StimulusProtocolPOD
{
    double stimStart;     // (ms) onset of S1 train
    double stimPeriodS1;  // (ms) S1-S1 cycle length
    double stimDuration;  // (ms) pulse width
    double stimAmplitude; // pulse amplitude (sign convention: applied as -amp)
    double stimPeriodS2;  // (ms) S2-S2 cycle length (optional)
    int    nStim1;        // number of S1 pulses
    int    nStim2;        // number of S2 pulses
    int    _pad0;         // reserved, keep 8-byte alignment
};

static_assert
(
    sizeof(StimulusProtocolPOD) == 5*sizeof(double) + 4*sizeof(int),
    "StimulusProtocolPOD size unexpected - check padding"
);

//- Bit-equivalent of `Foam::stimulusIO::computeStimulus(t, stim)`,
//  callable from host or device code. Pure function; no globals, no
//  Foam types.
//
//  NOTE: keep this implementation in lock-step with
//        `stimulusIO.C::computeStimulus(...)`. Any logic change there
//        must be mirrored here, otherwise host-vs-device runs will
//        diverge.
CARDIAC_HOST_DEVICE inline double computeStimulusPOD
(
    const double t,
    const StimulusProtocolPOD& s
)
{
    double Istim = 0.0;

    // ---------- S1 ----------
    if (s.stimPeriodS1 > 0.0 && s.nStim1 > 0)
    {
        const double tp = t - s.stimStart;
        const double trainEnd = s.stimPeriodS1 * static_cast<double>(s.nStim1);

        if (tp >= 0.0 && tp <= trainEnd)
        {
            const double phase =
                tp - std::floor(tp / s.stimPeriodS1) * s.stimPeriodS1;

            if (phase >= 0.0 && phase <= s.stimDuration)
            {
                Istim = -s.stimAmplitude;
            }
        }
    }

    // ---------- S2 ----------
    if (Istim == 0.0 && s.stimPeriodS2 > 0.0 && s.nStim2 > 0)
    {
        const double tS1End =
            s.stimStart + s.stimPeriodS1 * static_cast<double>(s.nStim1);
        const double tp = t - tS1End;
        const double trainEnd = s.stimPeriodS2 * static_cast<double>(s.nStim2);

        if (tp >= 0.0 && tp <= trainEnd)
        {
            const double phase =
                tp - std::floor(tp / s.stimPeriodS2) * s.stimPeriodS2;

            if (phase >= 0.0 && phase <= s.stimDuration)
            {
                Istim = -s.stimAmplitude;
            }
        }
    }

    return Istim;
}

} // End namespace Foam

#endif

// ************************************************************************* //
