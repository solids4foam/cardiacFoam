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
    Foam::BatchedTensionExecutor

Description
    Reusable batched-kernel execution layer for SoA active tension stepping.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef batchedActiveTensionExecution_H
#define batchedActiveTensionExecution_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Foam
{

template<class FieldBuffer>
struct BatchedTensionScratchT
{
    FieldBuffer stateValues;
    FieldBuffer rateValues;
    FieldBuffer algebraicValues;

    BatchedTensionScratchT
    (
        const label nStates,
        const label nAlgebraics
    )
    :
        stateValues(nStates, 0.0),
        rateValues(nStates, 0.0),
        algebraicValues(nAlgebraics, 0.0)
    {}

    void resetPrimary()
    {
        rateValues = 0.0;
        algebraicValues = 0.0;
    }
};

template<class Backend, class CellScratch>
class BatchedTensionExecutor
{
    const Backend& backend_;
    const bool parallelCellLoops_;
    PtrList<CellScratch>& threadScratch_;

    static label cellLoopThreadIndex()
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

public:

    BatchedTensionExecutor
    (
        const Backend& backend,
        const bool parallelCellLoops,
        PtrList<CellScratch>& threadScratch
    )
    :
        backend_(backend),
        parallelCellLoops_(parallelCellLoops),
        threadScratch_(threadScratch)
    {}

    template<class DriveSignalBuffer, class LambdaBuffer, class TaBuffer>
    void solveCells
    (
        const scalar t,
        const scalar dt,
        const DriveSignalBuffer& driveSignal,
        const LambdaBuffer& lambda,
        TaBuffer& Ta
    ) const
    {
        const label nSubsteps = backend_.nSubsteps();
        const scalar timeScale = backend_.timeScaleFactor();
        const scalar tModel = t*timeScale;
        const scalar dtSubstep = (dt / scalar(nSubsteps)) * timeScale;

#ifdef _OPENMP
#pragma omp parallel for if (parallelCellLoops_)
#endif
        for (label cellI = 0; cellI < backend_.nCells(); ++cellI)
        {
            CellScratch& scratch = threadScratch_[cellLoopThreadIndex()];
            backend_.gatherCellState(cellI, scratch.stateValues);

            const scalar cellDrive = driveSignal[cellI];
            const scalar cellLambda = lambda[cellI];

            for (label substep = 0; substep < nSubsteps; ++substep)
            {
                const scalar substepTime = tModel + substep * dtSubstep;

                scratch.resetPrimary();

                backend_.evaluateScratchAtTime
                (
                    cellI,
                    substepTime,
                    cellDrive,
                    cellLambda,
                    scratch
                );

                backend_.advanceSubstep(cellI, dtSubstep, scratch);
            }

            scratch.resetPrimary();

            backend_.evaluateScratchAtTime
            (
                cellI,
                tModel + scalar(nSubsteps) * dtSubstep,
                cellDrive,
                cellLambda,
                scratch
            );

            backend_.scatterCellState(cellI, scratch.stateValues);
            backend_.syncEvaluatedOutputs(cellI, scratch);

            // Extract the active tension and output it
            Ta[cellI] = backend_.activeTensionFromScratch(cellI, scratch);
        }
    }
};

} // End namespace Foam

#endif
