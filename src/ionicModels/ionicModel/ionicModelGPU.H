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
    Foam::ionicModelGPU

Description
    GPU-oriented ionic-model layer that preserves the existing solver-facing
    ionicModel API while separating direct SoA/device-friendly I/O access from
    host diagnostic/export fallback behavior.
\*---------------------------------------------------------------------------*/

#ifndef ionicModelGPU_H
#define ionicModelGPU_H

#include "ionicModel.H"

namespace Foam
{

class ionicModelGPU
:
    public ionicModel
{
        ionicModelGPU(const ionicModelGPU&) = delete;
        void operator=(const ionicModelGPU&) = delete;

        const ionicModelIO::SelectionPlan& selectedPlanForIO
        (
            const wordList& names,
            ionicModelIO::SelectedMapCache& cache
        ) const
        {
            return ionicModelIO::selectedPlanRef
            (
                names,
                ioStateNames(),
                ioNumStates(),
                ioAlgebraicNames(),
                ioNumAlgebraic(),
                cache
            );
        }

        bool canUseDirectSelectedIO(const wordList& names) const
        {
            return names.size() && !gpuSelectionNeedsAlgebraics(names);
        }

        scalar transformedVmForCellWrite(const label cellI) const
        {
            const auto transformVm = ioVmTransform();

            if (!transformVm)
            {
                return gpuStateValue(cellI, voltageStateIndex());
            }

            scalarField stateValues(ioNumStates(), 0.0);
            gpuGatherCellState(cellI, stateValues);
            return transformVm(stateValues);
        }

        scalar transformedVmForCell
        (
            const label cellI,
            scalarField& stateValues
        ) const
        {
            const auto transformVm = ioVmTransform();

            if (!transformVm)
            {
                return gpuStateValue(cellI, voltageStateIndex());
            }

            gpuGatherCellState(cellI, stateValues);
            return transformVm(stateValues);
        }

        scalar selectedCellValueForIO
        (
            const label cellI,
            const label source,
            const label index
        ) const
        {
            if (source == ionicModelIO::COL_STATE)
            {
                return gpuStateValue(cellI, index);
            }

            if (source == ionicModelIO::COL_RATE)
            {
                return gpuRateValue(cellI, index);
            }

            FatalErrorInFunction
                << "Direct GPU-aware I/O does not support source id "
                << source << " for ionic model " << type() << "."
                << exit(FatalError);

            return 0.0;
        }

protected:

        bool gpuSelectionNeedsAlgebraics(const wordList& names) const
        {
            if (!names.size())
            {
                return false;
            }

            return selectionNeedsAlgebraics(names);
        }

        ionicModelGPU
        (
            const dictionary& dict,
            const label nIntegrationPoints,
            const scalar initialDeltaT,
            const Switch solveVmWithinODESolver = false
        )
        :
            ionicModel
            (
                dict,
                nIntegrationPoints,
                initialDeltaT,
                solveVmWithinODESolver
            )
        {}

        virtual label gpuCellCount() const = 0;
        virtual scalar gpuStateValue(label cellI, label stateI) const = 0;
        virtual scalar gpuRateValue(label cellI, label stateI) const = 0;
        virtual void gpuGatherCellState
        (
            label cellI,
            scalarField& stateValues
        ) const = 0;
        virtual label voltageStateIndex() const = 0;

public:

        virtual void write(const scalar t, OFstream& os) const override
        {
            const wordList names = exportedFieldNames();

            if (!canUseDirectSelectedIO(names))
            {
                ionicModel::write(t, os);
                return;
            }

            const ionicModelIO::SelectionPlan& plan =
                selectedPlanForIO(names, writeSelectedPlanCache_);
            scalar transformedVm = 0.0;
            bool hasTransformedVm = false;

            os << t;

            forAll(plan.source, k)
            {
                const label source = plan.source[k];
                const label index = plan.index[k];

                if (source == ionicModelIO::COL_VM)
                {
                    if (!hasTransformedVm)
                    {
                        transformedVm = transformedVmForCellWrite(0);
                        hasTransformedVm = true;
                    }

                    os << " " << transformedVm;
                }
                else
                {
                    os << " " << selectedCellValueForIO(0, source, index);
                }
            }

            os << nl;
        }

        virtual void exportStates
        (
            PtrList<volScalarField>& outFields
        ) override
        {
            const wordList names = exportedFieldNames();

            if (!canUseDirectSelectedIO(names))
            {
                ionicModel::exportStates(outFields);
                return;
            }

            const ionicModelIO::SelectionPlan& plan =
                selectedPlanForIO(names, exportSelectedPlanCache_);
            scalarField stateValues(ioNumStates(), 0.0);

            if (outFields.size() != plan.source.size())
            {
                FatalErrorInFunction
                    << "Mismatch between selected export variables ("
                    << plan.source.size() << ") and allocated output fields ("
                    << outFields.size() << ")."
                    << exit(FatalError);
            }

            forAll(outFields, fieldI)
            {
                const label source = plan.source[fieldI];
                const label index = plan.index[fieldI];
                volScalarField& outField = outFields[fieldI];

                for (label cellI = 0; cellI < gpuCellCount(); ++cellI)
                {
                    outField[cellI] =
                        (source == ionicModelIO::COL_VM)
                      ? transformedVmForCell(cellI, stateValues)
                      : selectedCellValueForIO(cellI, source, index);
                }

                outField.correctBoundaryConditions();
            }
        }

        virtual void exportFields
        (
            const wordList& fieldNames,
            PtrList<volScalarField>& outFields
        ) const override
        {
            if (!canUseDirectSelectedIO(fieldNames))
            {
                ionicModel::exportFields(fieldNames, outFields);
                return;
            }

            const ionicModelIO::SelectionPlan& plan =
                selectedPlanForIO(fieldNames, exportTransferSelectedPlanCache_);
            scalarField stateValues(ioNumStates(), 0.0);

            if (outFields.size() != plan.source.size())
            {
                FatalErrorInFunction
                    << "Mismatch between selected export variables ("
                    << plan.source.size() << ") and allocated output fields ("
                    << outFields.size() << ")."
                    << exit(FatalError);
            }

            forAll(outFields, fieldI)
            {
                const label source = plan.source[fieldI];
                const label index = plan.index[fieldI];
                volScalarField& outField = outFields[fieldI];

                for (label cellI = 0; cellI < gpuCellCount(); ++cellI)
                {
                    outField[cellI] =
                        (source == ionicModelIO::COL_VM)
                      ? transformedVmForCell(cellI, stateValues)
                      : selectedCellValueForIO(cellI, source, index);
                }

                outField.correctBoundaryConditions();
            }
        }

        virtual void debugPrintFields
        (
            label cellI,
            scalar t1,
            scalar t2,
            scalar step
        ) const override
        {
            const wordList names = debugPrintedNames();

            if (!names.size())
            {
                return;
            }

            if (!canUseDirectSelectedIO(names))
            {
                ionicModel::debugPrintFields(cellI, t1, t2, step);
                return;
            }

            const ionicModelIO::SelectionPlan& plan =
                selectedPlanForIO(names, debugSelectedPlanCache_);
            scalarField stateValues(ioNumStates(), 0.0);

            Info<< "DEBUG cell=" << cellI
                << " t=" << t1 << "→" << t2;

            if (step >= 0)
            {
                Info<< " step=" << step;
            }

            forAll(names, k)
            {
                const label source = plan.source[k];
                const label index = plan.index[k];

                Info<< " " << names[k]
                    << "="
                    << (
                        source == ionicModelIO::COL_VM
                      ? transformedVmForCell(cellI, stateValues)
                      : selectedCellValueForIO(cellI, source, index)
                    );
            }

            Info<< nl;
        }

        virtual scalar signal
        (
            const label i,
            const CouplingSignal s
        ) const override
        {
            if (i < 0 || i >= gpuCellCount())
            {
                FatalErrorInFunction
                    << "Requested coupling signal index i=" << i
                    << " but valid integration-point range is [0, "
                    << (gpuCellCount() - 1) << "]."
                    << abort(FatalError);
            }

            if (s == CouplingSignal::VM)
            {
                scalarField stateValues(ioNumStates(), 0.0);
                return transformedVmForCell(i, stateValues);
            }

            if (s == CouplingSignal::CAI)
            {
                const label caiIdx =
                    ionicVariableCompatibility::findCaiStateIndex
                    (
                        ioStateNames(),
                        ioNumStates()
                    );

                if (caiIdx >= 0 && caiIdx < ioNumStates())
                {
                    return gpuStateValue(i, caiIdx);
                }
            }

            return ionicModel::signal(i, s);
        }
};

} // End namespace Foam

#endif
