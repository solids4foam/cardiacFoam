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
    Foam::ecgErrorStatistics

Description
    Per-electrode error statistics of a manufactured ECG verifier: the error
    of the numeric trace against a quadrature reference, against each coarser
    quadrature check, and the distance between each check and the reference.
    Each is accumulated as L1 sum, L2 sum of squares and L-infinity over the
    samples, and written as the per-electrode table of the summary file.

    Used by manufacturedPseudoECGVerifier and manufacturedEikonalECGVerifier.
    The lists are sized for the number of checks and electrodes the verifier
    reports, and add() must be given that number of checks.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef ecgErrorStatistics_H
#define ecgErrorStatistics_H

#include "List.H"
#include "Ostream.H"
#include "labelList.H"
#include "scalar.H"
#include "wordList.H"

namespace Foam
{

class ecgErrorStatistics
{
    // Private Data

        //- Samples accumulated
        label sampleCount_;

        //- Numeric against the reference, per electrode
        List<scalar> referenceErrorL1Sum_;
        List<scalar> referenceErrorL2Sum_;
        List<scalar> referenceErrorLinf_;

        //- Numeric against each check, per check and electrode
        List<List<scalar>> checkErrorL1Sum_;
        List<List<scalar>> checkErrorL2Sum_;
        List<List<scalar>> checkErrorLinf_;

        //- Each check against the reference, per check and electrode
        List<List<scalar>> referenceDeltaL1Sum_;
        List<List<scalar>> referenceDeltaL2Sum_;
        List<List<scalar>> referenceDeltaLinf_;


public:

    // Constructors

        //- Statistics of nElectrodes electrodes and no checks yet
        explicit ecgErrorStatistics(const label nElectrodes)
        :
            sampleCount_(0),
            referenceErrorL1Sum_(nElectrodes, scalar(0)),
            referenceErrorL2Sum_(nElectrodes, scalar(0)),
            referenceErrorLinf_(nElectrodes, scalar(0)),
            checkErrorL1Sum_(),
            checkErrorL2Sum_(),
            checkErrorLinf_(),
            referenceDeltaL1Sum_(),
            referenceDeltaL2Sum_(),
            referenceDeltaLinf_()
        {}


    // Member Functions

        //- Samples accumulated
        label sampleCount() const
        {
            return sampleCount_;
        }

        //- Size the lists for nChecks checks and nElectrodes electrodes. New
        //  entries are zero; accumulated entries are kept.
        void resize(const label nChecks, const label nElectrodes)
        {
            referenceErrorL1Sum_.setSize(nElectrodes, scalar(0));
            referenceErrorL2Sum_.setSize(nElectrodes, scalar(0));
            referenceErrorLinf_.setSize(nElectrodes, scalar(0));
            checkErrorL1Sum_.setSize(nChecks);
            checkErrorL2Sum_.setSize(nChecks);
            checkErrorLinf_.setSize(nChecks);
            referenceDeltaL1Sum_.setSize(nChecks);
            referenceDeltaL2Sum_.setSize(nChecks);
            referenceDeltaLinf_.setSize(nChecks);

            for (label checkI = 0; checkI < nChecks; ++checkI)
            {
                checkErrorL1Sum_[checkI].setSize(nElectrodes, scalar(0));
                checkErrorL2Sum_[checkI].setSize(nElectrodes, scalar(0));
                checkErrorLinf_[checkI].setSize(nElectrodes, scalar(0));
                referenceDeltaL1Sum_[checkI].setSize(nElectrodes, scalar(0));
                referenceDeltaL2Sum_[checkI].setSize(nElectrodes, scalar(0));
                referenceDeltaLinf_[checkI].setSize(nElectrodes, scalar(0));
            }
        }

        //- Add one sample: the numeric values, the value of each check and the
        //  reference value, per electrode. nChecks is the number of checks the
        //  lists were sized for.
        void add
        (
            const label nChecks,
            const List<scalar>& numericValues,
            const List<List<scalar>>& checkReferenceValues,
            const List<scalar>& referenceValues
        )
        {
            ++sampleCount_;

            forAll(numericValues, electrodeI)
            {
                const scalar referenceError =
                    Foam::mag
                    (
                        numericValues[electrodeI] - referenceValues[electrodeI]
                    );

                referenceErrorL1Sum_[electrodeI] += referenceError;
                referenceErrorL2Sum_[electrodeI] += referenceError*referenceError;
                referenceErrorLinf_[electrodeI] =
                    max(referenceErrorLinf_[electrodeI], referenceError);

                for (label checkI = 0; checkI < nChecks; ++checkI)
                {
                    const scalar checkError = Foam::mag
                    (
                        numericValues[electrodeI]
                      - checkReferenceValues[checkI][electrodeI]
                    );
                    const scalar referenceDelta = Foam::mag
                    (
                        referenceValues[electrodeI]
                      - checkReferenceValues[checkI][electrodeI]
                    );

                    checkErrorL1Sum_[checkI][electrodeI] += checkError;
                    checkErrorL2Sum_[checkI][electrodeI] += checkError*checkError;
                    checkErrorLinf_[checkI][electrodeI] =
                        max(checkErrorLinf_[checkI][electrodeI], checkError);

                    referenceDeltaL1Sum_[checkI][electrodeI] += referenceDelta;
                    referenceDeltaL2Sum_[checkI][electrodeI] +=
                        referenceDelta*referenceDelta;
                    referenceDeltaLinf_[checkI][electrodeI] =
                        max(referenceDeltaLinf_[checkI][electrodeI], referenceDelta);
                }
            }
        }

        //- The "qChecks" and "qReference" lines of the summary file
        static void writeQuadratureOrders
        (
            Ostream& os,
            const List<label>& checkQuadratureOrders,
            const label referenceQuadratureOrder
        )
        {
            os << "qChecks";
            forAll(checkQuadratureOrders, checkI)
            {
                os << " " << checkQuadratureOrders[checkI];
            }
            os << "\n";
            os << "qReference " << referenceQuadratureOrder << "\n";
        }

        //- The per-electrode table of the summary file: the column header and
        //  one row per electrode, the statistics normalised by the samples
        void writeTable
        (
            Ostream& os,
            const wordList& electrodeNames,
            const List<label>& checkQuadratureOrders
        ) const
        {
            os << "Electrode  L1_err_ref  L2_err_ref  Linf_err_ref";
            forAll(checkQuadratureOrders, checkI)
            {
                const label qCheck = checkQuadratureOrders[checkI];
                os << "  L1_err_q" << qCheck
                   << "  L2_err_q" << qCheck
                   << "  Linf_err_q" << qCheck
                   << "  L1_delta_q" << qCheck << "_ref"
                   << "  L2_delta_q" << qCheck << "_ref"
                   << "  Linf_delta_q" << qCheck << "_ref";
            }
            os << "\n";

            const scalar count = max(scalar(1), scalar(sampleCount_));

            forAll(electrodeNames, electrodeI)
            {
                os << electrodeNames[electrodeI] << " "
                   << referenceErrorL1Sum_[electrodeI]/count << " "
                   << Foam::sqrt(referenceErrorL2Sum_[electrodeI]/count) << " "
                   << referenceErrorLinf_[electrodeI];

                forAll(checkQuadratureOrders, checkI)
                {
                    os << " " << checkErrorL1Sum_[checkI][electrodeI]/count
                       << " "
                       << Foam::sqrt(checkErrorL2Sum_[checkI][electrodeI]/count)
                       << " " << checkErrorLinf_[checkI][electrodeI]
                       << " "
                       << referenceDeltaL1Sum_[checkI][electrodeI]/count
                       << " "
                       << Foam::sqrt(referenceDeltaL2Sum_[checkI][electrodeI]/count)
                       << " " << referenceDeltaLinf_[checkI][electrodeI];
                }

                os << "\n";
            }
        }

        //- Column names of the time-series file: per electrode the numeric
        //  value, the value of each check and of the reference, their errors
        //  against the reference, and each check's distance to the reference
        static wordList seriesColumns
        (
            const wordList& electrodeNames,
            const List<label>& checkQuadratureOrders,
            const label referenceQuadratureOrder
        )
        {
            wordList columns;

            forAll(electrodeNames, electrodeI)
            {
                const word& name = electrodeNames[electrodeI];
                columns.append("numeric_" + name);

                forAll(checkQuadratureOrders, checkI)
                {
                    columns.append
                    (
                        "refQ" + Foam::name(checkQuadratureOrders[checkI])
                      + "_" + name
                    );
                }

                columns.append
                (
                    "refQ" + Foam::name(referenceQuadratureOrder) + "_" + name
                );

                forAll(checkQuadratureOrders, checkI)
                {
                    columns.append
                    (
                        "errQ" + Foam::name(checkQuadratureOrders[checkI])
                      + "_" + name
                    );
                }

                columns.append
                (
                    "errQ" + Foam::name(referenceQuadratureOrder) + "_" + name
                );

                forAll(checkQuadratureOrders, checkI)
                {
                    columns.append
                    (
                        "deltaQuadratureQ"
                      + Foam::name(checkQuadratureOrders[checkI])
                      + "_Q"
                      + Foam::name(referenceQuadratureOrder)
                      + "_"
                      + name
                    );
                }
            }

            return columns;
        }
};

} // End namespace Foam

#endif

// ************************************************************************* //
