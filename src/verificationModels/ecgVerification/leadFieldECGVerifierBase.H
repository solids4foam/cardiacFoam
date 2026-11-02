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
    leadFieldECGVerifierBase

Description
    The channels and error norms shared by the lead-field ECG verifiers
    (manufacturedLeadFieldECGVerifier: slab and layered;
    manufacturedShellLeadFieldECGVerifier: the spherical shell). At each
    sample and electrode it compares:
    - leadField: the lead-field channel of the torsoECG domain with
      torsoCoupling leadField (z_e^T b(Vm_h) + c_e);
    - direct: a direct forward solve of the same conductor on Vm_h, with the
      same boundary values and reference, assembly and solve repeated until
      the relative change between passes is below the tolerance of the
      fvSolution dictionary leadFieldFixedPoint (as buildLeadFields), sampled
      in the electrode cell;
    - exactPoint, exactCentroid: the exact potential of the reference at the
      electrode and at its cell centroid, for the exact Vm;
    - leadFieldStatic, directStatic: the lead field and the direct solve
      applied to the exact Vm at the heart cell centres;
    - phiStatic: the error of the directStatic potential field against the
      exact potential at the cell centres, over the whole mesh, the heart,
      the torso, the interface band (verificationUtils::interfaceBandCells)
      and the surface (verificationUtils::boundaryAdjacentCells): mean |e|
      and volume-weighted RMS over the volume of the region (not of the
      mesh, as verificationUtils::computeMaskedNorms) and maximum, the
      largest over the samples; the cell of the largest error is logged.

    With a pinned reference the exact potential is shifted to equal the
    reference value at the centroid of the pinned cell.

    Writes postProcessing/manufacturedLeadFieldECG_<domain>.dat (samples)
    and manufacturedLeadFieldECGSummary_<domain>.dat (rows name L1 L2 Linf
    over samples, per electrode and as the per-sample maximum over the
    electrodes).

    The direct solve uses the fvSolution entry phiE.

SourceFiles
    leadFieldECGVerifierBase.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef leadFieldECGVerifierBase_H
#define leadFieldECGVerifierBase_H

#include "ecgVerificationModel.H"
#include "leadFieldECG.H"
#include "OFstream.H"
#include "FixedList.H"

namespace Foam
{

class leadFieldECGVerifierBase
:
    public ecgVerificationModel
{
    //- Running L1, L2 and Linf of an error over the samples
    struct errorNorms
    {
        scalar l1 = 0;
        scalar l2 = 0;
        scalar linf = 0;

        void add(const scalar e)
        {
            l1 += mag(e);
            l2 += sqr(e);
            linf = max(linf, mag(e));
        }
    };

    //- The channels and conductor of the ECG domain, evaluated on any Vm
    autoPtr<leadFieldECG> leadFieldPtr_;

    labelList electrodeCells_;
    List<point> electrodeCentroids_;

    //- Centroid of the pinned cell (gauge of the exact potential)
    point referenceCentroid_;

    //- Direct forward potentials for Vm_h and for the exact Vm, kept as
    //  the start of the next sample's solve
    autoPtr<volScalarField> directPhi_;
    autoPtr<volScalarField> directStaticPhi_;
    autoPtr<volScalarField> heartPhi_;

    autoPtr<OFstream> outputPtr_;

    label sampleCount_;
    scalar exactScale_;

    //- Per quantity: the norms per electrode and of the per-sample maximum
    //  over the electrodes
    List<List<errorNorms>> electrodeNorms_;
    List<errorNorms> maxNorms_;
    errorNorms relativeGapNorms_;

    //- phiStatic: per region (all, heart, torso, interfaceBand, surface), the
    //  largest mean |e|, RMS and maximum over the samples
    static const wordList fieldRegionNames_;
    List<FixedList<scalar, 3>> fieldErrors_;

    bool summaryWritten_;

    static const wordList quantityNames_;

    //- Build the channels and the reference at the first sample
    void initialise();

    //- Exact potential at p and time t, gauged; zone as referencePhi
    scalar exactPhi
    (
        const point& p,
        const scalar t,
        const label zone = -1
    ) const;

    //- Converged direct forward solve on heart Vm into phi, sampled in the
    //  electrode cells
    void solveDirect
    (
        const volScalarField& Vm,
        volScalarField& phi,
        scalarField& values,
        const volScalarField* extraSource = nullptr
    );

    //- Error of the directStatic field against the exact potential at t
    void recordFieldError(const scalar t);

    void writeSummary();


protected:

    //- True once the channels and the reference are built (first sample)
    bool initialised() const
    {
        return leadFieldPtr_.valid();
    }

    //- The channels and conductor of the ECG domain (once initialised)
    const leadFieldECG& leadField() const
    {
        return *leadFieldPtr_;
    }

    //- Base-mesh cells of the heart, of the interface band and of the
    //  surface
    boolList heartCells_;
    boolList interfaceBand_;
    boolList surfaceCells_;


    // Reference hooks

        //- Name of the reference in the summary and the log
        virtual word referenceName() const = 0;

        //- Check the conductor against the reference and set its
        //  parameters (heartCells_ is set)
        virtual void initialiseReference(const leadFieldConductor& lfc) = 0;

        //- Exact heart Vm at p and time t
        virtual scalar exactVm(const point& p, const scalar t) const = 0;

        //- Exact potential at p and time t, not gauged; zone is the cell
        //  zone of cellZone(), or -1 for the zone containing p
        virtual scalar referencePhi
        (
            const point& p,
            const scalar t,
            const label zone
        ) const = 0;

        //- Zone of base cell cellI passed to referencePhi (-1: none)
        virtual label cellZone(const label) const
        {
            return -1;
        }

        //- Extra source of both static channels at t, or null: its
        //  lead-field term is added to leadFieldStatic, the direct solve
        //  adds it to b
        virtual tmp<volScalarField> staticSource
        (
            const volScalarField&,
            const scalar,
            scalarField&
        ) const
        {
            return nullptr;
        }

        //- Reference parameters appended to the log line
        virtual void reportParameters(Ostream&) const
        {}

        //- Reference parameters in the summary header, one "# key value"
        //  line each
        virtual void writeParameters(Ostream&) const
        {}


public:

    leadFieldECGVerifierBase
    (
        const electroStateProvider& stateProvider,
        const word& domainName,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    leadFieldECGVerifierBase(const leadFieldECGVerifierBase&) = delete;
    void operator=(const leadFieldECGVerifierBase&) = delete;

    virtual ~leadFieldECGVerifierBase() = default;

    virtual Requirements requirements() const override;

    virtual void record(const List<scalar>& numericValues) override;

    virtual void end() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
