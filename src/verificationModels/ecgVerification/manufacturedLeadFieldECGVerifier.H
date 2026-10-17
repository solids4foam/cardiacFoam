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
    manufacturedLeadFieldECGVerifier

Description
    Lead-field ECG verifier. At each sample and electrode it compares:
    - leadField: the channel of the leadFieldECG domain (z_e^T b(Vm_h) + c_e);
    - direct: a direct forward solve of the same conductor on Vm_h, with the
      same boundary values and reference, assembly and solve repeated until
      the relative change between passes is below 1e-10, sampled in the
      electrode cell;
    - exactPoint, exactCentroid: the exact potential of the selected
      reference at the electrode and at its cell centroid, for the exact Vm;
    - leadFieldStatic, directStatic: the lead field and the direct solve
      applied to the exact Vm at the heart cell centres.

    With a pinned reference the exact potential is shifted to equal the
    reference value at the centroid of the pinned cell.

    Writes postProcessing/manufacturedLeadFieldECG_<domain>.dat (samples)
    and manufacturedLeadFieldECGSummary_<domain>.dat (rows name L1 L2 Linf
    over samples, per electrode and as the per-sample maximum over the
    electrodes).

    Keys of verificationModel:
        reference   slab;    // required: the exact solution, slab or layered
        layerInterfaces (-0.5 1.5);  // layered only, required: xL < 0, xR > 1

    layered takes alpha from the surface-current patches of the conductor
    (faces normal to x, current times n_x equal to alpha on every face, no
    net current; none for alpha = 0) and the layer conductivities from its
    torso cells, which must each lie in one layer with one isotropic value
    per layer.

    The direct solve uses the fvSolution entry phiE.

SourceFiles
    manufacturedLeadFieldECGVerifier.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedLeadFieldECGVerifier_H
#define manufacturedLeadFieldECGVerifier_H

#include "ecgVerificationModel.H"
#include "leadFieldECG.H"
#include "ecgVerification/manufacturedLeadFieldECGReference.H"
#include "OFstream.H"
#include "Enum.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedLeadFieldECGVerifier
:
    public ecgVerificationModel
{
public:

    //- Exact solutions
    enum class referenceType
    {
        slab,
        layered
    };

    static const Enum<referenceType> referenceTypeNames;


private:

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

    referenceType reference_;

    //- The channels and conductor of the ECG domain, evaluated on any Vm
    autoPtr<leadFieldECG> leadFieldPtr_;

    //- Slab and layered: Gi_xx/(Gi_xx + Ge_xx)
    scalar slabK_;

    //- Layered: the layers and the current density -alpha along x
    manufacturedLeadFieldLayers layers_;
    scalar alpha_;

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

    bool summaryWritten_;

    static const wordList quantityNames_;

    //- Build the channels and the reference at the first sample
    void initialise();

    //- Layered: alpha from the surface currents, the layer conductivities
    //  from the torso cells
    void initialiseLayered(const leadFieldConductor& lfc);

    //- Exact potential at p and time t, gauged
    scalar exactPhi(const point& p, const scalar t) const;

    //- Exact heart Vm at p and time t
    scalar exactVm(const point& p, const scalar t) const;

    //- Converged direct forward solve on heart Vm into phi, sampled in the
    //  electrode cells
    void solveDirect
    (
        const volScalarField& Vm,
        volScalarField& phi,
        scalarField& values
    );

    void writeSummary();

    manufacturedLeadFieldECGVerifier
    (
        const manufacturedLeadFieldECGVerifier&
    ) = delete;
    void operator=(const manufacturedLeadFieldECGVerifier&) = delete;


public:

    OverrideTypeName("manufacturedLeadFieldECGVerifier");

    manufacturedLeadFieldECGVerifier
    (
        const electroStateProvider& stateProvider,
        const word& domainName,
        const dictionary& dict,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual ~manufacturedLeadFieldECGVerifier() = default;

    virtual Requirements requirements() const override;

    virtual bool read(const dictionary& dict) override;

    virtual void record(const List<scalar>& numericValues) override;

    virtual void end() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
