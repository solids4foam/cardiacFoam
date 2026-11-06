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
    coupled1D3DMonodomainVerifier

Description
    Verifier for the coupled 1D-3D manufactured monodomain problem. The
    tissue field is V3 = sqrt(1 + t) F3 (tissueField box: the FDA cube
    field; shell: the sealed-shell field) and the graph field is
    V1 = sqrt(1 + t) cos(pi (x - xOrigin)).

    reference discrete (the original self-consistency check): the junction
    forcing is built with a second pvjMapper, so map errors cancel and only
    the feedback of the discretisation error through the junction is seen.

    reference continuum: the junction forcing is the continuum source of
    pvjContinuumReference, normalised over the tissue part of the support
    and averaged over each cell; the graph forcing uses the continuum
    gathered voltage. The discrete map's error then reaches the solution.
    The forcing is evaluated at the new time level t0 + dt of each phase,
    so the staggered lag between graph and tissue is not reproduced.

    Both references write, under verification/, the map and local files
    below; continuum mode also writes the steps file:
      coupled1D3DMonodomain_map.csv      discrete map vs continuum map of
                                         every diagnostic kernel at t = 0
      coupled1D3DMonodomain_steps.csv    per step: junction charge leaving
                                         the graph (from its own budget) and
                                         entering the tissue (from the
                                         deposited source), and the
                                         terminal currents
      coupled1D3DMonodomain_local.csv    final-time errors local to each
                                         support and at each terminal

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef coupled1D3DMonodomainVerifier_H
#define coupled1D3DMonodomainVerifier_H

#include "couplingVerificationModel.H"
#include "overrideTypeName.H"
#include "scalarField.H"
#include "OFstream.H"
#include "fvMesh.H"
#include "pvjContinuumReference.H"

namespace Foam
{

class coupled1D3DMonodomainVerifier : public couplingVerificationModel
{
    word reference_;
    word tissueField_;
    scalar xOrigin_;
    wordList diagnosticKernels_;

    mutable bool diagnosticsWritten_;
    mutable scalarField exactPrimarySource_;
    mutable scalarField exactSecondaryAppliedCurrent_;

    // Continuum map of the configured kernel, per terminal

        bool continuumReady_;
        scalar R_;
        word kernel_;
        scalarField Wc_;
        scalarField Fc_;
        List<labelList> refCells_;
        List<scalarField> refOmega_;
        List<scalarField> refOmegaF_;

    // Per-step charge balance

        scalarField graphVmOld_;
        scalarField sourceBeforeDeposit_;
        autoPtr<OFstream> stepFile_;
        scalar maxChargeImbalance_;
        scalarField maxCurrentError_;

    bool continuum() const
    {
        return reference_ == "continuum";
    }

    //- Manufactured tissue field F3 (V3 = sqrt(1 + t) F3)
    scalar tissueF(const point& x) const;

    //- Manufactured graph field F1 (V1 = sqrt(1 + t) F1)
    scalar graphF(const point& x) const;

    pvjContinuumReference::support supportGeometry
    (
        const fvMesh& mesh,
        const point& t,
        scalar R
    ) const;

    //- Continuum map of the configured kernel; the cell integrals of w F
    //  only for the implicit scheme's local operator
    void initialiseContinuum
    (
        const fvMesh& mesh,
        const pointField& terminalLocations,
        bool implicitCoupling
    );

    void writeMapDiagnostics
    (
        const fvMesh& mesh,
        const pointField& terminalLocations
    ) const;

    //- Continuum junction current at time t
    scalarField continuumCurrent
    (
        const pointField& terminalLocations,
        const scalarField& R_pvj,
        scalar t
    ) const;

    void updateContinuumSource
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        scalar evaluationTime,
        bool implicitCoupling,
        bool bidirectionalCoupling,
        const word& phaseName
    );

    void postProcessContinuum
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain
    );

public:

    OverrideTypeName("coupled1D3DMonodomainVerifier");

    explicit coupled1D3DMonodomainVerifier(const dictionary& dict);

    virtual void preProcess
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain
    ) override;

    virtual void updateManufacturedSource
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        scalar primaryTime,
        scalar secondaryTime,
        bool implicitCoupling,
        bool bidirectionalCoupling,
        const word& phaseName
    ) override;

    virtual void recordPrimaryDeposit
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain,
        scalar t0,
        scalar dt
    ) override;

    virtual void postProcess
    (
        tissueCouplingEndpoint& primaryDomain,
        electroDomainInterface& secondaryDomain
    ) override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
