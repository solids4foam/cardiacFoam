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
    leadFieldAdjoint

Description
    The explicit part C of the lead-field conductor's operator A = M + C,
    stored by coloured probing of passiveConductor::assemble, and the checks
    that make the adjoint solve A^T z = e exact for the assembled A.

    C x is minus the source of the operator assembled on x (base unknown
    phiE with the lead-field patch types, heart companion filled from x, no
    source, no reference). The sparsity pattern is the stencil the conductor
    publishes (passiveConductor::explicitStencil):
        face o face^p o S,
    with S the face stencil when every gradient scheme of
    passiveConductor::gradientNames() on the base and heart meshes is
    leastSquares or Gauss linear, the point stencil otherwise, and p the
    number of jump-aware gradient passes (nInterfaceGradientCorrectors with
    conormalHarmonic, 0 otherwise).

    Checks, each fatal above its tolerance:
    - linearity:    C(r1 + 2 r2) = C r1 + 2 C r2 on fields of the cell
                    centres (1e-11);
    - reproduction: the stored C equals the assembled C on another field
                    (1e-11);
    - transpose:    <z, C psi> = <C^T z, psi> (1e-12);
    - symmetry:     x^T M y = y^T M x, with the boundary diagonal and the
                    coupled-interface coefficients of M (1e-12);
    - residual:     e - M z - C^T z after each lead field (100 times the
                    leadFieldFixedPoint tolerance);
    - pin:          z_r d_r = 1 at a pinned reference cell, d_r its diagonal
                    before the pin (as the residual).
    The first four are round-off checks of exactly linear operators; the
    last two measure the convergence of the fixed point. The fields of the
    cell centres are the same in serial and in parallel.

    Gradient caching of the gradientNames() fields is fatal: a cached
    gradient can belong to another field of the same name.

SourceFiles
    leadFieldAdjoint.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef leadFieldAdjoint_H
#define leadFieldAdjoint_H

#include "leadFieldConductor.H"
#include "probedOperator.H"
#include "fvMatrices.H"

namespace Foam
{

class leadFieldAdjoint
{
    const leadFieldConductor& lfc_;

    // Tolerance of the residual and pin checks
    const scalar residualTolerance_;

    // Probe fields: base phiE with the lead-field patch types and its heart
    // companion
    mutable autoPtr<volScalarField> probePtr_;
    mutable autoPtr<volScalarField> heartProbePtr_;

    autoPtr<probedOperator> CPtr_;

    // C x assembled on x
    tmp<scalarField> assembledC(const scalarField& x) const;

    // Field of the cell centres (seed 1, 2, 3, ...)
    tmp<scalarField> centreField(const label seed) const;

    void checkLinearity() const;
    void checkReproduction() const;
    void checkTranspose() const;


public:

    //- Probe C of lfc's conductor and check it; the residual and pin
    //  tolerances follow the fixed-point tolerance
    leadFieldAdjoint
    (
        const leadFieldConductor& lfc,
        const leadFieldConductor::fixedPointControls& fixedPoint
    );

    leadFieldAdjoint(const leadFieldAdjoint&) = delete;
    void operator=(const leadFieldAdjoint&) = delete;


    //- C^T z
    tmp<scalarField> applyTranspose(const scalarField& z) const
    {
        return CPtr_->applyTranspose(z);
    }

    //- M x with the boundary diagonal and the coupled interfaces of psi
    //  (absolute coefficients and |x| when abs)
    static tmp<scalarField> matrixProduct
    (
        const fvScalarMatrix& M,
        const volScalarField& psi,
        const scalarField& x,
        const bool abs = false
    );

    //- Fatal unless M is symmetric
    void checkSymmetry
    (
        const fvScalarMatrix& M,
        const volScalarField& psi
    ) const;

    //- Fatal unless z solves M z + C^T z = e, and, with a pinned cell
    //  (refCell >= 0 on its processor, pinned true everywhere), unless
    //  z_r d_r = 1
    void checkSolution
    (
        const fvScalarMatrix& M,
        const volScalarField& z,
        const scalarField& e,
        const bool pinned,
        const label refCell,
        const scalar refDiag
    ) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
