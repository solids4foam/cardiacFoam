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
    manufacturedEikonalECGVerifier

Description
    Manufactured verifier for eikonal template-surrogate ECG.

    domain cube (default): Vm = U(t - psi) with psi = exp(k.x) on the unit
    cube and the uniform conductivity of the first cell, tensor-product
    Gauss-Legendre.

    domain shell: Vm = U(t - psi) with psi, grad psi and the conductivity
    G(x) of manufacturedShellEikonalReference.H on r1 < |x| < r2, spherical
    tensor rule of order q: q Gauss-Legendre nodes in r (Jacobian r^2), q in
    mu = cos(theta), 2q trapezoid nodes in phi. psi, G grad psi and the
    kernel are evaluated once per order; each sample is then the sum over
    the nodes of the kernel weight times U'(t - psi). Electrodes must lie
    outside the outer sphere.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef manufacturedEikonalECGVerifier_H
#define manufacturedEikonalECGVerifier_H

#include "ecgVerificationModel.H"
#include "OFstream.H"
#include "Switch.H"
#include "overrideTypeName.H"

namespace Foam
{

class manufacturedEikonalECGVerifier
:
    public ecgVerificationModel
{
    autoPtr<OFstream> outputPtr_;

    Switch enabled_;
    label dimension_;
    word domain_;
    label referenceQuadratureOrder_;
    List<label> checkQuadratureOrders_;
    vector k_;
    bool sineTemplate_;
    scalar templateDuration_;
    scalar templateRise_;

    label sampleCount_;
    List<scalar> referenceErrorL1Sum_;
    List<scalar> referenceErrorL2Sum_;
    List<scalar> referenceErrorLinf_;
    List<List<scalar>> checkErrorL1Sum_;
    List<List<scalar>> checkErrorL2Sum_;
    List<List<scalar>> checkErrorLinf_;
    List<List<scalar>> referenceDeltaL1Sum_;
    List<List<scalar>> referenceDeltaL2Sum_;
    List<List<scalar>> referenceDeltaLinf_;

    bool summaryWritten_;

    // Shell domain: per order (the checks, then the reference) the exact psi
    // at this rank's quadrature nodes and, per electrode, the node weight
    // times -(G grad psi).(x - e)/|x - e|^3
    bool shellCacheValid_;
    List<scalarField> shellPsi_;
    List<List<scalarField>> shellKernel_;

    void resizeStorage();
    void initialiseOutput();
    void updateStatistics
    (
        const List<scalar>& numericValues,
        const List<List<scalar>>& checkReferenceValues,
        const List<scalar>& referenceValues
    );
    void writeSummary();

    //- Shell domain: electrodes outside the outer sphere
    void checkShellElectrodes() const;

    //- Shell domain: chi and Cm of the header, conductivity G at every cell
    //  centre
    void checkShellConfiguration() const;

    void rebuildShellCache();

    //- Shell domain: the check references and the reference at sampleTime
    void shellReferences
    (
        const scalar sampleTime,
        List<List<scalar>>& checkReferenceValues,
        List<scalar>& referenceValues
    ) const;

    const dictionary& manufacturedDict(const dictionary& dict) const;

    manufacturedEikonalECGVerifier
    (
        const manufacturedEikonalECGVerifier&
    ) = delete;
    void operator=(const manufacturedEikonalECGVerifier&) = delete;

public:

    OverrideTypeName("manufacturedEikonalECGVerifier");

    manufacturedEikonalECGVerifier
    (
        const electroStateProvider& stateProvider,
        const word& domainName,
        const dictionary& dict,
        const wordList& electrodeNames,
        const List<vector>& electrodePositions
    );

    virtual ~manufacturedEikonalECGVerifier() = default;

    virtual Requirements requirements() const override;

    virtual bool read(const dictionary& dict) override;

    // Base class override for implicitly timed workflows
    virtual void record(const List<scalar>& numericValues) override;

    // Overload for eikonal ECG, which is driven by a post-processing pass
    // over activation times rather than a time-stepping loop, and thus
    // requires explicit time stamps.

    virtual void record
    (
        scalar sampleTime,
        const List<scalar>& numericValues
    ) override;

    virtual void end() override;
};

} // End namespace Foam

#endif

// ************************************************************************* //
