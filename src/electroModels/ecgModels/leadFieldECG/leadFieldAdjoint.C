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

\*---------------------------------------------------------------------------*/

#include "leadFieldAdjoint.H"
#include "fvMatrix.H"
#include "PstreamReduceOps.H"
#include <cmath>

namespace Foam
{

namespace
{

// The explicit part is linear and the probing exact, so the linearity,
// reproduction and transpose defects are round-off of the assembly and of
// the probed sums; a scheme that is not linear or whose stencil exceeds the
// probing pattern gives defects far above round-off
const scalar linearityTolerance = 1e-11;
const scalar reproductionTolerance = 1e-11;
const scalar transposeTolerance = 1e-12;

// The implicit matrix is assembled symmetric, so its defect is round-off of
// the products
const scalar symmetryTolerance = 1e-12;

// The relative residual of A^T z = e and the pin defect of a fixed point
// stopped at a relative change t between passes are of order t; the factor
// leaves room for the linear solver tolerance and round-off
const scalar fixedPointResidualFactor = 100;


// The resolved gradient schemes of the base and heart meshes, as text
std::string gradSchemesText(const fvMesh& base, const fvMesh& heart)
{
    OStringStream os;
    for (const fvMesh* meshPtr : {&base, &heart})
    {
        os  << (meshPtr == &base ? "base" : "; heart");
        for (const word& name : passiveConductor::gradientNames())
        {
            os  << ' ' << name << ':';
            for (const token& t : ITstream(meshPtr->gradScheme(name)))
            {
                os  << ' ' << t;
            }
        }
    }
    return os.str();
}


void report(const word& check, const scalar value, const scalar tolerance)
{
    Info<< "leadFieldAdjoint: " << check << ' ' << value
        << " (tolerance " << tolerance << ')' << endl;
}

} // End anonymous namespace


leadFieldAdjoint::leadFieldAdjoint
(
    const leadFieldConductor& lfc,
    const leadFieldConductor::fixedPointControls& fixedPoint
)
:
    lfc_(lfc),
    residualTolerance_(fixedPointResidualFactor*fixedPoint.tolerance),
    probePtr_(lfc.leadFieldPotential().ptr()),
    heartProbePtr_(lfc.heartPotential().ptr()),
    CPtr_()
{
    for (const fvMesh* meshPtr : {&lfc.baseMesh(), &lfc.heartMesh()})
    {
        for (const word& name : passiveConductor::gradientNames())
        {
            if (meshPtr->solution().cache(name))
            {
                FatalErrorInFunction
                    << "Lead fields need " << name << " recomputed for "
                    << "every field: remove it from the cache of "
                    << "fvSolution of mesh '" << meshPtr->name() << "'."
                    << exit(FatalError);
            }
        }
    }

    checkLinearity();

    const List<probedOperator::hop> hops(lfc.conductor().explicitStencil());
    CPtr_.reset
    (
        new probedOperator
        (
            lfc.baseMesh(),
            hops,
            [this](const scalarField& x) { return assembledC(x); }
        )
    );
    Info<< "leadFieldAdjoint: " << hops.size() << " pattern hops, "
        << CPtr_->nColours() << " probes, " << CPtr_->nNonZeros()
        << " entries" << endl;

    checkReproduction();
    checkTranspose();
}


tmp<scalarField> leadFieldAdjoint::assembledC(const scalarField& x) const
{
    volScalarField& probe = probePtr_.ref();
    volScalarField& heartProbe = heartProbePtr_.ref();
    const passiveConductor& conductor = lfc_.conductor();

    probe.primitiveFieldRef() = x;
    probe.correctBoundaryConditions();
    conductor.fillHeart(probe, heartProbe);

    tmp<fvScalarMatrix> tEqn(conductor.assemble(probe, heartProbe));
    return tmp<scalarField>(new scalarField(-tEqn().source()));
}


tmp<scalarField> leadFieldAdjoint::centreField(const label seed) const
{
    const fvMesh& mesh = lfc_.baseMesh();
    const vectorField& C = mesh.C().primitiveField();
    const vector origin(gMin(C));
    const scalar h = cbrt(gMin(mesh.V().field()));

    // r = fract(43758.5453 sin(x.a + 1.618 seed)) - 1/2, x the cell centre
    // in units of h from the global origin: a value in [-1/2, 1/2) without
    // a pattern from cell to cell, the same for a cell in serial and in
    // parallel. a and 43758.5453 are arbitrary incommensurate hash
    // constants (those of the common sin-hash one-liner, a third component
    // added for 3D); 1.618 shifts the phase so that each seed gives
    // another field
    const vector a(12.9898, 78.233, 37.719);

    tmp<scalarField> tr(new scalarField(C.size()));
    scalarField& r = tr.ref();
    forAll(C, cellI)
    {
        const scalar s =
            Foam::sin(((C[cellI] - origin)/h & a) + 1.618*seed)*43758.5453;
        r[cellI] = s - std::floor(s) - 0.5;
    }
    return tr;
}


void leadFieldAdjoint::checkLinearity() const
{
    const scalarField r1(centreField(1));
    const scalarField r2(centreField(2));
    const scalarField C1(assembledC(r1));
    const scalarField C2(assembledC(r2));
    const scalarField C12(assembledC(r1 + 2*r2));

    const scalar scale = gMax(mag(C1) + 2*mag(C2));
    const scalar value = gMax(mag(C12 - C1 - 2*C2))/max(scale, VSMALL);
    report("linearity", value, linearityTolerance);

    if (!(value <= linearityTolerance))
    {
        FatalErrorInFunction
            << "The explicit part of the lead-field operator is not linear "
            << "(relative defect " << value << ", tolerance "
            << linearityTolerance << "). Lead fields need linear schemes: "
            << "laplacian snGrad corrected, uncorrected or orthogonal, a "
            << "linear interpolation of the cross term and unlimited "
            << "gradients. Gradient schemes: "
            << gradSchemesText(lfc_.baseMesh(), lfc_.heartMesh()).c_str()
            << "; see also laplacianSchemes and interpolationSchemes of "
            << "fvSchemes."
            << exit(FatalError);
    }
}


void leadFieldAdjoint::checkReproduction() const
{
    const scalarField r(centreField(3));
    const scalarField Cr(assembledC(r));
    const scalarField Cpr(CPtr_->apply(r));

    const scalar scale = gMax(CPtr_->applyAbs(r));
    const scalar value = gMax(mag(Cpr - Cr))/max(scale, VSMALL);
    report("reproduction", value, reproductionTolerance);

    if (!(value <= reproductionTolerance))
    {
        FatalErrorInFunction
            << "The probed explicit operator does not reproduce the "
            << "assembled one (relative difference " << value
            << ", tolerance " << reproductionTolerance << "): a scheme is "
            << "nonlinear or has a stencil wider than the probing pattern. "
            << "Gradient schemes: "
            << gradSchemesText(lfc_.baseMesh(), lfc_.heartMesh()).c_str()
            << "; interface option "
            << passiveConductor::interfaceInterpolationNames
               [
                   lfc_.conductor().torso().interface
               ]
            << ", nInterfaceGradientCorrectors "
            << lfc_.conductor().nInterfaceGradientCorrectors() << '.'
            << exit(FatalError);
    }
}


void leadFieldAdjoint::checkTranspose() const
{
    const scalarField z(centreField(4));
    const scalarField psi(centreField(5));

    const scalar zCpsi = gSum(z*assembledC(psi));
    const scalar Ctzpsi = gSum(CPtr_->applyTranspose(z)*psi);
    const scalar scale = gSum(mag(z)*CPtr_->applyAbs(psi));
    const scalar value = mag(zCpsi - Ctzpsi)/max(scale, VSMALL);
    report("transpose", value, transposeTolerance);

    if (!(value <= transposeTolerance))
    {
        FatalErrorInFunction
            << "<z, C psi> = " << zCpsi << " and <C^T z, psi> = " << Ctzpsi
            << " differ by " << value << " relative (tolerance "
            << transposeTolerance << ")."
            << exit(FatalError);
    }
}


tmp<scalarField> leadFieldAdjoint::matrixProduct
(
    const fvScalarMatrix& M,
    const volScalarField& psi,
    const scalarField& x,
    const bool abs
)
{
    // The diagonal the linear solver sees: with every patch's internal
    // coefficients, coupled ones included
    lduMatrix A(M);
    scalarField& diag = A.diag();
    forAll(M.internalCoeffs(), patchI)
    {
        const labelUList& cells = M.lduAddr().patchAddr(patchI);
        const scalarField& coeffs = M.internalCoeffs()[patchI];
        forAll(cells, faceI)
        {
            diag[cells[faceI]] += coeffs[faceI];
        }
    }

    // Amul subtracts the interface coefficient times the neighbour value
    FieldField<Field, scalar> interfaceCoeffs(M.boundaryCoeffs());
    if (abs)
    {
        A.diag() = mag(A.diag());
        A.upper() = mag(A.upper());
        forAll(interfaceCoeffs, patchI)
        {
            interfaceCoeffs[patchI] = -mag(interfaceCoeffs[patchI]);
        }
    }

    tmp<scalarField> tAx(new scalarField(x.size()));
    A.Amul
    (
        tAx.ref(),
        tmp<scalarField>(abs ? scalarField(mag(x)) : x),
        interfaceCoeffs,
        psi.boundaryField().scalarInterfaces(),
        0
    );
    return tAx;
}


void leadFieldAdjoint::checkSymmetry
(
    const fvScalarMatrix& M,
    const volScalarField& psi
) const
{
    const scalarField x(centreField(6));
    const scalarField y(centreField(7));

    const scalar xMy = gSum(x*matrixProduct(M, psi, y));
    const scalar yMx = gSum(y*matrixProduct(M, psi, x));
    const scalar scale = gSum(mag(x)*matrixProduct(M, psi, y, true));
    const scalar value = mag(xMy - yMx)/max(scale, VSMALL);
    report("symmetry", value, symmetryTolerance);

    if (!(value <= symmetryTolerance))
    {
        FatalErrorInFunction
            << "The implicit lead-field matrix is not symmetric: x^T M y = "
            << xMy << ", y^T M x = " << yMx << ", relative " << value
            << " (tolerance " << symmetryTolerance << ")."
            << exit(FatalError);
    }
}


void leadFieldAdjoint::checkSolution
(
    const fvScalarMatrix& M,
    const volScalarField& z,
    const scalarField& e,
    const bool pinned,
    const label refCell,
    const scalar refDiag
) const
{
    const scalarField& zI = z.primitiveField();
    const scalarField r
    (
        e - matrixProduct(M, z, zI) - CPtr_->applyTranspose(zI)
    );
    const scalar scale = gMax
    (
        matrixProduct(M, z, zI, true) + CPtr_->applyTransposeAbs(zI)
    );
    const scalar value = gMax(mag(r))/max(scale, VSMALL);
    report("residual", value, residualTolerance_);

    if (!(value <= residualTolerance_))
    {
        FatalErrorInFunction
            << "Lead field does not solve A^T z = e: relative residual "
            << value << " (tolerance " << residualTolerance_ << ")."
            << exit(FatalError);
    }

    if (pinned)
    {
        const scalar zd = returnReduce
        (
            refCell >= 0 ? zI[refCell]*refDiag : 0, sumOp<scalar>()
        );
        const scalar pinValue = mag(zd - 1);
        report("pin", pinValue, residualTolerance_);

        if (!(pinValue <= residualTolerance_))
        {
            FatalErrorInFunction
                << "Lead field at the reference cell times its diagonal is "
                << zd << ", not 1 (tolerance " << residualTolerance_ << ")."
                << exit(FatalError);
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
