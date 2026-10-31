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

const scalar linearityTolerance = 1e-11;
const scalar reproductionTolerance = 1e-11;
const scalar transposeTolerance = 1e-12;
const scalar symmetryTolerance = 1e-12;
const scalar residualTolerance = 1e-8;
const scalar pinTolerance = 1e-8;

// Gradient entries whose schemes enter C
const wordList gradientNames({"grad(phiE)", "grad(phiE.component(0))"});


// The resolved gradient schemes of the base and heart meshes, as text
std::string gradSchemesText(const fvMesh& base, const fvMesh& heart)
{
    OStringStream os;
    for (const fvMesh* meshPtr : {&base, &heart})
    {
        os  << (meshPtr == &base ? "base" : "; heart");
        for (const word& name : gradientNames)
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


// True for a gradient whose stencil is the cell and its face neighbours
bool faceStencilGradient(const fvMesh& mesh, const word& name)
{
    ITstream is(mesh.gradScheme(name));
    const word scheme(is);
    if (scheme == "leastSquares")
    {
        return true;
    }
    if (scheme == "Gauss" && !is.eof())
    {
        const word interpolation(is);
        return interpolation == "linear";
    }
    return false;
}


void report(const word& check, const scalar value, const scalar tolerance)
{
    Info<< "leadFieldAdjoint: " << check << ' ' << value
        << " (tolerance " << tolerance << ')' << endl;
}

} // End anonymous namespace


leadFieldAdjoint::leadFieldAdjoint(const leadFieldConductor& lfc)
:
    lfc_(lfc),
    probePtr_(lfc.leadFieldPotential().ptr()),
    heartProbePtr_(lfc.heartPotential().ptr()),
    CPtr_()
{
    for (const fvMesh* meshPtr : {&lfc.baseMesh(), &lfc.heartMesh()})
    {
        for (const word& name : gradientNames)
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

    const List<probedOperator::hop> hops(patternHops());
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


List<probedOperator::hop> leadFieldAdjoint::patternHops() const
{
    bool faceStencil = true;
    for (const fvMesh* meshPtr : {&lfc_.baseMesh(), &lfc_.heartMesh()})
    {
        for (const word& name : gradientNames)
        {
            faceStencil = faceStencil && faceStencilGradient(*meshPtr, name);
        }
    }

    // Divergence and face interpolation, interface-gradient passes, and the
    // gradient stencil
    List<probedOperator::hop> hops(1, probedOperator::hop::face);
    if
    (
        lfc_.conductor().interfaceConductivityInterpolation()
     == "conormalHarmonic"
    )
    {
        hops.append(probedOperator::hop::face);
        hops.append(probedOperator::hop::face);
    }
    hops.append
    (
        faceStencil ? probedOperator::hop::face : probedOperator::hop::point
    );
    return hops;
}


tmp<scalarField> leadFieldAdjoint::centreField(const label seed) const
{
    const fvMesh& mesh = lfc_.baseMesh();
    const vectorField& C = mesh.C().primitiveField();
    const vector origin(gMin(C));
    const scalar h = cbrt(gMin(mesh.V().field()));
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
            << lfc_.conductor().interfaceConductivityInterpolation() << '.'
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
    report("residual", value, residualTolerance);

    if (!(value <= residualTolerance))
    {
        FatalErrorInFunction
            << "Lead field does not solve A^T z = e: relative residual "
            << value << " (tolerance " << residualTolerance << ")."
            << exit(FatalError);
    }

    if (pinned)
    {
        const scalar zd = returnReduce
        (
            refCell >= 0 ? zI[refCell]*refDiag : 0, sumOp<scalar>()
        );
        const scalar pinValue = mag(zd - 1);
        report("pin", pinValue, pinTolerance);

        if (!(pinValue <= pinTolerance))
        {
            FatalErrorInFunction
                << "Lead field at the reference cell times its diagonal is "
                << zd << ", not 1 (tolerance " << pinTolerance << ")."
                << exit(FatalError);
        }
    }
}

} // End namespace Foam

// ************************************************************************* //
