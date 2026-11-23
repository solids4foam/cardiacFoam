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

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "boundedCorrectionGaussLaplacianScheme.H"
#include "centredCPCCellToCellStencilObject.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class GType>
const Foam::Enum
<
    typename Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::limiterType
>
Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::limiterTypeNames
({
    { limiterType::MC, "MC" },
    { limiterType::vanLeer, "vanLeer" },
    { limiterType::minmod, "minmod" },
    { limiterType::linear, "linear" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GType>
inline Foam::scalar Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::limit
(
    const scalar a,
    const scalar b,
    const scalar w
) const
{
    if (limiter_ == limiterType::linear)
    {
        return w*a + (1 - w)*b;
    }
    if (a*b <= 0)
    {
        return 0;
    }
    switch (limiter_)
    {
        case limiterType::MC:
        {
            return sign(a)
               *min(2*min(mag(a), mag(b)), mag(w*a + (1 - w)*b));
        }
        case limiterType::vanLeer:
        {
            return 2*a*b/(a + b);
        }
        case limiterType::minmod:
        {
            return sign(a)*min(mag(a), mag(b));
        }
        default:
        {
            return 0;
        }
    }
}


namespace Foam
{
namespace fv
{

template<>
tmp<surfaceVectorField> boundedCorrectionGaussLaplacianScheme<scalar>::SfGamma
(
    const surfaceScalarField& gamma
)
{
    return gamma.mesh().Sf()*gamma;
}

}
}


template<class GType>
Foam::tmp<Foam::surfaceVectorField>
Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::SfGamma
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma
)
{
    return gamma.mesh().Sf() & gamma;
}


template<class GType>
Foam::tmp<Foam::surfaceVectorField>
Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::correctionVector
(
    const surfaceVectorField& SfG,
    const surfaceScalarField& SfGammaSn
) const
{
    const fvMesh& mesh = this->mesh();

    // Tensor part: the component of S & gamma in the face plane
    tmp<surfaceVectorField> tc
    (
        new surfaceVectorField
        (
            "boundedCorrectionVector",
            SfG - SfGammaSn*mesh.Sf()/mesh.magSf()
        )
    );

    // Non-orthogonal part, as correctedSnGrad: k & grad, k = n - delta/(n & delta)
    if (this->tsnGradScheme_().corrected())
    {
        tc.ref() += SfGammaSn*mesh.nonOrthCorrectionVectors();
    }

    return tc;
}


template<class GType>
Foam::tmp<Foam::surfaceScalarField>
Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::limitedCorrection
(
    const surfaceVectorField& c,
    const volScalarField& vf
) const
{
    const fvMesh& mesh = this->mesh();

    volVectorField g
    (
        "boundedCorrectionGrad(" + vf.name() + ')',
        fvc::grad(vf)
    );
    vectorField& gi = g.primitiveFieldRef();
    label nExtrema = 0;

    if (limiter_ != limiterType::linear)
    {
        // Values of the cell-point-cell stencil, boundary faces included
        const extendedCentredCellToCellStencil& stencil =
            centredCPCCellToCellStencilObject::New(mesh);

        List<scalar> flatVf(stencil.map().constructSize(), Zero);
        forAll(vf, celli)
        {
            flatVf[celli] = vf[celli];
        }
        forAll(vf.boundaryField(), patchi)
        {
            const fvPatchScalarField& psf = vf.boundaryField()[patchi];
            label nCompact =
                psf.patch().start() - mesh.nInternalFaces() + mesh.nCells();
            forAll(psf, i)
            {
                flatVf[nCompact++] = psf[i];
            }
        }
        stencil.map().distribute(flatVf);

        // Gradient scaled by phi = min(1, r/band), r = distance of psi from
        // the nearer of the stencil minimum and maximum over the stencil
        // range: zero at a local extremum, continuous in psi for band > 0;
        // band = 0 zeroes it at extrema only
        const labelListList& stencilAddr = stencil.stencil();
        forAll(gi, celli)
        {
            const scalar v = vf[celli];
            scalar vMin = v;
            scalar vMax = v;
            for (const label j : stencilAddr[celli])
            {
                vMin = min(vMin, flatVf[j]);
                vMax = max(vMax, flatVf[j]);
            }
            const scalar r = min(v - vMin, vMax - v);
            scalar phi = 1;
            if (r <= 0)
            {
                phi = 0;
                ++nExtrema;
            }
            else if (r < band_*(vMax - vMin))
            {
                phi = r/(band_*(vMax - vMin));
            }
            gi[celli] *= phi;
        }
        g.correctBoundaryConditions();
    }

    tmp<surfaceScalarField> tcorr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "boundedCorrection(" + vf.name() + ')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<scalar>
            (
                c.dimensions()*vf.dimensions()/dimLength, Zero
            )
        )
    );
    surfaceScalarField& corr = tcorr.ref();
    corr.setOriented();

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& w = mesh.weights();
    label nLimited = 0;

    forAll(own, facei)
    {
        const scalar a = c[facei] & gi[own[facei]];
        const scalar b = c[facei] & gi[nei[facei]];
        corr[facei] = limit(a, b, w[facei]);
        if (debug && corr[facei] != w[facei]*a + (1 - w[facei])*b)
        {
            ++nLimited;
        }
    }

    surfaceScalarField::Boundary& corrBf = corr.boundaryFieldRef();
    forAll(corrBf, patchi)
    {
        const fvsPatchVectorField& pc = c.boundaryField()[patchi];
        const fvPatchVectorField& pg = g.boundaryField()[patchi];
        fvsPatchScalarField& pcorr = corrBf[patchi];
        const vectorField gP(pg.patchInternalField());

        if (pg.coupled())
        {
            const vectorField gN(pg.patchNeighbourField());
            const scalarField& pw = mesh.weights().boundaryField()[patchi];
            forAll(pcorr, i)
            {
                pcorr[i] = limit(pc[i] & gP[i], pc[i] & gN[i], pw[i]);
            }
        }
        else
        {
            forAll(pcorr, i)
            {
                pcorr[i] = pc[i] & gP[i];
            }
        }
    }

    if (debug)
    {
        Info<< typeName << '(' << vf.name() << "): "
            << returnReduce(nExtrema, sumOp<label>()) << " extremum cells of "
            << returnReduce(mesh.nCells(), sumOp<label>()) << ", "
            << returnReduce(nLimited, sumOp<label>())
            << " limited internal faces of "
            << returnReduce(mesh.nInternalFaces(), sumOp<label>()) << endl;
    }

    return tcorr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class GType>
Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::
boundedCorrectionGaussLaplacianScheme
(
    const fvMesh& mesh,
    Istream& is
)
:
    gaussLaplacianScheme<scalar, GType>(mesh, is),
    limiter_(limiterTypeNames.read(is)),
    band_(0)
{
    const ITstream* itsPtr = dynamic_cast<const ITstream*>(&is);
    if (itsPtr && itsPtr->nRemainingTokens())
    {
        band_ = readScalar(is);
    }
    if (band_ < 0 || band_ > 0.5)
    {
        FatalIOErrorInFunction(is)
            << typeName << ": the extremum band must lie in [0, 0.5], not "
            << band_ << exit(FatalIOError);
    }

    const word snGradType(this->tsnGradScheme_().type());
    if (snGradType != "corrected" && snGradType != "uncorrected")
    {
        FatalIOErrorInFunction(is)
            << typeName << " needs the corrected or uncorrected snGrad scheme,"
            << " not " << snGradType << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GType>
Foam::tmp<Foam::fvScalarMatrix>
Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::fvmLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const volScalarField& vf
)
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField SfG(SfGamma(gamma));
    const surfaceScalarField SfGammaSn(SfG & (mesh.Sf()/mesh.magSf()));

    tmp<fvScalarMatrix> tfvm = this->fvmLaplacianUncorrected
    (
        SfGammaSn,
        this->tsnGradScheme_().deltaCoeffs(vf),
        vf
    );
    fvScalarMatrix& fvm = tfvm.ref();

    tmp<surfaceScalarField> tcorr
    (
        limitedCorrection(correctionVector(SfG, SfGammaSn), vf)
    );

    fvm.source() -= mesh.V()*fvc::div(tcorr())().primitiveField();

    if (mesh.fluxRequired(vf.name()))
    {
        fvm.faceFluxCorrectionPtr(tcorr.ptr());
    }

    return tfvm;
}


template<class GType>
Foam::tmp<Foam::volScalarField>
Foam::fv::boundedCorrectionGaussLaplacianScheme<GType>::fvcLaplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const volScalarField& vf
)
{
    const fvMesh& mesh = this->mesh();

    const surfaceVectorField SfG(SfGamma(gamma));
    const surfaceScalarField SfGammaSn(SfG & (mesh.Sf()/mesh.magSf()));

    tmp<volScalarField> tLaplacian
    (
        fvc::div
        (
            SfGammaSn
           *snGradScheme<scalar>::snGrad
            (
                vf, this->tsnGradScheme_().deltaCoeffs(vf), "snGrad"
            )
          + limitedCorrection(correctionVector(SfG, SfGammaSn), vf)
        )
    );

    tLaplacian.ref().rename
    (
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * Selection * * * * * * * * * * * * * * * * //

#define makeBoundedCorrectionGaussLaplacianScheme(GType)                       \
    typedef Foam::fv::boundedCorrectionGaussLaplacianScheme<Foam::GType>       \
        boundedCorrectionGaussLaplacianScheme##GType;                          \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        boundedCorrectionGaussLaplacianScheme##GType, 0                        \
    );                                                                         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fv                                                           \
        {                                                                      \
            laplacianScheme<scalar, GType>::addIstreamConstructorToTable       \
            <                                                                  \
                boundedCorrectionGaussLaplacianScheme<GType>                   \
            > addBoundedCorrectionGauss##GType##IstreamConstructorToTable_;    \
        }                                                                      \
    }

makeBoundedCorrectionGaussLaplacianScheme(scalar)
makeBoundedCorrectionGaussLaplacianScheme(symmTensor)
makeBoundedCorrectionGaussLaplacianScheme(tensor)


// ************************************************************************* //
