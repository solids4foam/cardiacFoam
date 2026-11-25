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
    Foam::fv::boundedCorrectionGaussLaplacianScheme

Description
    Gauss Laplacian of a scalar with a tensor (or scalar) face conductivity
    whose explicit correction is limited so that a local extremum of the
    field is not deepened by diffusion.

    The face flux of the stock scheme (Gauss linear corrected) is

        F_f = alpha_f (psi_N - psi_P) + c_f & grad(psi)_f,

        alpha_f = (S_f & gamma_f & n_f) deltaCoeffs_f,
        c_f     = S_f & gamma_f - (S_f & gamma_f & n_f) (n_f - k_f),

    with k_f the non-orthogonal correction vector (zero with an uncorrected
    snGrad) and grad(psi)_f the linear interpolate of the cell gradients.
    The first term is the implicit two-point flux; c_f & grad(psi)_f holds
    both explicit corrections (tensor and non-orthogonal).

    This scheme keeps the two-point flux and replaces the explicit correction
    with a limited value of the two one-sided values

        q_P = c_f & g_P,  q_N = c_f & g_N,

    where g is the cell gradient (gradSchemes entry grad(psi)) scaled by

        phi_P = min(1, r_P/band),
        r_P   = min(psi_P - min_S psi, max_S psi - psi_P)/(max_S psi - min_S psi)

    over the cell-point-cell stencil S of P (boundary faces included):
    zero at a local minimum or maximum, one away from it, continuous in psi
    for band > 0. band = 0 (the default) only zeroes it at extrema. The
    limiter is zero when q_P q_N <= 0:

        MC       sign(q_P) min(2|q_P|, 2|q_N|, |w q_P + (1-w) q_N|)
        vanLeer  2 q_P q_N/(q_P + q_N)
        minmod   sign(q_P) min(|q_P|, |q_N|)
        linear   w q_P + (1-w) q_N with no extremum test: the stock
                 correction, for checking

    with w the linear interpolation weight; MC returns the stock value where
    q_P/q_N lies between about 1/3 and 3. On a non-coupled boundary face
    the correction is c_f & g_P.

    At a local extremum every correction flux of the cell is zero and the
    two-point fluxes point towards it, so diffusion does not deepen it. With
    band = 0 the switch is discontinuous and a cell next to an extremum can
    still pass it in one step; band > 0 makes the correction go to zero
    continuously as the cell approaches the extremum
    (after Sharma and Hammett, J. Comput. Phys. 227 (2007) 123, who limit
    the transverse flux only). With implicit Euler and converged correctors
    the global minimum and maximum of psi are not exceeded by diffusion;
    with a lagged correction or a two-step time scheme this is not proved.

    Selected in fvSchemes:

    laplacianSchemes
    {
        default         boundedCorrectionGauss linear corrected MC 0.25;
    }

    The last number is the band, in [0, 0.5]; it may be left out (0).

    The snGrad scheme must be corrected or uncorrected.

SourceFiles
    boundedCorrectionGaussLaplacianScheme.C

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#ifndef boundedCorrectionGaussLaplacianScheme_H
#define boundedCorrectionGaussLaplacianScheme_H

#include "gaussLaplacianScheme.H"
#include "Enum.H"

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
            Class boundedCorrectionGaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class GType>
class boundedCorrectionGaussLaplacianScheme
:
    public gaussLaplacianScheme<scalar, GType>
{
public:

    //- Limiter of the two one-sided correction values
    enum class limiterType { MC, vanLeer, minmod, linear };

    static const Enum<limiterType> limiterTypeNames;


private:

    limiterType limiter_;

    //- Extremum band: the cell gradient is scaled by min(1, r/band)
    scalar band_;


    //- Limited value of the one-sided values a, b; w weights a
    inline scalar limit(const scalar a, const scalar b, const scalar w) const;

    //- Face vector S & gamma
    static tmp<surfaceVectorField> SfGamma
    (
        const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma
    );

    //- Correction vector c_f of the face flux
    tmp<surfaceVectorField> correctionVector
    (
        const surfaceVectorField& SfGamma,
        const surfaceScalarField& SfGammaSn
    ) const;

    //- Limited explicit correction flux c_f & grad(psi)
    tmp<surfaceScalarField> limitedCorrection
    (
        const surfaceVectorField& c,
        const volScalarField& vf
    ) const;

    boundedCorrectionGaussLaplacianScheme
    (
        const boundedCorrectionGaussLaplacianScheme&
    ) = delete;

    void operator=(const boundedCorrectionGaussLaplacianScheme&) = delete;


public:

    TypeName("boundedCorrectionGauss");


    // Constructors

        boundedCorrectionGaussLaplacianScheme(const fvMesh& mesh, Istream& is);


    virtual ~boundedCorrectionGaussLaplacianScheme() = default;


    // Member Functions

        virtual tmp<fvScalarMatrix> fvmLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const volScalarField&
        );

        virtual tmp<volScalarField> fvcLaplacian
        (
            const GeometricField<GType, fvsPatchField, surfaceMesh>&,
            const volScalarField&
        );
};


} // End namespace fv
} // End namespace Foam

#endif

// ************************************************************************* //
