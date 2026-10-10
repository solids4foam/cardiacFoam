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

#include "interfaceGradient.H"
#include "exposedPhiETrace.H"
#include "surfaceFields.H"
#include "fvc.H"

namespace Foam
{

tmp<volVectorField> interfaceGradient
(
    const volScalarField& phi,
    const volTensorField& sigmaE,
    const volTensorField& sigmaI,
    const volVectorField& gradPrev
)
{
    const fvMesh& m = phi.mesh();
    const vectorField& C = m.C().primitiveField();
    const surfaceVectorField& Cf = m.Cf();
    const surfaceVectorField& Sf = m.Sf();
    const surfaceScalarField& magSf = m.magSf();
    const labelUList& own = m.owner();
    const labelUList& nei = m.neighbour();
    const scalarField& w = m.weights().primitiveField();

    surfaceScalarField phiF(linearInterpolate(phi));
    scalarField& phiFI = phiF.primitiveFieldRef();

    forAll(own, faceI)
    {
        const label o = own[faceI];
        const label nb = nei[faceI];
        const bool oHeart = magSqr(sigmaI[o]) > SMALL;
        if (oHeart == (magSqr(sigmaI[nb]) > SMALL))
        {
            continue;
        }
        const label h = oHeart ? o : nb;
        const label b = oHeart ? nb : o;
        const vector n((oHeart ? 1 : -1)*Sf[faceI]/magSf[faceI]);
        phiFI[faceI] = interfaceValue
        (
            n, Cf[faceI], C[h], C[b],
            sigmaE[h],
            sigmaE[b],
            phi[h],
            phi[b],
            w[faceI]*gradPrev[o] + (1 - w[faceI])*gradPrev[nb]
        );
    }

    forAll(m.boundary(), patchI)
    {
        const fvPatch& p = m.boundary()[patchI];
        if (!p.coupled())
        {
            continue;
        }
        const scalarField phiN
        (
            phi.boundaryField()[patchI].patchNeighbourField()
        );
        const tensorField sigmaEN
        (
            sigmaE.boundaryField()[patchI].patchNeighbourField()
        );
        const tensorField sigmaIN
        (
            sigmaI.boundaryField()[patchI].patchNeighbourField()
        );
        const vectorField gradN
        (
            gradPrev.boundaryField()[patchI].patchNeighbourField()
        );
        const vectorField nf(p.nf());
        const vectorField delta(p.delta());
        const scalarField& pw = m.weights().boundaryField()[patchI];
        const labelUList& fc = p.faceCells();
        scalarField& phiFp = phiF.boundaryFieldRef()[patchI];
        forAll(phiFp, faceI)
        {
            const label c = fc[faceI];
            const bool cHeart = magSqr(sigmaI[c]) > SMALL;
            if (cHeart == (magSqr(sigmaIN[faceI]) > SMALL))
            {
                continue;
            }
            const point& CfP = Cf.boundaryField()[patchI][faceI];
            const point CN(C[c] + delta[faceI]);
            const vector gradF
            (
                pw[faceI]*gradPrev[c] + (1 - pw[faceI])*gradN[faceI]
            );
            phiFp[faceI] =
                cHeart
              ? interfaceValue
                (
                    nf[faceI], CfP, C[c], CN, sigmaE[c], sigmaEN[faceI],
                    phi[c], phiN[faceI], gradF
                )
              : interfaceValue
                (
                    -nf[faceI], CfP, CN, C[c], sigmaEN[faceI], sigmaE[c],
                    phiN[faceI], phi[c], gradF
                );
        }
    }

    tmp<volVectorField> tgrad
    (
        new volVectorField
        (
            IOobject
            (
                "interfaceGrad(" + phi.name() + ")", m.time().timeName(), m
            ),
            fvc::grad(phi)
        )
    );
    vectorField& g = tgrad.ref().primitiveFieldRef();
    const scalarField& phiI = phi.primitiveField();

    labelList slot(m.nCells(), -1);
    label nTouched = 0;
    forAll(own, faceI)
    {
        if
        (
            (magSqr(sigmaI[own[faceI]]) > SMALL)
         != (magSqr(sigmaI[nei[faceI]]) > SMALL)
        )
        {
            for (const label c : {own[faceI], nei[faceI]})
            {
                if (slot[c] < 0)
                {
                    slot[c] = nTouched++;
                }
            }
        }
    }
    forAll(m.boundary(), patchI)
    {
        const fvPatch& p = m.boundary()[patchI];
        if (!p.coupled())
        {
            continue;
        }
        const tensorField sigmaIN
        (
            sigmaI.boundaryField()[patchI].patchNeighbourField()
        );
        forAll(p, faceI)
        {
            const label c = p.faceCells()[faceI];
            if
            (
                slot[c] < 0
             && (magSqr(sigmaI[c]) > SMALL)
             != (magSqr(sigmaIN[faceI]) > SMALL)
            )
            {
                slot[c] = nTouched++;
            }
        }
    }

    symmTensorField dd(nTouched, Zero);
    vectorField rhs(nTouched, Zero);
    const auto add = [&]
    (
        const label c,
        const vector& d,
        const scalar dPhi,
        const scalar area
    )
    {
        const scalar wt = area/magSqr(d);
        dd[slot[c]] += wt*sqr(d);
        rhs[slot[c]] += wt*dPhi*d;
    };

    forAll(own, faceI)
    {
        const label o = own[faceI];
        const label nb = nei[faceI];
        const bool jump =
            (magSqr(sigmaI[o]) > SMALL) != (magSqr(sigmaI[nb]) > SMALL);
        if (slot[o] >= 0)
        {
            jump
          ? add(o, Cf[faceI] - C[o], phiFI[faceI] - phiI[o], magSf[faceI])
          : add(o, C[nb] - C[o], phiI[nb] - phiI[o], magSf[faceI]);
        }
        if (slot[nb] >= 0)
        {
            jump
          ? add(nb, Cf[faceI] - C[nb], phiFI[faceI] - phiI[nb], magSf[faceI])
          : add(nb, C[o] - C[nb], phiI[o] - phiI[nb], magSf[faceI]);
        }
    }
    forAll(m.boundary(), patchI)
    {
        const fvPatch& p = m.boundary()[patchI];
        const labelUList& fc = p.faceCells();
        const vectorField nf(p.nf());
        const vectorField delta(p.delta());
        const scalarField& area = magSf.boundaryField()[patchI];
        const scalarField& phiFp = phiF.boundaryField()[patchI];
        if (p.coupled())
        {
            const scalarField phiN
            (
                phi.boundaryField()[patchI].patchNeighbourField()
            );
            const tensorField sigmaIN
            (
                sigmaI.boundaryField()[patchI].patchNeighbourField()
            );
            forAll(p, faceI)
            {
                const label c = fc[faceI];
                if (slot[c] < 0)
                {
                    continue;
                }
                (magSqr(sigmaI[c]) > SMALL)
             != (magSqr(sigmaIN[faceI]) > SMALL)
              ? add
                (
                    c, Cf.boundaryField()[patchI][faceI] - C[c],
                    phiFp[faceI] - phiI[c], area[faceI]
                )
              : add(c, delta[faceI], phiN[faceI] - phiI[c], area[faceI]);
            }
        }
        else
        {
            const scalarField& phiB = phi.boundaryField()[patchI];
            forAll(p, faceI)
            {
                const label c = fc[faceI];
                if (slot[c] >= 0)
                {
                    add
                    (
                        c, nf[faceI]*(nf[faceI] & delta[faceI]),
                        phiB[faceI] - phiI[c], area[faceI]
                    );
                }
            }
        }
    }

    const symmTensorField invDd(inv(dd));
    forAll(slot, c)
    {
        if (slot[c] >= 0)
        {
            g[c] = invDd[slot[c]] & rhs[slot[c]];
        }
    }
    tgrad.ref().correctBoundaryConditions();
    return tgrad;
}

}

// ************************************************************************* //
