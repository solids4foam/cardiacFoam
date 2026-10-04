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

Function
    exposedPhiETrace

Description
    Heart phiE exposed-face trace from the global phiE.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

#include "exposedPhiETrace.H"
#include "fixedValueFvPatchFields.H"
#include "surfaceFields.H"
#include "fvc.H"

namespace Foam
{

namespace
{

bool exposedBaseFace(const fvMesh& baseMesh, const label f)
{
    return
        f < baseMesh.nInternalFaces()
     || baseMesh.boundary()[baseMesh.boundaryMesh().whichPatch(f)].coupled();
}



}


scalar interfaceValue
(
    const vector& n,
    const point& Cf,
    const point& CH,
    const point& CB,
    const tensor& aH,
    const tensor& aB,
    const scalar phiH,
    const scalar phiB,
    const vector& gradF
)
{
    const scalar dH = mag((Cf - CH) & n);
    const scalar dB = mag((Cf - CB) & n);
    const scalar aHNN = n & aH & n;
    const scalar aBNN = n & aB & n;
    const scalar kH = aHNN/max(dH, SMALL);
    const scalar kB = aBNN/max(dB, SMALL);
    if (kH + kB <= VSMALL)
    {
        return 0.5*(phiH + phiB);
    }
    const scalar tH = ((aH & n) - aHNN*n) & gradF;
    const scalar tB = ((aB & n) - aBNN*n) & gradF;
    const point yH(CH + n*(n & (Cf - CH)));
    const point yB(CB + n*(n & (Cf - CB)));
    const vector offset(Cf - (kH*yH + kB*yB)/(kH + kB));
    return
        (kH*phiH + kB*phiB + tB - tH)/(kH + kB)
      + (gradF & offset);
}


label useExposedPhiETrace
(
    volScalarField& heartPhiE,
    const labelUList& faceMap,
    const fvMesh& baseMesh
)
{
    label nExposed = 0;
    volScalarField::Boundary& bf = heartPhiE.boundaryFieldRef();
    forAll(bf, patchI)
    {
        const fvPatch& p = heartPhiE.mesh().boundary()[patchI];
        if (p.coupled())
        {
            continue;
        }
        bool exposed = false;
        forAll(p, i)
        {
            exposed =
                exposed || exposedBaseFace(baseMesh, faceMap[p.start() + i]);
        }
        if (exposed)
        {
            bf.set
            (
                patchI,
                fvPatchField<scalar>::New
                (
                    fixedValueFvPatchScalarField::typeName, p, heartPhiE
                )
            );
            bf[patchI] = bf[patchI].patchInternalField();
            nExposed += p.size();
        }
    }
    return returnReduce(nExposed, sumOp<label>());
}


void setExposedPhiETrace
(
    volScalarField& heartPhiE,
    const volScalarField& globalPhiE,
    const labelUList& cellMap,
    const labelUList& faceMap,
    const volTensorField& Ge,
    const volTensorField& sigmaGlobal
)
{
    scalarField& heartI = heartPhiE.primitiveFieldRef();
    forAll(cellMap, cellI)
    {
        heartI[cellI] = globalPhiE[cellMap[cellI]];
    }

    bool anyTrace = false;
    forAll(heartPhiE.boundaryField(), patchI)
    {
        anyTrace =
            anyTrace
         || isA<fixedValueFvPatchScalarField>
            (
                heartPhiE.boundaryField()[patchI]
            );
    }
    if (!returnReduce(anyTrace, orOp<bool>()))
    {
        heartPhiE.correctBoundaryConditions();
        return;
    }

    const fvMesh& baseMesh = globalPhiE.mesh();
    const volVectorField& C = baseMesh.C();
    const surfaceVectorField& Cf = baseMesh.Cf();
    const surfaceVectorField& Sf = baseMesh.Sf();
    const surfaceScalarField& magSf = baseMesh.magSf();
    const labelUList& own = baseMesh.owner();
    const labelUList& nei = baseMesh.neighbour();
    const polyBoundaryMesh& basePatches = baseMesh.boundaryMesh();

    const volVectorField gradG(fvc::grad(globalPhiE));
    PtrList<scalarField> phiNbr(basePatches.size());
    PtrList<tensorField> sigmaNbr(basePatches.size());
    PtrList<vectorField> deltaNbr(basePatches.size());
    PtrList<vectorField> nfNbr(basePatches.size());
    PtrList<vectorField> gradNbr(basePatches.size());
    forAll(basePatches, bp)
    {
        if (baseMesh.boundary()[bp].coupled())
        {
            phiNbr.set
            (
                bp, globalPhiE.boundaryField()[bp].patchNeighbourField()
            );
            sigmaNbr.set
            (
                bp, sigmaGlobal.boundaryField()[bp].patchNeighbourField()
            );
            deltaNbr.set(bp, baseMesh.boundary()[bp].delta());
            nfNbr.set(bp, baseMesh.boundary()[bp].nf());
            gradNbr.set(bp, gradG.boundaryField()[bp].patchNeighbourField());
        }
    }

    volScalarField::Boundary& bf = heartPhiE.boundaryFieldRef();
    forAll(bf, patchI)
    {
        if (!isA<fixedValueFvPatchScalarField>(bf[patchI]))
        {
            continue;
        }
        const fvPatch& hp = heartPhiE.mesh().boundary()[patchI];
        scalarField& trace = bf[patchI];
        forAll(trace, i)
        {
            const label f = faceMap[hp.start() + i];
            const label hc = hp.faceCells()[i];
            const label hb = cellMap[hc];
            if (f < baseMesh.nInternalFaces())
            {
                const label bb = own[f] == hb ? nei[f] : own[f];
                const vector n((own[f] == hb ? 1 : -1)*Sf[f]/magSf[f]);
                const scalar wf = baseMesh.weights()[f];
                trace[i] = interfaceValue
                (
                    n, Cf[f], C[hb], C[bb],
                    Ge[hc], sigmaGlobal[bb],
                    globalPhiE[hb], globalPhiE[bb],
                    own[f] == hb
                  ? wf*gradG[hb] + (1 - wf)*gradG[bb]
                  : wf*gradG[bb] + (1 - wf)*gradG[hb]
                );
                continue;
            }
            const label bp = basePatches.whichPatch(f);
            const label pf = f - basePatches[bp].start();
            if (phiNbr.set(bp))
            {
                const vector& n = nfNbr[bp][pf];
                const scalar wp = baseMesh.weights().boundaryField()[bp][pf];
                trace[i] = interfaceValue
                (
                    n,
                    Cf.boundaryField()[bp][pf],
                    C[hb],
                    C[hb] + deltaNbr[bp][pf],
                    Ge[hc], sigmaNbr[bp][pf],
                    globalPhiE[hb], phiNbr[bp][pf],
                    wp*gradG[hb] + (1 - wp)*gradNbr[bp][pf]
                );
            }
            else
            {
                trace[i] = globalPhiE.boundaryField()[bp][pf];
            }
        }
    }
    heartPhiE.correctBoundaryConditions();
}

}

// ************************************************************************* //
