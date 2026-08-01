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

Application
    bathBidomainInterfaceMetrics

Description
    Compute region-separated potential errors and one-sided heart--bath
    interface flux diagnostics on a reconstructed serial mesh.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "mathematicalConstants.H"
#include "HashSet.H"
#include "fvMeshSubset.H"
#include "extracellularFaceConductivity.H"
#include "conductivityFieldIO.H"

using namespace Foam;

namespace
{

struct WeightedNorms
{
    scalar weight{0.0};
    scalar sumAbs{0.0};
    scalar sumSq{0.0};
    scalar maxAbs{0.0};

    void add(const scalar error, const scalar w)
    {
        const scalar magnitude = mag(error);
        weight += w;
        sumAbs += w*magnitude;
        sumSq += w*sqr(magnitude);
        maxAbs = max(maxAbs, magnitude);
    }

    scalar l1() const { return sumAbs/max(weight, VSMALL); }
    scalar l2() const { return Foam::sqrt(sumSq/max(weight, VSMALL)); }
};


scalar phiEExact
(
    const scalar x,
    const scalar t,
    const scalar k,
    const scalar alpha,
    const scalar sigmaE
)
{
    const scalar sqrt1t = Foam::sqrt(1.0 + t);
    const scalar sigmaB = 0.5*sigmaE;
    const scalar C = k*sqrt1t + alpha/sigmaB;

    if (x <= 0.0)
    {
        return -k*sqrt1t + alpha*x/sigmaB + C;
    }
    if (x <= 1.0)
    {
        return
            -k*sqrt1t*Foam::cos(constant::mathematical::pi*x)
          + alpha*x/sigmaE + C;
    }

    return
        k*sqrt1t + alpha/sigmaE + alpha*(x - 1.0)/sigmaB + C;
}


void writeNorms(OFstream& output, const WeightedNorms& norms)
{
    output << ',' << norms.l1() << ',' << norms.l2() << ',' << norms.maxAbs;
}


struct FaceFit
{
    scalar value;
    scalar normalDerivative;
};


FaceFit oneSidedFaceFit
(
    const label cellI,
    const vector& faceCentre,
    const vector& outwardNormal,
    const boolList& isHeart,
    const volScalarField& field,
    const fvMesh& mesh
)
{
    HashSet<label> candidates;
    candidates.insert(cellI);

    // Two point-neighbour expansions provide at least three normal layers on
    // the structured reference and a well-populated one-sided tet stencil.
    for (label expansion = 0; expansion < 2; ++expansion)
    {
        const labelList current(candidates.toc());
        for (const label currentCell : current)
        {
            for (const label pointI : mesh.cellPoints()[currentCell])
            {
                for (const label candidate : mesh.pointCells()[pointI])
                {
                    if (isHeart[candidate] == isHeart[cellI])
                    {
                        candidates.insert(candidate);
                    }
                }
            }
        }
    }

    tensor normal(tensor::zero);
    vector rhs(vector::zero);
    forAllConstIter(HashSet<label>, candidates, iter)
    {
        const label candidate = iter.key();
        const scalar s =
            (mesh.C()[candidate] - faceCentre) & outwardNormal;
        const vector basis(1.0, s, sqr(s));

        normal.xx() += basis.x()*basis.x();
        normal.xy() += basis.x()*basis.y();
        normal.xz() += basis.x()*basis.z();
        normal.yx() += basis.y()*basis.x();
        normal.yy() += basis.y()*basis.y();
        normal.yz() += basis.y()*basis.z();
        normal.zx() += basis.z()*basis.x();
        normal.zy() += basis.z()*basis.y();
        normal.zz() += basis.z()*basis.z();
        rhs += basis*field[candidate];
    }

    if (mag(det(normal)) <= VSMALL)
    {
        FatalErrorInFunction
            << "Cannot construct a same-material quadratic face fit for cell "
            << cellI
            << " from " << candidates.size() << " point-connected cells."
            << exit(FatalError);
    }

    const vector coefficients(inv(normal) & rhs);
    return FaceFit{coefficients.x(), coefficients.y()};
}

} // End anonymous namespace


int main(int argc, char* argv[])
{
    argList::noParallel();
    timeSelector::addOptions();
    argList::addNote
    (
        "Report manufactured bath-bidomain region and interface metrics."
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    if (timeDirs.empty())
    {
        FatalErrorInFunction << "No time directory selected." << exit(FatalError);
    }
    runTime.setTime(timeDirs.last(), timeDirs.size() - 1);

    #include "createMesh.H"

    Info<< "Reading phiE" << nl;

    volScalarField phiE
    (
        IOobject
        (
            "phiE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    Info<< "Reading sigmaTotal" << nl;
    volTensorField sigmaTotal
    (
        IOobject
        (
            "sigmaTotal",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    Info<< "Reading VmGlobal" << nl;
    volScalarField VmGlobal
    (
        IOobject
        (
            "VmGlobal",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Info<< "Reading electroProperties" << nl;

    IOdictionary electroProperties
    (
        IOobject
        (
            "electroProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const dictionary& coeffs = electroProperties.subDict("bidomainSolverCoeffs");
    const dictionary& verification = coeffs.subDict("verificationModel");
    const dictionary& bath = coeffs.subDict("bathPotentialDomain");
    const scalar k = verification.lookupOrDefault<scalar>
    (
        "k",
        1.0/Foam::sqrt(2.0)
    );
    const scalar alpha = verification.lookupOrDefault<scalar>("alpha", 0.01);
    tmp<volTensorField> sigmaEFieldTmp = readConductivityField
    (
        mesh,
        mesh,
        nullptr,
        coeffs,
        conductivityFieldSpec
        {
            "ConductivityExtracellular",
            "conductivityExtracellular",
            "conductivityExtracellularDiagnostic"
        }
    );
    tmp<volTensorField> sigmaIFieldTmp = readConductivityField
    (
        mesh,
        mesh,
        nullptr,
        coeffs,
        conductivityFieldSpec
        {
            "ConductivityIntracellular",
            "conductivityIntracellular",
            "conductivityIntracellularDiagnostic"
        }
    );
    const tensor sigmaE(sigmaEFieldTmp().primitiveField()[0]);
    const tensor sigmaI(sigmaIFieldTmp().primitiveField()[0]);

    const scalar sigmaEVariation =
        gMax(mag(sigmaEFieldTmp().primitiveField() - sigmaE));
    const scalar sigmaIVariation =
        gMax(mag(sigmaIFieldTmp().primitiveField() - sigmaI));
    if (sigmaEVariation > SMALL || sigmaIVariation > SMALL)
    {
        FatalErrorInFunction
            << "bathBidomainInterfaceMetrics is a manufactured-solution "
            << "diagnostic and requires spatially uniform Gi and Ge. "
            << "max variations are Gi=" << sigmaIVariation
            << ", Ge=" << sigmaEVariation << '.'
            << exit(FatalError);
    }
    const word method = bath.lookupOrDefault<word>
    (
        "interfaceConductivityInterpolation",
        "unweightedHarmonic"
    );
    const word assembly("matchedSubmesh");


    volScalarField phiIHeart
    (
        IOobject
        (
            "phiIHeartDiagnostic",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phiE + VmGlobal
    );

    Info<< "Resolving cellZones" << nl;

    const label heartZone = mesh.cellZones().findZoneID("myocardium");
    const label bathZone = mesh.cellZones().findZoneID("bath");
    if (heartZone < 0 || bathZone < 0)
    {
        FatalErrorInFunction
            << "Required myocardium/bath cellZones are missing."
            << exit(FatalError);
    }

    boolList isHeart(mesh.nCells(), false);
    for (const label cellI : mesh.cellZones()[heartZone])
    {
        isHeart[cellI] = true;
    }

    surfaceTensorField sigmaExtracellularf
    (
        IOobject
        (
            "sigmaExtracellularfDiagnostic",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", sigmaTotal.dimensions(), tensor::zero)
    );
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const scalarField& weights = mesh.weights().primitiveField();
    const tensorField& sigma = sigmaTotal.primitiveField();
    tensorField& sigmaExtracellularFace =
        sigmaExtracellularf.primitiveFieldRef();
    forAll(neighbour, faceI)
    {
        const label ownCell = owner[faceI];
        const label neiCell = neighbour[faceI];
        const tensor sigmaIOwn = isHeart[ownCell] ? sigmaI : tensor::zero;
        const tensor sigmaINei = isHeart[neiCell] ? sigmaI : tensor::zero;
        const tensor ownSigmaE = sigma[ownCell] - sigmaIOwn;
        const tensor neiSigmaE = sigma[neiCell] - sigmaINei;

        if (isHeart[ownCell] == isHeart[neiCell])
        {
            sigmaExtracellularFace[faceI] =
                extracellularFaceConductivity::linear
                (
                    ownSigmaE, neiSigmaE, weights[faceI]
                );
        }
        else if (method == "unweightedHarmonic")
        {
            sigmaExtracellularFace[faceI] =
                extracellularFaceConductivity::unweightedHarmonic
                (
                    ownSigmaE, neiSigmaE
                );
        }
        else
        {
            sigmaExtracellularFace[faceI] =
                extracellularFaceConductivity::distanceWeightedHarmonic
                (
                    ownSigmaE,
                    neiSigmaE,
                    1.0 - weights[faceI],
                    weights[faceI]
                );
        }
    }
    forAll(sigmaExtracellularf.boundaryField(), patchI)
    {
        sigmaExtracellularf.boundaryFieldRef()[patchI] =
            sigmaTotal.boundaryField()[patchI].patchInternalField();
    }

    fvScalarMatrix assembledLaplacian
    (
        fvm::laplacian(sigmaExtracellularf, phiE)
    );
    if (!assembledLaplacian.hasFaceFluxCorrection())
    {
        FatalErrorInFunction
            << "phiE must be listed in fvSchemes fluxRequired to inspect the "
            << "non-orthogonal face-flux correction."
            << exit(FatalError);
    }
    const surfaceScalarField assembledFlux(assembledLaplacian.flux());
    const surfaceScalarField& correctionFlux =
        *assembledLaplacian.faceFluxCorrectionPtr();


    Info<< "Computing region norms" << nl;

    WeightedNorms heartPotential;
    WeightedNorms bathPotential;
    const scalarField& volumes = mesh.V();
    const vectorField& centres = mesh.C();
    forAll(phiE, cellI)
    {
        const scalar error =
            phiE[cellI]
          - phiEExact(centres[cellI].x(), runTime.value(), k, alpha, sigmaE.xx());
        (isHeart[cellI] ? heartPotential : bathPotential).add
        (
            error,
            volumes[cellI]
        );
    }

    const vectorField& faceCentres = mesh.Cf();
    const vectorField& faceAreas = mesh.Sf();

    WeightedNorms x0Potential;
    WeightedNorms x1Potential;
    WeightedNorms x0TraceJump;
    WeightedNorms x1TraceJump;
    WeightedNorms x0HeartFlux;
    WeightedNorms x0BathFlux;
    WeightedNorms x0FluxJump;
    WeightedNorms x1HeartFlux;
    WeightedNorms x1BathFlux;
    WeightedNorms x1FluxJump;
    WeightedNorms x0IntracellularLeak;
    WeightedNorms x1IntracellularLeak;
    WeightedNorms x0AssembledFlux;
    WeightedNorms x1AssembledFlux;
    WeightedNorms x0CorrectionFlux;
    WeightedNorms x1CorrectionFlux;
    WeightedNorms x0OrthogonalFlux;
    WeightedNorms x1OrthogonalFlux;
    scalar x0IntracellularLeakIntegral = 0.0;
    scalar x1IntracellularLeakIntegral = 0.0;
    label interfaceFaces = 0;

    Info<< "Computing interface metrics" << nl;

    forAll(neighbour, faceI)
    {
        const label ownCell = owner[faceI];
        const label neiCell = neighbour[faceI];
        if (isHeart[ownCell] == isHeart[neiCell])
        {
            continue;
        }

        ++interfaceFaces;
        const label heartCell = isHeart[ownCell] ? ownCell : neiCell;
        const label bathCell = isHeart[ownCell] ? neiCell : ownCell;
        const vector nHeart =
            (isHeart[ownCell] ? 1.0 : -1.0)*faceAreas[faceI]
           /mag(faceAreas[faceI]);
        const scalar area = mag(faceAreas[faceI]);
        const FaceFit heartFit = oneSidedFaceFit
        (
            heartCell,
            faceCentres[faceI],
            nHeart,
            isHeart,
            phiE,
            mesh
        );
        const FaceFit bathFit = oneSidedFaceFit
        (
            bathCell,
            faceCentres[faceI],
            -nHeart,
            isHeart,
            phiE,
            mesh
        );
        const FaceFit phiIFit = oneSidedFaceFit
        (
            heartCell,
            faceCentres[faceI],
            nHeart,
            isHeart,
            phiIHeart,
            mesh
        );
        const scalar qHeart =
            (nHeart & (sigmaE & nHeart))*heartFit.normalDerivative;
        const scalar qBath =
            ((-nHeart) & (sigma[bathCell] & (-nHeart)))
           *bathFit.normalDerivative;
        const scalar qIntracellular =
            (nHeart & (sigmaI & nHeart))*phiIFit.normalDerivative;
        const bool leftInterface = faceCentres[faceI].x() < 0.5;
        const scalar exactHeartFlux = leftInterface ? -alpha : alpha;
        const scalar heartOrientation = isHeart[ownCell] ? 1.0 : -1.0;
        const scalar qAssembled =
            heartOrientation*assembledFlux[faceI]/area;
        const scalar qCorrection =
            heartOrientation*correctionFlux[faceI]/area;
        const scalar qOrthogonal = qAssembled - qCorrection;
        const scalar facePotential =
            weights[faceI]*phiE[ownCell]
          + (1.0 - weights[faceI])*phiE[neiCell];
        const scalar exactPotential = phiEExact
        (
            faceCentres[faceI].x(),
            runTime.value(),
            k,
            alpha,
            sigmaE.xx()
        );
        const scalar heartTrace = heartFit.value;
        const scalar bathTrace = bathFit.value;

        WeightedNorms& potential = leftInterface ? x0Potential : x1Potential;
        WeightedNorms& traceJump = leftInterface ? x0TraceJump : x1TraceJump;
        WeightedNorms& heartFlux = leftInterface ? x0HeartFlux : x1HeartFlux;
        WeightedNorms& bathFlux = leftInterface ? x0BathFlux : x1BathFlux;
        WeightedNorms& fluxJump = leftInterface ? x0FluxJump : x1FluxJump;
        WeightedNorms& intracellularLeak =
            leftInterface ? x0IntracellularLeak : x1IntracellularLeak;
        WeightedNorms& assembledFluxError =
            leftInterface ? x0AssembledFlux : x1AssembledFlux;
        WeightedNorms& correctionFluxMagnitude =
            leftInterface ? x0CorrectionFlux : x1CorrectionFlux;
        WeightedNorms& orthogonalFluxError =
            leftInterface ? x0OrthogonalFlux : x1OrthogonalFlux;

        potential.add(facePotential - exactPotential, area);
        traceJump.add(heartTrace - bathTrace, area);
        heartFlux.add(qHeart - exactHeartFlux, area);
        bathFlux.add(qBath + exactHeartFlux, area);
        fluxJump.add(qHeart + qBath, area);
        intracellularLeak.add(qIntracellular, area);
        assembledFluxError.add(qAssembled - exactHeartFlux, area);
        correctionFluxMagnitude.add(qCorrection, area);
        orthogonalFluxError.add(qOrthogonal - exactHeartFlux, area);
        if (leftInterface)
        {
            x0IntracellularLeakIntegral += qIntracellular*area;
        }
        else
        {
            x1IntracellularLeakIntegral += qIntracellular*area;
        }
    }

    WeightedNorms xMinPotential;
    WeightedNorms xMaxFlux;
    WeightedNorms sidesFlux;
    scalar exteriorFluxIntegral = 0.0;
    Info<< "Computing boundary metrics" << nl;
    forAll(mesh.boundary(), patchI)
    {
        const word& patchName = mesh.boundary()[patchI].name();
        const vectorField normals(mesh.boundary()[patchI].nf());
        const scalarField areas(mesh.boundary()[patchI].magSf());
        const scalarField snGrad(phiE.boundaryField()[patchI].snGrad());
        const tensorField patchSigma
        (
            sigmaTotal.boundaryField()[patchI].patchInternalField()
        );

        forAll(areas, faceI)
        {
            const scalar flux =
                normals[faceI]
              & (patchSigma[faceI] & (snGrad[faceI]*normals[faceI]));
            exteriorFluxIntegral += flux*areas[faceI];

            if (patchName == "xMin")
            {
                xMinPotential.add(phiE.boundaryField()[patchI][faceI], areas[faceI]);
            }
            else if (patchName == "xMax")
            {
                xMaxFlux.add(flux - alpha, areas[faceI]);
            }
            else if (patchName == "sides")
            {
                sidesFlux.add(flux, areas[faceI]);
            }
        }
    }

    const fileName outputDir(runTime.globalPath()/"postProcessing");
    mkDir(outputDir);
    OFstream output(outputDir/"bathBidomainInterfaceMetrics.csv");
    output
        << "method,assembly,fieldSource,time,interfaceFaces"
        << ",heartPhiE_L1,heartPhiE_L2,heartPhiE_Linf"
        << ",bathPhiE_L1,bathPhiE_L2,bathPhiE_Linf"
        << ",x0PhiE_L1,x0PhiE_L2,x0PhiE_Linf"
        << ",x1PhiE_L1,x1PhiE_L2,x1PhiE_Linf"
        << ",x0TraceJump_L1,x0TraceJump_L2,x0TraceJump_Linf"
        << ",x1TraceJump_L1,x1TraceJump_L2,x1TraceJump_Linf"
        << ",x0HeartFlux_L1,x0HeartFlux_L2,x0HeartFlux_Linf"
        << ",x0BathFlux_L1,x0BathFlux_L2,x0BathFlux_Linf"
        << ",x0FluxJump_L1,x0FluxJump_L2,x0FluxJump_Linf"
        << ",x1HeartFlux_L1,x1HeartFlux_L2,x1HeartFlux_Linf"
        << ",x1BathFlux_L1,x1BathFlux_L2,x1BathFlux_Linf"
        << ",x1FluxJump_L1,x1FluxJump_L2,x1FluxJump_Linf"
        << ",x0IntracellularLeak_L1,x0IntracellularLeak_L2,x0IntracellularLeak_Linf"
        << ",x1IntracellularLeak_L1,x1IntracellularLeak_L2,x1IntracellularLeak_Linf"
        << ",x0IntracellularLeakIntegral,x1IntracellularLeakIntegral"
        << ",x0AssembledFlux_L1,x0AssembledFlux_L2,x0AssembledFlux_Linf"
        << ",x1AssembledFlux_L1,x1AssembledFlux_L2,x1AssembledFlux_Linf"
        << ",x0CorrectionFlux_L1,x0CorrectionFlux_L2,x0CorrectionFlux_Linf"
        << ",x1CorrectionFlux_L1,x1CorrectionFlux_L2,x1CorrectionFlux_Linf"
        << ",x0OrthogonalFlux_L1,x0OrthogonalFlux_L2,x0OrthogonalFlux_Linf"
        << ",x1OrthogonalFlux_L1,x1OrthogonalFlux_L2,x1OrthogonalFlux_Linf"
        << ",xMinPhiE_L1,xMinPhiE_L2,xMinPhiE_Linf"
        << ",xMaxFlux_L1,xMaxFlux_L2,xMaxFlux_Linf"
        << ",sidesFlux_L1,sidesFlux_L2,sidesFlux_Linf"
        << ",exteriorFluxIntegral\n";
    output
        << method << ',' << assembly << ','
        << "numerical" << ','
        << runTime.value() << ',' << interfaceFaces;
    writeNorms(output, heartPotential);
    writeNorms(output, bathPotential);
    writeNorms(output, x0Potential);
    writeNorms(output, x1Potential);
    writeNorms(output, x0TraceJump);
    writeNorms(output, x1TraceJump);
    writeNorms(output, x0HeartFlux);
    writeNorms(output, x0BathFlux);
    writeNorms(output, x0FluxJump);
    writeNorms(output, x1HeartFlux);
    writeNorms(output, x1BathFlux);
    writeNorms(output, x1FluxJump);
    writeNorms(output, x0IntracellularLeak);
    writeNorms(output, x1IntracellularLeak);
    output
        << ',' << x0IntracellularLeakIntegral
        << ',' << x1IntracellularLeakIntegral;
    writeNorms(output, x0AssembledFlux);
    writeNorms(output, x1AssembledFlux);
    writeNorms(output, x0CorrectionFlux);
    writeNorms(output, x1CorrectionFlux);
    writeNorms(output, x0OrthogonalFlux);
    writeNorms(output, x1OrthogonalFlux);
    writeNorms(output, xMinPotential);
    writeNorms(output, xMaxFlux);
    writeNorms(output, sidesFlux);
    output << ',' << exteriorFluxIntegral << nl;

    Info<< "Wrote " << outputDir/"bathBidomainInterfaceMetrics.csv" << nl
        << "interfaceFaces=" << interfaceFaces << nl
        << "heartPhiE_L2=" << heartPotential.l2() << nl
        << "bathPhiE_L2=" << bathPotential.l2() << nl
        << "exteriorFluxIntegral=" << exteriorFluxIntegral << nl
        << "End" << nl << endl;

    return 0;
}

// ************************************************************************* //
