/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "greensFunctionECGElectro.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace electroModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(greensFunctionECGElectro, 0);
addToRunTimeSelectionTable(electroModel, greensFunctionECGElectro, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

tmp<volTensorField> greensFunctionECGElectro::initialiseGi() const
{
    tmp<volTensorField> tresult
    (
        new volTensorField
        (
            IOobject
            (
                "Gi",
                runTime().timeName(),
                mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor
            (
                "zero",
                pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                tensor::zero
            )
        )
    );
    volTensorField& result = tresult.ref();

    if (!result.headerOk())
    {
        Info<< "\nGi not found on disk, using Gi from ECG subdict" << nl
            << endl;

        const dictionary& ecgDict = electroProperties().subDict("ECG");

        result =
            dimensionedTensor
            (
                dimensionedSymmTensor
                (
                    "Gi",
                    pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
                    ecgDict
                ) & tensor(I)
            );

        if (result.size() > 0)
        {
            Info<< "Gi tensor (cell 0): " << result[0] << nl;
        }
    }
    else
    {
        Info<< "Gi field read from " << runTime().timeName() << nl << endl;
    }

    return tresult;
}


void greensFunctionECGElectro::readElectrodes(const dictionary& ecgDict)
{
    const dictionary& eDict = ecgDict.subDict("electrodes");

    // toc() returns keys in insertion order
    const wordList names(eDict.toc());

    electrodeNames_.setSize(names.size());
    electrodePositions_.setSize(names.size());

    forAll(names, i)
    {
        electrodeNames_[i]     = names[i];
        electrodePositions_[i] = eDict.get<vector>(names[i]);
    }

    Info<< "ECG electrodes (" << electrodeNames_.size() << "):" << nl;
    forAll(electrodeNames_, i)
    {
        Info<< "  " << electrodeNames_[i]
            << "  @  " << electrodePositions_[i] << nl;
    }
    Info<< endl;
}


void greensFunctionECGElectro::writeHeader()
{
    OFstream& os = outputPtr_.ref();
    os.setf(std::ios::scientific);
    os.precision(8);

    os << "# time";
    forAll(electrodeNames_, i)
    {
        os << "  " << electrodeNames_[i];
    }
    os << nl;
}


void greensFunctionECGElectro::loadTorsoSurface(const dictionary& ecgDict)
{
    const fileName stlPath(ecgDict.get<fileName>("torsoSurface"));
    const fileName fullPath(runTime().path()/stlPath);

    Info<< "greensFunctionECGElectro: loading torso surface from '"
        << fullPath << "'" << nl << endl;

    torsoSurfacePtr_.reset(new triSurface(fullPath));
    torsoFaceCentres_ = torsoSurfacePtr_().faceCentres();

    Info<< "Torso surface: " << torsoFaceCentres_.size()
        << " faces loaded." << nl << endl;

    torsoVtkDir_ = runTime().path()/"postProcessing"/"ECG";
    mkDir(torsoVtkDir_);
}


void greensFunctionECGElectro::writeTorsoVtk(const scalarList& phiTorso) const
{
    const triSurface& surf = torsoSurfacePtr_();
    const pointField& pts  = surf.points();
    const label nPoints    = pts.size();
    const label nTris      = surf.size();

    const fileName vtkFile
    (
        torsoVtkDir_/"phiT_" + runTime().timeName() + ".vtk"
    );

    OFstream os(vtkFile);
    os.setf(std::ios::scientific);
    os.precision(8);

    os << "# vtk DataFile Version 3.0\n";
    os << "Torso ECG phiT t=" << runTime().timeName() << "\n";
    os << "ASCII\n";
    os << "DATASET POLYDATA\n";

    os << "POINTS " << nPoints << " float\n";
    forAll(pts, pI)
    {
        os << pts[pI].x() << " " << pts[pI].y() << " " << pts[pI].z() << "\n";
    }

    os << "POLYGONS " << nTris << " " << 4*nTris << "\n";
    forAll(surf, tI)
    {
        const triFace& tri = surf[tI];
        os << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }

    os << "CELL_DATA " << nTris << "\n";
    os << "SCALARS phiT float 1\n";
    os << "LOOKUP_TABLE default\n";
    forAll(phiTorso, fI)
    {
        os << phiTorso[fI] << "\n";
    }

    Info<< "ECG surface written to " << vtkFile << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

greensFunctionECGElectro::greensFunctionECGElectro
(
    Time& runTime,
    const word& region
)
:
    monoDomainElectro(typeName, runTime, region),
    Is_
    (
        IOobject
        (
            "Is",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimCurrent/pow3(dimLength), 0.0),
        "zeroGradient"
    ),
    Gi_(initialiseGi()),
    sigmaT_
    (
        "sigmaT",
        pow3(dimTime)*sqr(dimCurrent)/(dimMass*dimVolume),
        electroProperties().subDict("ECG")
    ),
    electrodeNames_(),
    electrodePositions_(),
    outputPtr_(),
    torsoSurfacePtr_(),
    torsoFaceCentres_(),
    torsoVtkDir_()
{
    const dictionary& ecgDict = electroProperties().subDict("ECG");

    const bool hasElectrodes = ecgDict.found("electrodes");
    const bool hasSurface    = ecgDict.found("torsoSurface");

    if (!hasElectrodes && !hasSurface)
    {
        FatalErrorInFunction
            << "ECG subdict must contain 'electrodes', 'torsoSurface', or both."
            << abort(FatalError);
    }

    if (hasElectrodes)
    {
        readElectrodes(ecgDict);
        const fileName outDir(runTime.path()/"postProcessing");
        mkDir(outDir);
        outputPtr_.reset(new OFstream(outDir/"ECG.dat"));
        writeHeader();
    }

    if (hasSurface)
    {
        loadTorsoSurface(ecgDict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool greensFunctionECGElectro::evolve()
{
    // 1) Advance monodomain — updates Vm
    monoDomainElectro::evolve();

    // 2) Current source density: Is = -div( Gi . grad(Vm) )
    Is_ = -fvc::div(Gi_ & fvc::grad(Vm()));
    Is_.correctBoundaryConditions();

    // 3) Electrode potentials via Green's function volume integral:
    //    phi(P) = 1/(4*pi*sigmaT) * sum_c [ Is_c * V_c / |C_c - P| ]
    const scalarField& IsI  = Is_.primitiveField();
    const scalarField& Vols = mesh().V();
    const vectorField& Ctrs = mesh().C().primitiveField();
    const scalar invCoeff =
        1.0/(4.0*constant::mathematical::pi*sigmaT_.value());

    scalarList phiElectrodes(electrodePositions_.size(), scalar(0));

    forAll(electrodePositions_, eI)
    {
        const vector& P = electrodePositions_[eI];
        scalar localSum = 0.0;

        forAll(Ctrs, cI)
        {
            const scalar r = mag(Ctrs[cI] - P);
            if (r > VSMALL)
            {
                localSum += IsI[cI]*Vols[cI]/r;
            }
        }

        // Sum across all processors (parallel-safe)
        reduce(localSum, sumOp<scalar>());
        phiElectrodes[eI] = invCoeff*localSum;
    }

    // 4) Write one row to ECG.dat every time step
    OFstream& os = outputPtr_.ref();
    os << runTime().value();
    forAll(phiElectrodes, eI)
    {
        os << "  " << phiElectrodes[eI];
    }
    os << nl;

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electroModels
} // End namespace Foam

// ************************************************************************* //
