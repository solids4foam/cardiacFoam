/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    Foam

Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "openHFDIBDEM.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#include "scalarMatrices.H"
#include "OFstream.H"
#include <iostream>
#include "defineExternVars.H"
// #include "solverInfo.H"
//#include "virtualMeshLevel.H"
//#include "parameters.H"

#define ORDER 2

using namespace Foam;
//using namespace contactModel;

//---------------------------------------------------------------------------//
openHFDIBDEM::openHFDIBDEM(const Foam::fvMesh& mesh)
:
mesh_(mesh),
HFDIBDEMDict_
(
    IOobject
    (
        "HFDIBDEMDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
transportProperties_
(
    IOobject
    (
        "transportProperties",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
bodyNames_(HFDIBDEMDict_.lookup("bodyNames")),
stepDEM_(readScalar(HFDIBDEMDict_.lookup("stepDEM"))),
recordSimulation_(readBool(HFDIBDEMDict_.lookup("recordSimulation")))
{
    if(HFDIBDEMDict_.found("recordFirstTimeStep"))
    {
        recordFirstTimeStep_ = readBool(HFDIBDEMDict_.lookup("recordFirstTimeStep"));
    }

    dictionary demDic = HFDIBDEMDict_.subDict("DEM");
    (void)demDic;

    if (HFDIBDEMDict_.found("geometricD"))
    {
        geometricD = vector(HFDIBDEMDict_.lookup("geometricD"));
    }
    else
    {
        geometricD = mesh_.geometricD();
    }

    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDir[direction] = 1;
            emptyDim = direction;
            break;
        }
    }

    recordOutDir_ = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/bodiesInfo";
}
//---------------------------------------------------------------------------//
openHFDIBDEM::~openHFDIBDEM()
{}
//---------------------------------------------------------------------------//
void openHFDIBDEM::initialize
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF,
    label recomputeM0,
    word runTime
)
{
    if(HFDIBDEMDict_.found("outputSetup"))
    {
        dictionary outputDic = HFDIBDEMDict_.subDict("outputSetup");
        bool basicOutput = readBool(outputDic.lookup("basic"));
        bool iBoutput = readBool(outputDic.lookup("iB"));
        bool DEMoutput = readBool(outputDic.lookup("DEM"));
        bool addModelOutput = readBool(outputDic.lookup("addModel"));
        bool parallelDEMOutput = readBool(outputDic.lookup("parallelDEM"));
        writeIBPoints_ =
            outputDic.lookupOrDefault<Switch>("writeIBPoints", false);
        InfoH.setOutput(
            basicOutput,
            iBoutput,
            DEMoutput,
            addModelOutput,
            parallelDEMOutput
        );
    }

    preCalculateCellPoints();

    if(HFDIBDEMDict_.found("interpolationSchemes"))
    {
        HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");

        if(HFDIBinterpDict_.found("method"))
        {
            const word intMethod(HFDIBinterpDict_.lookup("method"));

            if(intMethod == "leastSquares")
            {
                const dictionary& lsCoeffsDict =
                    HFDIBinterpDict_.subDict("leastSquaresCoeffs");
                ibInterp_.reset
                (
                    new leastSquaresInt
                    (
                        mesh_,
                        readScalar(lsCoeffsDict.lookup("distFactor")),
                        readScalar(lsCoeffsDict.lookup("radiusFactor")),
                        readScalar(lsCoeffsDict.lookup("angleFactor")),
                        readScalar(lsCoeffsDict.lookup("maxCCRows"))
                    )
                );
            }
            else if(intMethod == "line")
            {
                ibInterp_.reset(new lineInt(HFDIBinterpDict_));
            }
        }
    }

    bool startTime0(runTime == "0");

    // initialize addModels
    addModels_.setSize(bodyNames_.size());
    immersedBodies_.setSize(0);                                         //on the fly creation
    refineF *= 0.0;
    recomputeM0_ = recomputeM0;

    if(!startTime0)
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            fileNameList entries(readDir(recordOutDir_,fileName::Type::DIRECTORY));
            scalar runTimeS(stod(runTime));
            forAll(entries,entry)
            {
                scalar dirTime(stod(entries[entry].name()));
                if(dirTime > runTimeS)
                {
                    word pathI(recordOutDir_ + "/" + entries[entry]);
                    rmDir(pathI);
                }
            }
        }

        restartSimulation(body, refineF, runTime);
    }
    else
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            rmDir(recordOutDir_);
            mkDir(recordOutDir_);
        }
    }

    #include "initializeAddModels.H"

    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);
        InfoH << basic_Info << "Creating immersed body based on: " << bodyName << endl;

        label maxAdditions(1000);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions) // and immersedBodies_.size() < solverInfo::getNSolidsTreshnold())
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));
            cAddition++;

            // initialize the immersed bodies
            if (addModels_[modelI].getBodyAdded())
            {
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                InfoH << addModel_Info << "Trying to set immersedBodies" << endl;
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );
                immersedBodies_[addIBPos].createImmersedBody(body,refineF);
                immersedBodies_[addIBPos].computeBodyCharPars();
                if (immersedBodies_[addIBPos].getStartSynced())
                {
                    immersedBodies_[addIBPos].initSyncWithFlow(U);
                }
                InfoH << addModel_Info << "Body based on: " << bodyName << " successfully added" << endl;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info << "Body based on: "
                    << bodyName << " should have been added but was not "
                    << "(probably overlap with an already existing body)"
                    << endl;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::createBodies(volScalarField& body,volScalarField& refineF)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].postContactUpdateBodyField(body,refineF);
        }
    }

    DynamicList<scalar> particleMasses;
    DynamicList<label> particleCells;
    DynamicList<symmTensor> particleInertiaTensors;

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncImmersedBodyParralell1(body,refineF);
            {
                particleMasses.append(immersedBodies_[bodyId].getGeomModel().getM());
                particleCells.append(immersedBodies_[bodyId].getGeomModel().getNCells());
                particleInertiaTensors.append(immersedBodies_[bodyId].getGeomModel().getI());
            }
        }
    }
    reduce(particleMasses,sumOp<List<scalar>>());
    reduce(particleCells,sumOp<List<label>>());
    reduce(particleInertiaTensors,sumOp<List<symmTensor>>());

    label bodyIndex(0);
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            {
                immersedBodies_[bodyId].getGeomModel().setM(particleMasses[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setNCells(particleCells[bodyIndex]);
                immersedBodies_[bodyId].getGeomModel().setI(particleInertiaTensors[bodyIndex]);
                bodyIndex++;
            }

            immersedBodies_[bodyId].syncImmersedBodyParralell2(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
            immersedBodies_[bodyId].updateOldMovementVars();
        }
    }

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].chceckBodyOp();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // create body or compute body-fluid coupling and estimate
            // potential contacts with walls
            immersedBodies_[bodyId].inContactWithStatic(false);

            immersedBodies_[bodyId].updateOldMovementVars();
            immersedBodies_[bodyId].printStats();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::postUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].clearIntpInfo();
            immersedBodies_[bodyId].postPimpleUpdateImmersedBody(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::recreateBodies
(
    volScalarField& body,
    volScalarField& refineF
)
{
    refineF *= 0.0;
    preCalculateCellPoints();
    forAll (addModels_,modelI)
    {
        addModels_[modelI].recreateBoundBox();
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].recreateBodyField(body,refineF);
        }
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncCreateImmersedBody(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
            if(immersedBodies_[bodyId].getrecomputeM0() > 0)
            {
                immersedBodies_[bodyId].computeBodyCharPars();
                immersedBodies_[bodyId].recomputedM0();
            }
            InfoH << iB_Info << "-- body "
                << immersedBodies_[bodyId].getBodyId() << " Re-created" << endl;
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::interpolateIB
(
    volVectorField& V,
    volVectorField& Vs,
    volScalarField& body
)
{
    if (ibInterp_.valid())
    {
        ibInterp_->resetInterpolator(V);
    }

    // Reset imposed field
    Vs = V;

    // Loop over all immersed bodies
    forAll(immersedBodies_, bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // Update imposed field according to body
            immersedBodies_[bodyId].updateVectorField(Vs, V.name(), body);

            if (ibInterp_.valid())
            {
                vectorField uIb = immersedBodies_[bodyId].getUatIbPoints();

                // Optional diagnostic output shared by lineInt and
                // leastSquares. uIb is evaluated immediately above so the VTK
                // file contains the target velocity used this write time.
                if (writeIBPoints_ && mesh_.time().writeTime())
                {
                    writeIBPoints
                    (
                        bodyId,
                        immersedBodies_[bodyId].getIntpInfo(),
                        uIb
                    );
                }

                ibInterp_->ibInterpolate
                (
                    immersedBodies_[bodyId].getIntpInfo(),
                    Vs,
                    uIb,
                    mesh_
                );
            }
        }
    }
}
//--------------------------------------------------------------------------//
void openHFDIBDEM::writeIBPoints
(
    const label bodyId,
    interpolationInfo& ibInfo,
    const vectorField& ibPointsVal
) const
{
    const List<point>& ibPts = ibInfo.getIbPoints();
    const List<vector>& ibNormals = ibInfo.getIbNormals();
    const auto& surfCells = ibInfo.getSurfCells();

    if (ibPts.empty())
    {
        return;
    }

    const word& bodyName = immersedBodies_[bodyId].getBodyName();

    fileName outDir =
        mesh_.time().rootPath()
      / mesh_.time().globalCaseName()
      / "IBPoints"
      / bodyName
      / mesh_.time().timeName();

    mkDir(outDir);

    fileName fName =
        outDir
      / (
            "ibPoints_body"
          + Foam::name(immersedBodies_[bodyId].getBodyId())
          + "_proc"
          + Foam::name(Pstream::myProcNo())
          + ".vtk"
        );

    OFstream os(fName);

    os  << "# vtk DataFile Version 3.0\n"
        << "HFDIB IB points body=" << bodyName
        << " bodyId=" << immersedBodies_[bodyId].getBodyId()
        << " time=" << mesh_.time().value()
        << " proc=" << Pstream::myProcNo() << "\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << ibPts.size() << " float\n";

    forAll(ibPts, i)
    {
        os << ibPts[i].x() << " "
           << ibPts[i].y() << " "
           << ibPts[i].z() << "\n";
    }

    os  << "\nVERTICES " << ibPts.size() << " "
        << 2*ibPts.size() << "\n";

    forAll(ibPts, i)
    {
        os << "1 " << i << "\n";
    }

    os  << "\nPOINT_DATA " << ibPts.size() << "\n"
        << "VECTORS ibNormal float\n";

    forAll(ibPts, i)
    {
        const vector n =
            (i < ibNormals.size() ? ibNormals[i] : vector::zero);
        os << n.x() << " " << n.y() << " " << n.z() << "\n";
    }

    os << "\nVECTORS ibTargetVelocity float\n";

    forAll(ibPts, i)
    {
        const vector u =
            (i < ibPointsVal.size() ? ibPointsVal[i] : vector::zero);
        os << u.x() << " " << u.y() << " " << u.z() << "\n";
    }

    os  << "\nSCALARS bodyId int 1\n"
        << "LOOKUP_TABLE default\n";

    forAll(ibPts, i)
    {
        os << immersedBodies_[bodyId].getBodyId() << "\n";
    }

    os  << "\nSCALARS surfaceCellId int 1\n"
        << "LOOKUP_TABLE default\n";

    forAll(ibPts, i)
    {
        const label surfCell = (i < surfCells.size() ? surfCells[i] : -1);
        os << surfCell << "\n";
    }

    os  << "\nSCALARS processor int 1\n"
        << "LOOKUP_TABLE default\n";

    forAll(ibPts, i)
    {
        os << Pstream::myProcNo() << "\n";
    }
}
//--------------------------------------------------------------------------//
void openHFDIBDEM::writeBodiesInfo()
{
    if(!recordSimulation_)
        return;

    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());


    mkDir(curOutDir);
    mkDir(curOutDir +"/stlFiles");
    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }
    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateDEM(volScalarField& body,volScalarField& refineF)
{
    scalar deltaTime(mesh_.time().deltaT().value());
    scalar pos(0.0);
    scalar step(stepDEM_);
    // scalar timeStep(step*deltaTime);
    List<DynamicList<pointField>> bodiesPositionList(Pstream::nProcs());
    // Infos <<bodiesPositionList.size() << endl;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> syncOutForceKeyTable;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> contactResolvedKeyTable;
    HashTable <label,label,Hash<label>> wallContactIBTable;
    while( pos < 1)
    {
        bodiesPositionList[Pstream::myProcNo()].clear();

        InfoH << DEM_Info << " Start DEM pos: " << pos
            << " DEM step: " << step << endl;

        InfoH << basic_Info << " DEM - CFD Time: "
            << mesh_.time().value() + deltaTime*pos << endl;

        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);

            if(Pstream::myProcNo() == 0 )
            {
	      //immersedBodies_[ib].moveImmersedBody(deltaTime*step);

		const scalar tSubEnd = mesh_.time().value() + deltaTime*(pos + step);
		immersedBodies_[ib].moveImmersedBody(deltaTime*step, tSubEnd);
		
                bodiesPositionList[Pstream::myProcNo()].append
                (
                    immersedBodies_[ib].getGeomModel().getBodyPoints()
                );
            }
        }

        Pstream::gatherList(bodiesPositionList,0);
        Pstream::broadcastList(bodiesPositionList,0);

        label bodyIndex(0);
        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].getGeomModel().setBodyPosition
            (
                bodiesPositionList[0][bodyIndex++]
            );
        }

        bodiesPositionList[Pstream::myProcNo()].clear();

        forAll (immersedBodies_,ib)
        {
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);
            immersedBodies_[ib].printBodyInfo();
        }
        pos += step;

        if (pos + step + SMALL >= 1)
        {
            step = 1 - pos;
        }
    }
}

//---------------------------------------------------------------------------//
// function to either add or remove bodies from the simulation
void openHFDIBDEM::addRemoveBodies
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF
)
{
    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);

        label maxAdditions(50);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));

            cAddition++;

            if (addModels_[modelI].getBodyAdded())
            {
                InfoH << addModel_Info << "STL file correctly generated, registering the new body" << endl;

                // prepare pointer list for IBs (increase its size)
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                // create the new body
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );

                // get reference for further processing
                immersedBody& nBody(immersedBodies_[addIBPos]);
                nBody.createImmersedBody(body,refineF);
                nBody.computeBodyCharPars();
                if (nBody.getStartSynced())
                {
                    nBody.initSyncWithFlow(U);
                }
                // verletList_.addBodyToVList(nBody);

                InfoH << addModel_Info
                    << "new body included into the simulation" << endl;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info
                    << "new body should have been added but was not "
                    << "(probably overlap with an existing body)"
                    << endl;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateFSCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].pimpleUpdate(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::restartSimulation
(
    volScalarField& body,
    volScalarField& refineF,
    word runTime
)
{
    word timePath(recordOutDir_+"/"+runTime);
    fileNameList files(readDir(timePath));

    forAll(files,f)
    {
        IOdictionary bodyDict
        (
            IOobject
            (
                timePath + "/" + files[f],
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word bodyId(std::to_string(readLabel(bodyDict.lookup("bodyId"))));
        word bodyName(bodyDict.lookup("bodyName"));
        vector Vel(bodyDict.lookup("Vel"));
        scalar omega(readScalar(bodyDict.lookup("omega")));
        vector Axis(bodyDict.lookup("Axis"));
        bool isStatic(readBool(bodyDict.lookup("static")));
        label timeStepsInContWStatic(readLabel(bodyDict.lookup("timeStepsInContWStatic")));

        std::shared_ptr<geomModel> bodyGeomModel;
        word bodyGeom;
        if (HFDIBDEMDict_.subDict(bodyName).found("bodyGeom"))
        {
            const word input(HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom"));
            bodyGeom = input;
            InfoH << iB_Info << "Found bodyGeom for "
                << bodyName << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            notImplemented("STOP");
        }

        if(bodyGeom == "convex")
        {
            // Convex-body restart support is not implemented here.
        }
        else if(bodyGeom == "nonConvex")
        {
            notImplemented("STOP");
        }
        else if(bodyGeom == "sphere")
        {
            notImplemented("STOP");
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            notImplemented("STOP");
        }

        label newIBSize(immersedBodies_.size()+1);
        label addIBPos(newIBSize - 1);
        immersedBodies_.setSize(newIBSize);

        InfoH << iB_Info << "Restarting body: " << bodyId << " as "
            << addIBPos << " bodyName: " << bodyName << endl;
        immersedBodies_.set
        (
            addIBPos,
            new immersedBody
            (
                bodyName,
                mesh_,
                HFDIBDEMDict_,
                transportProperties_,
                addIBPos,
                recomputeM0_,
                bodyGeomModel,
                ibInterp_,
                cellPoints_
            )
        );

        immersedBodies_[addIBPos].createImmersedBody(body,refineF);
        immersedBodies_[addIBPos].computeBodyCharPars();
        immersedBodies_[addIBPos].setRestartSim(Vel,omega,Axis,isStatic,timeStepsInContWStatic);
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preCalculateCellPoints()
{
    cellPoints_.clear();
    cellPoints_.setSize(mesh_.nCells());
    forAll(mesh_.C(), cellI)
    {
        cellPoints_[cellI] = mesh_.cellPoints()[cellI];
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeFirtsTimeBodiesInfo()
{
    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());
    bool checkExistance(false);
    if(!recordSimulation_ || isDir(curOutDir))
        return;
    reduce(checkExistance,orOp<bool>());
    if(Pstream::myProcNo() == 0)
    {
        mkDir(curOutDir);
        mkDir(curOutDir +"/stlFiles");
    }
    reduce(checkExistance,orOp<bool>());

    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }

    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
// void openHFDIBDEM::setSolverInfo()
// {
//     solverInfo::setOnlyDEM(true);
// }
//---------------------------------------------------------------------------//
