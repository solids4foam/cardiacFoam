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
Description
    class for immersed bodies representation.
SourceFiles
    immersedBodies.C
Contributors
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "immersedBody.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"
#include "meshSearch.H"
#include "List.H"
#include "ListOps.H"

#include "OFstream.H"

#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"
#include "PstreamReduceOps.H"

#include "fvcSmooth.H"
#include "fvMeshSubset.H"
#include "processorPolyPatch.H"

#define ORDER 2

using namespace Foam;

//---------------------------------------------------------------------------//
immersedBody::immersedBody
(
    word bodyName,
    const Foam::fvMesh& mesh,
    dictionary& HFDIBDEMDict,
    dictionary& transportProperties,
    label bodyId,
    label recomputeM0,
    std::shared_ptr<geomModel> bodyGeomModel,
    autoPtr<ibInterpolation>& ibIntp,
    List<labelList>& cellPoints
)
:
bodyName_(bodyName),
isActive_(true),
immersedDict_(HFDIBDEMDict.subDict(bodyName_)),
mesh_(mesh),
transportProperties_(transportProperties),
geomModel_(std::move(bodyGeomModel)),
cellPoints_(cellPoints),
Axis_(vector::one),
AxisOld_(vector::one),
omega_(0.0),
omegaOld_(0.0),
Vel_(vector::zero),
VelOld_(vector::zero),
a_(vector::zero),
alpha_(vector::zero),
totalAngle_(vector::zero),
CoNum_(0.0),
rhoF_("rho", transportProperties_),
bodyId_(bodyId),
updateTorque_(false),
bodyOperation_(0),
octreeField_(mesh_.nCells(), 0),
cellToStartInCreateIB_(0),
startSynced_(false),
totRotMatrix_(tensor::I),
sdBasedLambda_(false),
intSpan_(2.0),
charCellSize_(1e3),
refineBuffers_(0),
recomputeM0_(recomputeM0),
timesToSetStatic_(-1),
staticContactPost_(vector::zero)
{
    #include "initializeIB.H"

  {
    // If the underlying geom model is STL-based and a valve dict is present, parse it once.
    stlBased* sb = dynamic_cast<stlBased*>(geomModel_.get());

    //     if (sb && immersedDict_.found("deformingBody_mitral_ref"))
    // {
    //     const dictionary& d = immersedDict_.subDict("deformingBody_mitral_ref");
    //     sb->readValveSliceAxisRef(d);
    //     InfoH << basic_Info
    //           << "Valve motion enabled (deformingBody_mitral_ref) for body "
    //           << bodyName_ << endl;
    // }
    if (sb && immersedDict_.found("deformingBody_mitral_ref"))
    {
        const dictionary& d = immersedDict_.subDict("deformingBody_mitral_ref");
        sb->readValveSliceAxisRef(d);

        word mitralRefVelocityMode
        (
            d.lookupOrDefault<word>("mitralRefVelocityMode", word("off"))
        );

        Foam::string modeStr = mitralRefVelocityMode;
        for (auto& ch : modeStr)
        {
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        }
        mitralRefVelocityMode = word(modeStr);

        mitralRefVelocityOff_ = true;
        mitralRefVelocityIbPointsOnly_ = false;

        if (mitralRefVelocityMode == "ibpointsonly")
        {
            mitralRefVelocityOff_ = false;
            mitralRefVelocityIbPointsOnly_ = true;
        }
        else if (mitralRefVelocityMode == "off")
        {
            mitralRefVelocityOff_ = true;
            mitralRefVelocityIbPointsOnly_ = false;
        }
        else
        {
            WarningInFunction
                << "Unknown mitralRefVelocityMode = " << mitralRefVelocityMode
                << " for body " << bodyName_
                << ". Using 'off'." << endl;

            mitralRefVelocityOff_ = true;
            mitralRefVelocityIbPointsOnly_ = false;
        }

        InfoH << basic_Info
              << "Valve motion enabled (deformingBody_mitral_ref) for body "
              << bodyName_
              << " ; mitralRefVelocityMode = " << mitralRefVelocityMode
              << endl;
    }
     if (sb && immersedDict_.found("deformingBody_mitral"))
      {
        const dictionary& d = immersedDict_.subDict("deformingBody_mitral");
        sb->readValveSliceAxis(d);   // also sets mitralEnabled_ = true
        InfoH << basic_Info << "Valve motion enabled (deformingBody_mitral) for body "
              << bodyName_ << endl;
      }
     if (sb && immersedDict_.found("deformingBody_breathingSphere"))
        {
            const dictionary& d = immersedDict_.subDict("deformingBody_breathingSphere");
            sb->readBreathingSphere(d);
            InfoH << basic_Info
                  << "Breathing-sphere motion enabled (deformingBody_breathingSphere) for body "
                  << bodyName_ << endl;
        }
     if (sb && immersedDict_.found("deformingBody_sinusoidalTranslation"))
       {
	 const dictionary& d = immersedDict_.subDict("deformingBody_sinusoidalTranslation");
	 sb->readSinusoidalTranslation(d);
	 InfoH << basic_Info
	       << "Sinusoidal translation enabled (deformingBody_sinusoidalTranslation) for body "
	       << bodyName_ << endl;
       }
       if (sb && immersedDict_.found("deformingBody_quadraticBend"))
       {
         const dictionary& d =
             immersedDict_.subDict("deformingBody_quadraticBend");
         sb->readQuadraticBend(d);
         InfoH << basic_Info
               << "Quadratic bending enabled (deformingBody_quadraticBend) for body "
               << bodyName_ << endl;
       }

       if (sb && immersedDict_.found("deformingBody_customProfileBend"))
       {
         const dictionary& d =
             immersedDict_.subDict("deformingBody_customProfileBend");
         sb->readCustomProfileBend(d);
         InfoH << basic_Info
               << "Custom profile bending enabled (deformingBody_customProfileBend) for body "
               << bodyName_ << endl;
       }
  }

    InfoH << iB_Info << "Finished body initialization" << endl;
    InfoH << basic_Info << "New bodyID: " << bodyId_ << " name: "
        << bodyName_ << " rhoS: " << geomModel_->getRhoS()
        << " dC: " << getDC() << endl;
}
//---------------------------------------------------------------------------//
immersedBody::~immersedBody()
{
}
//---------------------------------------------------------------------------//
// Create immersed body info
void immersedBody::createImmersedBody
(
    volScalarField& body,
    volScalarField& refineF,
    bool synchCreation
)
{
    geomModel_->createImmersedBody(
        body,
        octreeField_,
        cellPoints_
    );

    if(synchCreation)
    {
        syncCreateImmersedBody(body, refineF);
    }
}
//---------------------------------------------------------------------------//
void immersedBody::syncCreateImmersedBody
(
    volScalarField& body,
    volScalarField& refineF
)
{
    geomModel_->setOwner();
    InfoH << iB_Info << "body: " << bodyId_
        << " owner: " << geomModel_->getOwner() << endl;

    InfoH << iB_Info << "Computing geometrical properties" << endl;
    geomModel_->calculateGeometricalProperties(body);

    // update body courant number
    computeBodyCoNumber();

    InfoH << iB_Info << "-- body: " << bodyId_
        << " current center of mass position: " << geomModel_->getCoM() << endl;

    const List<DynamicLabelList>& surfCells = geomModel_->getSurfaceCellList();
    DynamicLabelList zeroList(surfCells[Pstream::myProcNo()].size(), 0);

    constructRefineField(
        body,
        refineF,
        surfCells[Pstream::myProcNo()],
        zeroList
    );

    scalarList charCellSizeL(Pstream::nProcs(),1e4);
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells[Pstream::myProcNo()][sCellI];
        charCellSizeL[Pstream::myProcNo()] =
            min(charCellSizeL[Pstream::myProcNo()],
                Foam::pow(mesh_.V()[cellI],0.3333)
            );
    }
    forAll(charCellSizeL,indl)
    {
        if(charCellSizeL[indl] > 5e3)
        {
            charCellSizeL[indl] = -1.0;
        }
    }

    charCellSize_ = gMax(charCellSizeL);
    InfoH << iB_Info << "Body characteristic cell size: "
        << charCellSize_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::syncImmersedBodyParralell1
(
    volScalarField& body,
    volScalarField& refineF
)
{
    geomModel_->setOwner();
    InfoH << iB_Info << "body: " << bodyId_
        << " owner: " << geomModel_->getOwner() << endl;

    InfoH << iB_Info << "Computing geometrical properties" << endl;
    geomModel_->calculateGeometricalPropertiesParallel(body);
}
//---------------------------------------------------------------------------//
void immersedBody::syncImmersedBodyParralell2
(
    volScalarField& body,
    volScalarField& refineF
)
{
    // update body courant number
    // computeBodyCoNumber();

    InfoH << iB_Info << "-- body: " << bodyId_
        << " current center of mass position: " << geomModel_->getCoM() << endl;

    const List<DynamicLabelList>& surfCells = geomModel_->getSurfaceCellList();
    DynamicLabelList zeroList(surfCells[Pstream::myProcNo()].size(), 0);

    constructRefineField(
        body,
        refineF,
        surfCells[Pstream::myProcNo()],
        zeroList
    );

    scalarList charCellSizeL(Pstream::nProcs(),1e4);
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI = surfCells[Pstream::myProcNo()][sCellI];
        charCellSizeL[Pstream::myProcNo()] =
            min(charCellSizeL[Pstream::myProcNo()],
                Foam::pow(mesh_.V()[cellI],0.3333)
            );
    }
    forAll(charCellSizeL,indl)
    {
        if(charCellSizeL[indl] > 5e3)
        {
            charCellSizeL[indl] = -1.0;
        }
    }

    charCellSize_ = gMax(charCellSizeL);
    InfoH << iB_Info << "Body characteristic cell size: "
        << charCellSize_ << endl;
}
//---------------------------------------------------------------------------//
void immersedBody::constructRefineField
(
    volScalarField& body,
    volScalarField& refineF,
    DynamicLabelList cellsToIterate,
    DynamicLabelList startLevel
)
{
    if(refineBuffers_ == 0)
        return;

    DynamicLabelList cellsToIterateC;
    DynamicLabelList cellsToIterateF;

    List<DynamicLabelList> cellsToSendToProcs;
    cellsToSendToProcs.setSize(Pstream::nProcs());
    List<DynamicLabelList> cellsToSendToProcsLevel;
    cellsToSendToProcsLevel.setSize(Pstream::nProcs());

    for(label i = 0; i < refineBuffers_; i++)
    {
        forAll(cellsToIterate, cellI)
        {
            if(startLevel[cellI] == i)
                cellsToIterateC.append(cellsToIterate[cellI]);
        }

        forAll(cellsToIterateC, cellI)
        {
            labelList cellFaces(mesh_.cells()[cellsToIterateC[cellI]]);
            forAll(cellFaces, faceI)
            {
                if (mesh_.isInternalFace(cellFaces[faceI]))
                {
                    label nCell(mesh_.owner()[cellFaces[faceI]]);
                    if(nCell == cellsToIterateC[cellI])
                    {
                        nCell = mesh_.neighbour()[cellFaces[faceI]];
                    }

                    if(refineF[nCell] == 0)
                    {
                        if(i > 0)
                        {
                            if(body[nCell] < SMALL)
                            {
                                refineF[nCell] = 1;
                                cellsToIterateF.append(nCell);
                            }
                        }
                        else
                        {
                            refineF[nCell] = 1;
                            cellsToIterateF.append(nCell);
                        }
                    }
                }
                else
                {
                    label facePatchId(mesh_.boundaryMesh().whichPatch(
                        cellFaces[faceI]
                    ));

                    const polyPatch& cPatch = mesh_.boundaryMesh()[facePatchId];
                    if (cPatch.type() == "processor")
                    {
                        const processorPolyPatch& procPatch =
                            refCast<const processorPolyPatch>(cPatch);
                        if (procPatch.myProcNo() == Pstream::myProcNo())
                        {
                            cellsToSendToProcs[procPatch.neighbProcNo()].append(
                                cPatch.whichFace(cellFaces[faceI])
                            );
                            cellsToSendToProcsLevel[procPatch.neighbProcNo()]
                                .append(i+1);
                        }
                        else
                        {
                            cellsToSendToProcs[procPatch.myProcNo()].append(
                                cPatch.whichFace(cellFaces[faceI])
                            );
                            cellsToSendToProcsLevel[procPatch.myProcNo()]
                                .append(i+1);
                        }
                    }
                }
            }
        }
        cellsToIterateC = cellsToIterateF;
        cellsToIterateF.clear();
    }

    List<DynamicLabelList> facesReceivedFromProcs;
    List<DynamicLabelList> cellsReceivedFromProcsLevel;

    // send points that are not on this proc to other proc
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UOPstream send(proci, pBufs);
            send << cellsToSendToProcs[proci];
        }
    }
    pBufs.finishedSends();
    // recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            facesReceivedFromProcs.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            facesReceivedFromProcs.append(recList);
        }
    }
    pBufs.clear();

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {

            UOPstream send(proci, pBufs);
            send << cellsToSendToProcsLevel[proci];
        }
    }
    pBufs.finishedSends();

    // recieve points from other procs
    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        if (proci != Pstream::myProcNo())
        {
            UIPstream recv(proci, pBufs);
            DynamicLabelList recList (recv);
            cellsReceivedFromProcsLevel.append(recList);
        }
        else
        {
            DynamicLabelList recList;
            cellsReceivedFromProcsLevel.append(recList);
        }
    }
    pBufs.clear();


    DynamicLabelList newCellsToIterate;
    DynamicLabelList newCellsToIterateStartLevel;

    // check if some point from other proc is on this processor
    for (label otherProci = 0;
        otherProci < facesReceivedFromProcs.size();
        otherProci++
    )
    {
        for (label faceI = 0;
            faceI < facesReceivedFromProcs[otherProci].size();
            faceI++
        )
        {
            label cellProcI(0);
            forAll (mesh_.boundaryMesh(), patchi)
            {
                const polyPatch& cPatch = mesh_.boundaryMesh()[patchi];
                if (cPatch.type() == "processor")
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(cPatch);

                    if (procPatch.myProcNo() == Pstream::myProcNo()
                        && procPatch.neighbProcNo() == otherProci)
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start()
                            + facesReceivedFromProcs[otherProci][faceI]];
                        if(refineF[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(
                                cellsReceivedFromProcsLevel[otherProci][faceI]
                            );
                        }
                        break;
                    }
                    else if (procPatch.myProcNo() == otherProci
                        && procPatch.neighbProcNo() == Pstream::myProcNo())
                    {
                        cellProcI = mesh_.faceOwner()[cPatch.start()
                            + facesReceivedFromProcs[otherProci][faceI]];
                        if(refineF[cellProcI] == 0)
                        {
                            newCellsToIterate.append(cellProcI);
                            newCellsToIterateStartLevel.append(
                                cellsReceivedFromProcsLevel[otherProci][faceI]
                            );
                        }
                        break;
                    }
                }
            }
        }
    }

    bool contBool(false);
    if(newCellsToIterate.size() > 0)
    {
        contBool = true;
    }

    reduce(contBool, orOp<bool>());

    if(contBool)
    {
        constructRefineField
        (
            body,
            refineF,
            newCellsToIterate,
            newCellsToIterateStartLevel
        );
    }
}
//---------------------------------------------------------------------------//
void immersedBody::postPimpleUpdateImmersedBody
(
    volScalarField& body,
    volVectorField& f
)
{
    // update Vel_, Axis_ and omega_
    updateCoupling(body,f);

    Vel_ = VelOld_;
    Axis_ = AxisOld_;
    omega_ = omegaOld_;
}
//---------------------------------------------------------------------------//
void immersedBody::updateCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
    vector FV(vector::zero);
    vector TA(vector::zero);

    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists(
        intLists,
        surfLists,
        refCoMList
    );

  // calcualate viscous force and torque

    forAll (intLists, i)
    {
        DynamicLabelList& intListI = intLists[i];
        forAll (intListI, intCell)
        {
            label cellI = intListI[intCell];

            FV -=  f[cellI]*mesh_.V()[cellI];
            TA -=  ((mesh_.C()[cellI] - refCoMList[i])^f[cellI])
                *mesh_.V()[cellI];
        }
    }

    forAll (surfLists, i)
    {
        DynamicLabelList& surfListI = surfLists[i];
        forAll (surfListI, surfCell)
        {
            label cellI = surfListI[surfCell];

            FV -=  f[cellI]*mesh_.V()[cellI];
            TA -=  ((mesh_.C()[cellI] - refCoMList[i])^f[cellI])
                *mesh_.V()[cellI];
        }
    }

  reduce(FV, sumOp<vector>());
  reduce(TA, sumOp<vector>());
  FV *= rhoF_.value();
  TA *= rhoF_.value();

  FCoupling_ = forces(FV, TA);
}
//---------------------------------------------------------------------------//
// update movement variables of the body
void immersedBody::updateMovement
(
    scalar deltaT
)
{
    InfoH << iB_Info
        << "immersedBody::updateMovement: Vel_ = " << Vel_ << endl;
    updateMovementComp(deltaT,Vel_,Axis_,omega_);
}
void immersedBody::updateMovement
(
    vector Vel,
    vector Axis,
    scalar omega
)
{
    scalar deltaT = mesh_.time().deltaT().value();
    updateMovementComp(deltaT,Vel,Axis,omega);
}
void immersedBody::updateMovementComp
(
    scalar deltaT,
    vector Vel,
    vector Axis,
    scalar omega
)
{
    auto updateTranslation = [&]()
    {

        const uniformDimensionedVectorField& g =
            mesh_.lookupObject<uniformDimensionedVectorField>("g");

        vector FG(vector::zero);
        // if(!solverInfo::getOnlyDEM())
        //     FG = geomModel_->getM0()*(1.0-rhoF_.value()
        //     /geomModel_->getRhoS().value())*g.value();
        // else
            FG = geomModel_->getM0()*g.value();

        vector F(FCoupling_.F);
        // F += FContact_.F;
        F += FG;

        if(!case3D)
        {
            F[emptyDim] *= 0;
            FG[emptyDim] *= 0;
        }
        if(geomModel_->getM0() > 0)
        {
            // compute current acceleration (assume constant over timeStep)

            InfoH << iB_Info <<"-- body "<< bodyId_ <<" ParticelMass  : " << geomModel_->getM0() << endl;
            InfoH << iB_Info <<"-- body "<< bodyId_ <<" Acting Force  : " << F << endl;
            a_  = F/(geomModel_->getM0());
            // update body linear velocity
            Vel_ = Vel + deltaT*a_;
            InfoH << iB_Info <<"-- body "<< bodyId_ <<" accelaration  : " << a_ << endl;
        }
    };

    auto updateRotation = [&]()
    {
        if(mag(geomModel_->getI()) > 0)
        {
            vector T(FCoupling_.T);
            // T += FContact_.T;

            // update body angular acceleration
            alpha_ = inv(geomModel_->getI()) & T;
            // update body angular velocity
            vector Omega(Axis*omega + deltaT*alpha_);
            // split Omega into Axis_ and omega_
            omega_ = mag(Omega);

            if (omega_ < SMALL)
            {
                Axis_ = vector::one;
                if (!case3D)
                {
                    const vector validDirs = (geometricD + vector::one)/2;
                    Axis_ -= validDirs;
                }
            }
            else
            {
                Axis_ =  Omega/(omega_+SMALL);
                if (!case3D)
                {// in 2D, I need to keep only the part of the rotation axis
                    const vector validDirs = (geometricD + vector::one)/2;
                    Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
                }
            }
            Axis_ /= mag(Axis_);
        }
    };

    auto updateRotationFixedAxis = [&]()
    {
        vector T(FCoupling_.T);
        // T += FContact_.T;

        // update body angular velocity
        vector Omega(Axis*omega + deltaT * (inv(geomModel_->getI()) & T));

        // split Omega into Axis_ and omega_
        omega_ = mag(Omega);

        vector newAxis = Omega/(omega_+SMALL);
        if ((newAxis & Axis_) < 0) Axis_ *= (-1.0);;
    };

    if (bodyOperation_ == 0 or bodyOperation_ == 3)
    {
        return;
    }
    else if (bodyOperation_ == 1)
    {
        updateRotation();
        return;
    }
    else if (bodyOperation_ == 2)
    {
        updateTranslation();
        return;
    }
    else if (bodyOperation_ == 4)
    {
        updateRotationFixedAxis();
        return;
    }

    updateTranslation();
    if (updateTorque_)
    {
        updateRotation();
    }

    return;
}
//---------------------------------------------------------------------------//
// move immersed body according to body operation
//void immersedBody::moveImmersedBody(scalar deltaT)

//-----------------06 March 2026, changes---------------------------
void immersedBody::moveImmersedBody(scalar deltaT, scalar motionTime)
{
    // NOTE: This is called ONLY on proc 0 from updateDEM.
    // All parallel sync is done via bodiesPositionList + setBodyPosition.
    // So: NO Pstream, NO synchronPos() in here.

  //------------- 06 march 2026, new changes----------------

  scalar tEval = motionTime;
  if (tEval <= -0.5*GREAT)
    {
      tEval = mesh_.time().value();
    }

  lastMotionTime_ = tEval;
  //--------------------------------------------

    // 1) Apply STL-based intrinsic deformation (valves, breathing)
    if (stlBased* sb = dynamic_cast<stlBased*>(&getGeomModel()))
    {
        if (sb->mitralRefEnabled())
        {
            sb->applyValveSliceAxisRef(tEval);
            InfoH << iB_Info
                << "Applied slice-axis REF valve motion for body "
                << bodyName_ << endl;
        }
        if (sb->mitralEnabled())
        {
	  sb->applyValveSliceAxis(mesh_.time());
	  //sb->applyValveSliceAxis(tEval);
            InfoH << iB_Info
                << "Applied slice-axis valve motion for body "
                << bodyName_ << endl;
        }

        if (sb->breathingEnabled())
        {
            sb->applyBreathingSphere(mesh_.time());
            InfoH << iB_Info
                << "Applied breathing-sphere motion for body "
                << bodyName_ << endl;
        }
	// if (sb->sinusoidalEnabled())
	// {
	//     sb->applySinusoidalTranslation(mesh_.time());
	//     Info << "Applied sinusoidal translation for body " << bodyName_ << endl;
	// }

	//---------------- changes 06 march 2026---------------

	if (sb->sinusoidalEnabled())
	  {
	    sb->applySinusoidalTranslation(tEval);

	    const vector vSin = sb->sinusoidalVelocity(tEval);

	    InfoH << iB_Info
		 << "Applied sinusoidal translation for body "
		 << bodyName_
		 << " at t = " << tEval
		 << " ; Vel_ = " << Vel_
		 << " ; analyticalSinV = " << vSin
		 << " ; |analyticalSinV| = " << mag(vSin)
		 << endl;
	  }

	//------------------------------------------------

	//------------- change 10 march 2026-----------
        if (sb->quadraticBendEnabled())
        {
            sb->applyQuadraticBend(tEval);

	       InfoH << iB_Info
                 << "Applied quadratic bend for body "
                 << bodyName_
                 << " at t = " << tEval
                 << endl;
        }
	//  if (sb->quadraticBendEnabled())
        //{
        //    sb->applyQuadraticBend(mesh_.time());

	    //  Info << "Applied quadratic bend for body "
	    //    << bodyName_
	    //   << " at CFD time = " << mesh_.time().value()
	//   << endl;
        //}
	//-----------------------------------------------


	if (sb->customProfileBendEnabled())
        {
            sb->applyCustomProfileBend(tEval);

            scalar maxUx = 0.0;
            scalar maxUy = 0.0;
            scalar maxUmag = 0.0;

            sb->customProfileBendMaxVelocity(tEval, maxUx, maxUy, maxUmag);

            InfoH << iB_Info
                << "Applied custom profile bend for body "
                << bodyName_
                << " at t = " << tEval
                << " ; max|Ux| = " << maxUx
                << " ; max|Uy| = " << maxUy
                << " ; max|U| = " << maxUmag
                << endl;
        }
    }

    // 2) Static body: can still deform (we just did), but no rigid motion.
    if (bodyOperation_ == 0)
    {
        return;
    }

    // 3) Rigid-body motion (unchanged)
    if (mag(deltaT + 1.0) < SMALL)
    {
        deltaT = mesh_.time().deltaT().value();
    }

    scalar angle     = omega_*deltaT;
    vector transIncr = Vel_*deltaT;

    tensor rotMatrix(Foam::cos(angle)*tensor::I);
    rotMatrix += Foam::sin(angle)*tensor(
        0.0,        -Axis_.z(),  Axis_.y(),
        Axis_.z(),   0.0,        -Axis_.x(),
        -Axis_.y(),  Axis_.x(),   0.0
    );
    rotMatrix += (1.0 - Foam::cos(angle))*(Axis_ * Axis_);

    totRotMatrix_ = rotMatrix & totRotMatrix_;

    geomModel_->bodyRotatePoints(angle, Axis_);
    geomModel_->bodyMovePoints(transIncr);

    // 4) DO NOT CALL synchronPos() HERE.
    // updateDEM will:
    //   - read getBodyPoints() on proc 0
    //   - broadcast via gather/scatter
    //   - setBodyPosition(...) on all ranks.
}
//---------------------------------------------------------------------------//
void immersedBody::printBodyInfo()
{
    InfoH << iB_Info;
    InfoH << "-- body " << bodyId_ << " CoM                  : "
        << geomModel_->getCoM() << endl;
    InfoH << "-- body " << bodyId_ << " linear velocity      : "
        << Vel_ << endl;
    InfoH << "-- body " << bodyId_ << " angluar velocity     : "
        << omega_ << endl;
    InfoH << "-- body " << bodyId_ << " axis of rotation     : "
        << Axis_ << endl;
    InfoH << "-- body " << bodyId_ << " total rotation matrix: "
        << totRotMatrix_ << endl;
}
void immersedBody::updateVectorField
(
    volVectorField& VS,
    word VName,
    volScalarField& body
)
{
    word BC(immersedDict_.subDict(VName).lookup("BC"));

    List<DynamicLabelList> intLists;
    List<DynamicLabelList> surfLists;
    DynamicVectorList refCoMList;

    geomModel_->getReferencedLists
    (
        intLists,
        surfLists,
        refCoMList
    );

    if (BC == "noSlip")
    {
        stlBased* sb = dynamic_cast<stlBased*>(geomModel_.get());

        const bool useAnalyticalCustomV =
            (sb && sb->hasAnalyticalCustomProfileBendVelocity());

        const bool useAnalyticalQuadV =
            (sb && sb->hasAnalyticalQuadraticBendVelocity());


        const bool useAnalyticalMitralV =
            (sb && sb->hasAnalyticalMitralVelocity());

        const bool useAnalyticalMitralRefV =
            (
                sb
             && sb->hasAnalyticalMitralRefVelocity()
             && !mitralRefVelocityOff_
             && !mitralRefVelocityIbPointsOnly_
            );	

        const bool useAnalyticalSinV =
            (sb && sb->hasAnalyticalSinusoidalVelocity());

        const scalar tNow =
            (lastMotionTime_ > -0.5*GREAT ? lastMotionTime_ : mesh_.time().value());

        const vector analyticalSinV =
            (useAnalyticalSinV ? sb->sinusoidalVelocity(tNow) : vector::zero);

        scalar localMaxVSmag = -GREAT;
        vector localMaxVSvec = vector::zero;
        label  localMaxVScell = -1;
        point  localMaxVSpt = point::zero;
        label  localOwnerRank = Pstream::myProcNo();

        InfoH << iB_Info
             << "updateVectorField body " << bodyName_
             << " t = " << tNow
             << " useAnalyticalCustomV = " << useAnalyticalCustomV
             << " useAnalyticalQuadV = " << useAnalyticalQuadV
             << " useAnalyticalMitralV = " << useAnalyticalMitralV
             << " useAnalyticalMitralRefV = " << useAnalyticalMitralRefV
             << " useAnalyticalSinV = " << useAnalyticalSinV
             << endl;

        auto updateLocalMax =
        [&](const label cellI, const vector& VSvalue)
        {
            const scalar vmag = mag(VSvalue);
            if (vmag > localMaxVSmag)
            {
                localMaxVSmag = vmag;
                localMaxVSvec = VSvalue;
                localMaxVScell = cellI;
                localMaxVSpt = mesh_.C()[cellI];
                localOwnerRank = Pstream::myProcNo();
            }
        };

        if (bodyOperation_ == 0)
        {
            forAll(intLists, i)
            {
                DynamicLabelList& intListI = intLists[i];
                forAll(intListI, intCell)
                {
                    const label cellI = intListI[intCell];

                    if (useAnalyticalCustomV)
                    {
                        VS[cellI] =
                            sb->customProfileBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalQuadV)
                    {
                        VS[cellI] =
                            sb->quadraticBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralRefV)
                    {
                        VS[cellI] =
                            sb->valveSliceAxisRefVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralV)
                    {
                        VS[cellI] =
                            sb->valveSliceAxisVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalSinV)
                    {
                        VS[cellI] = analyticalSinV;
                    }
                    else
                    {
                        VS[cellI] = Vel_;
                    }

                    updateLocalMax(cellI, VS[cellI]);
                }
            }

            forAll(surfLists, i)
            {
                DynamicLabelList& surfListI = surfLists[i];
                forAll(surfListI, surfCell)
                {
                    const label cellI = surfListI[surfCell];

                    if (useAnalyticalCustomV)
                    {
                        VS[cellI] =
                            sb->customProfileBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalQuadV)
                    {
                        VS[cellI] =
                            sb->quadraticBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralRefV)
                    {
                        VS[cellI] =
                            sb->valveSliceAxisRefVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralV)
                    {
                        VS[cellI] =
                            sb->valveSliceAxisVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalSinV)
                    {
                        VS[cellI] = analyticalSinV;
                    }
                    else
                    {
                        VS[cellI] = Vel_;
                    }

                    updateLocalMax(cellI, VS[cellI]);
                }
            }
        }
        else
        {
            forAll(intLists, i)
            {
                DynamicLabelList& intListI = intLists[i];
                forAll(intListI, intCell)
                {
                    const label cellI = intListI[intCell];

                    vector planarVec =
                        mesh_.C()[cellI] - refCoMList[i]
                      - Axis_*((mesh_.C()[cellI] - refCoMList[i]) & Axis_);

                    vector VSvalue = vector::zero;

                    if (useAnalyticalCustomV)
                    {
                        VSvalue =
                            sb->customProfileBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalQuadV)
                    {
                        VSvalue =
                            sb->quadraticBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralRefV)
                    {
                        VSvalue =
                            sb->valveSliceAxisRefVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralV)
                    {
                        VSvalue =
                            sb->valveSliceAxisVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalSinV)
                    {
                        VSvalue = analyticalSinV;
                    }
                    else
                    {
                        VSvalue = (-(planarVec ^ Axis_)*omega_ + Vel_);
                    }

                    VS[cellI] = VSvalue;
                    updateLocalMax(cellI, VSvalue);
                }
            }

            forAll(surfLists, i)
            {
                DynamicLabelList& surfListI = surfLists[i];
                forAll(surfListI, surfCell)
                {
                    const label cellI = surfListI[surfCell];

                    vector planarVec =
                        mesh_.C()[cellI] - refCoMList[i]
                      - Axis_*((mesh_.C()[cellI] - refCoMList[i]) & Axis_);

                    vector VSvalue = vector::zero;

                    if (useAnalyticalCustomV)
                    {
                        VSvalue =
                            sb->customProfileBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalQuadV)
                    {
                        VSvalue =
                            sb->quadraticBendVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralRefV)
                    {
                        VSvalue =
                            sb->valveSliceAxisRefVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalMitralV)
                    {
                        VSvalue =
                            sb->valveSliceAxisVelocityAtPoint(mesh_.C()[cellI], tNow);
                    }
                    else if (useAnalyticalSinV)
                    {
                        VSvalue = analyticalSinV;
                    }
                    else
                    {
                        VSvalue = (-(planarVec ^ Axis_)*omega_ + Vel_);
                    }

                    VS[cellI] = VSvalue;
                    updateLocalMax(cellI, VSvalue);
                }
            }
        }

        // -------- global reduction --------
        const scalar globalMaxVSmag = returnReduce(localMaxVSmag, maxOp<scalar>());

        if (mag(globalMaxVSmag - localMaxVSmag) < SMALL && localMaxVScell >= 0)
        {
            InfoH << iB_Info
                 << "GLOBAL-MAX-CANDIDATE body " << bodyName_
                 << " rank " << localOwnerRank
                 << " t = " << tNow
                 << " maxVSvec = " << localMaxVSvec
                 << " max|VS| = " << localMaxVSmag
                 << " at cell " << localMaxVScell
                 << " point " << localMaxVSpt
                 << endl;
        }

        if (Pstream::master())
        {
            InfoH << iB_Info
                 << "GLOBAL MAX |VS| body " << bodyName_
                 << " at t = " << tNow
                 << " is " << globalMaxVSmag
                 << endl;
        }
    }
}


//---------------------------------------------------------------------------//
// Evaluates the target velocity at each IB interpolation point.  STL-based
// deformation laws take precedence; rigid-body velocity is the fallback.
vectorField immersedBody::getUatIbPoints()
{
    if (!intpInfo_.valid())
    {
        FatalErrorInFunction
            << "intpInfo_ is null for bodyId=" << bodyId_
            << " bodyName=" << bodyName_
            << abort(FatalError);
    }

    const List<point>& ibPoints = intpInfo_->getIbPoints();

    vectorField ibPointsVal(ibPoints.size(), vector::zero);

    stlBased* sb = dynamic_cast<stlBased*>(geomModel_.get());

    const bool useAnalyticalCustomV =
        (sb && sb->hasAnalyticalCustomProfileBendVelocity());

    const bool useAnalyticalQuadV =
        (sb && sb->hasAnalyticalQuadraticBendVelocity());

    const bool useAnalyticalMitralV =
        (sb && sb->hasAnalyticalMitralVelocity());

    const bool useAnalyticalMitralRefV =
        (
            sb
         && sb->hasAnalyticalMitralRefVelocity()
         && !mitralRefVelocityOff_
        );

    const bool useAnalyticalSinV =
        (sb && sb->hasAnalyticalSinusoidalVelocity());

    const scalar tNow =
    (lastMotionTime_ > -0.5*GREAT ? lastMotionTime_ : mesh_.time().value());
    
    const vector analyticalSinV =
        (useAnalyticalSinV ? sb->sinusoidalVelocity(tNow) : vector::zero);

    forAll(ibPoints, pointI)
    {
        if (useAnalyticalCustomV)
        {
            ibPointsVal[pointI] =
                sb->customProfileBendVelocityAtPoint(ibPoints[pointI], tNow);
        }
        else if (useAnalyticalQuadV)
        {
            ibPointsVal[pointI] =
                sb->quadraticBendVelocityAtPoint(ibPoints[pointI], tNow);
        }
        else if (useAnalyticalMitralRefV)
        {
            ibPointsVal[pointI] =
                sb->valveSliceAxisRefVelocityAtPoint(ibPoints[pointI], tNow);
        }
        else if (useAnalyticalMitralV)
        {
            ibPointsVal[pointI] =
                sb->valveSliceAxisVelocityAtPoint(ibPoints[pointI], tNow);
        }
        else if (useAnalyticalSinV)
        {
            ibPointsVal[pointI] = analyticalSinV;
        }
        else
        {
            vector planarVec =
                ibPoints[pointI] - geomModel_->getCoM()
              - Axis_*((ibPoints[pointI] - geomModel_->getCoM()) & Axis_);

            ibPointsVal[pointI] = (-(planarVec^Axis_)*omega_ + Vel_);
        }
    }

    scalar localMaxMag = -GREAT;
    label  localNZ     = 0;

    forAll(ibPointsVal, pointI)
    {
        const scalar magUi = mag(ibPointsVal[pointI]);

        if (magUi > localMaxMag)
        {
            localMaxMag = magUi;
        }

        if (magUi > SMALL)
        {
            localNZ++;
        }
    }

    const scalar globalMaxMag = returnReduce(localMaxMag, maxOp<scalar>());
    const label  globalNZ     = returnReduce(localNZ, sumOp<label>());

    if (Pstream::master())
    {
        InfoH << iB_Info
             << "IBVEL-SUMMARY body " << bodyName_
             << " t = " << tNow
             << " nIbPts = " << ibPoints.size()
             << " nonZeroIbPts = " << globalNZ
             << " globalMax|U| = " << globalMaxMag
             << endl;
    }

    return ibPointsVal;
}

//---------------------------
//---------------------------------------------------------------------------//
// function to move the body after the contact
void immersedBody::postContactUpdateBodyField
(
    volScalarField& body,
    volScalarField& refineF
)
{
    geomModel_->resetBody(body);

    createImmersedBody(body,refineF,false);
}
//---------------------------------------------------------------------------//
void immersedBody::recreateBodyField
(
    volScalarField& body,
    volScalarField& refineF
)
{
    octreeField_ = Field<label>(mesh_.nCells(), 0);
    geomModel_->getSurfaceCellList()[Pstream::myProcNo()].clear();
    geomModel_->getInternalCellList()[Pstream::myProcNo()].clear();
    createImmersedBody(body,refineF,false);

    label nSurf = 0;
    label nInt  = 0;

    forAll(geomModel_->getSurfaceCellList(), procI)
    {
        nSurf += geomModel_->getSurfaceCellList()[procI].size();
    }

    forAll(geomModel_->getInternalCellList(), procI)
    {
        nInt += geomModel_->getInternalCellList()[procI].size();
    }

    reduce(nSurf, sumOp<label>());
    reduce(nInt,  sumOp<label>());

    Info << "Body " << bodyName_
         << " : surface cells = " << nSurf
         << " , internal cells = " << nInt
         << endl;


}
//---------------------------------------------------------------------------//
// function to compute maximal and mean courant number of the body
void immersedBody::computeBodyCoNumber()
{
    label auxCntr(0);
    scalar VelMag(mag(Vel_));

    meanCoNum_ = 0.0;

    // rotation body courant number
    scalar rotCoNumB(omega_*getDC()*0.5*mesh_.time().deltaT().value());

    const List<DynamicLabelList>& surfCells = geomModel_->getSurfaceCellList();
    forAll (surfCells[Pstream::myProcNo()],sCellI)
    {
        label cellI(surfCells[Pstream::myProcNo()][sCellI]);

        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);

        CoNumCell+=rotCoNumB/dCell;

        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }
    const List<DynamicLabelList>& intCells = geomModel_->getInternalCellList();
    forAll (intCells[Pstream::myProcNo()],iCellI)
    {
        label cellI(intCells[Pstream::myProcNo()][iCellI]);

        scalar dCell(Foam::pow(mesh_.V()[cellI],0.3333));
        scalar CoNumCell(VelMag*mesh_.time().deltaT().value()/dCell);

        CoNum_      = max(CoNum_,CoNumCell);
        meanCoNum_ += CoNumCell;
        auxCntr    += 1;
    }

    reduce(meanCoNum_, sumOp<scalar>());
    reduce(auxCntr, sumOp<scalar>());
    reduce(CoNum_, maxOp<scalar>());

    if(auxCntr > 0)
    {
        meanCoNum_ /= auxCntr;
    }

    InfoH << iB_Info << "-- body " << bodyId_
        << " Courant Number mean: " << meanCoNum_
        << " max: " << CoNum_ << endl;

}

//---------------------------------------------------------------------------//
// print out body linear and angular momentum
void immersedBody::printMomentum()
{
    vector L(geomModel_->getI()&(Axis_*omega_));
    vector p(geomModel_->getM()*Vel_);

    InfoH << iB_Info;
    InfoH << "-- body " << bodyId_ << "  linear momentum:" << p
         << " magnitude: " << mag(p) <<endl;
    InfoH << "-- body " << bodyId_ << " angular momentum:" << L
         << " magnitude: " << mag(L) <<endl;
}
//---------------------------------------------------------------------------//
// print out body statistics
void immersedBody::printStats()
{
    vector L(geomModel_->getI()&(Axis_*omega_));
    vector p(geomModel_->getM()*Vel_);

    InfoH << iB_Info << "-- body " << bodyId_ << "  linear momentum:" << p
        << " magnitude: " << mag(p) <<endl;
    InfoH << "-- body " << bodyId_ << " angular momentum:" << L
        << " magnitude: " << mag(L) <<endl;
    InfoH << basic_Info << "-- body " << bodyId_ << " CoM :"
        << geomModel_->getCoM() << endl;
    InfoH << basic_Info << "-- body " << bodyId_ << "  linear velocity:"
        << Vel_ << " magnitude: " << mag(Vel_) <<endl;
    InfoH << "-- body " << bodyId_ << " angular velocity:" << omega_
        << " magnitude: " << mag(omega_) <<endl;
    InfoH << "-- body " << bodyId_ << "    rotation axis:" << Axis_
        << " magnitude: " << mag(Axis_) <<endl;
}
//---------------------------------------------------------------------------//
// switch the particle off (remove it from the simulation)
void immersedBody::switchActiveOff
(
    volScalarField& body
)
{
    // turn of the particle
    isActive_ = false;

    // rewrite the body field
    geomModel_->resetBody(body);
}
//---------------------------------------------------------------------------//
void immersedBody::initSyncWithFlow(const volVectorField& U)
{
    // auxiliary computation (unnecessarily expensive)
    volVectorField curlU(fvc::curl(U));
    // Note (MI): if this initialization proves OK, than this needs to
    //            be computed only ONCE for all the bodies and re-used

    // computation itself
    vector meanV(vector::zero);
    scalar totVol(0);
    vector meanC(vector::zero);
    label  cellI;
    const List<DynamicLabelList>& intCells = geomModel_->getInternalCellList();
    forAll (intCells[Pstream::myProcNo()],iCellI)
    {
        cellI   = intCells[Pstream::myProcNo()][iCellI];
        meanV  += U[cellI]*mesh_.V()[cellI];
        meanC  += curlU[cellI]*mesh_.V()[cellI];
        totVol += mesh_.V()[cellI];
    }
    reduce(meanV, sumOp<vector>());
    reduce(meanC, sumOp<vector>());
    reduce(totVol, sumOp<scalar>());
    if(totVol > 0)
    {
        Vel_ = meanV/(totVol);
        meanC/=(totVol);
    }
    vector Omega(0.5*meanC);
    if(updateTorque_)
    {
        omega_ = mag(Omega);
        if (omega_ < SMALL)
        {
            Axis_ = vector::one;
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs =
                    (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ -= validDirs;
            }
        }
        else
        {
            Axis_ =  Omega/(omega_+SMALL);
            if (mesh_.nGeometricD() < 3)
            {
                const vector validDirs =
                    (mesh_.geometricD() + Vector<label>::one)/2;
                Axis_ = cmptMultiply(vector::one-validDirs,Axis_);
            }
        }
        Axis_ /= mag(Axis_);
    }
    // update old storage
    VelOld_     = Vel_;
    omegaOld_   = omega_;
    AxisOld_    = Axis_;
    // print data:
    InfoH << basic_Info << "-- body " << bodyId_
        << "initial movement variables:" << endl;
    printStats();
}
//---------------------------------------------------------------------------//
void immersedBody::pimpleUpdate
(
    volScalarField& body,
    volVectorField& f
)
{
    updateCoupling(body, f);
    updateMovement(VelOld_, AxisOld_, omegaOld_);
}
//---------------------------------------------------------------------------//
void immersedBody::checkIfInDomain(volScalarField& body)
{
    if(geomModel_->getM0() < SMALL)
    {
        switchActiveOff(body);
        geomModel_->resetBody(body);
    }

    InfoH << iB_Info << "-- body " << bodyId_ << " current M/M0: "
        << geomModel_->getM()/geomModel_->getM0() << endl;
    // if only 1% of the initial particle mass remains in the domain, switch it off
    if (geomModel_->getM()/(geomModel_->getM0()+SMALL) < 1e-2 && case3D)
    {
        switchActiveOff(body);
        geomModel_->resetBody(body);
    }
    else if (!case3D && geomModel_->getNCells() <= 1 && !geomModel_->isCluster())
    {
        switchActiveOff(body);
        geomModel_->resetBody(body);
        InfoH << iB_Info << "-- body " << bodyId_ << " switched off" << endl;
    }
}
//---------------------------------------------------------------------------//
void immersedBody::setRestartSim(vector vel, scalar angVel, vector axisRot, bool setStatic, label timesInContact)
{
    Vel_ = vel;
    omega_ = angVel;
    Axis_ = axisRot;
    // ibContactClass_->setTimeStepsInContWStatic(timesInContact);
    // InfoH << iB_Info << "-- body " << bodyId_
    //     << " timeStepsInContWStatic_: "
    //     << ibContactClass_->getTimeStepsInContWStatic() << endl;
    if(setStatic)
    {
        bodyOperation_ = 0;
        omega_ = 0;
        Vel_ *= 0;
        InfoH << basic_Info << "-- body " << bodyId_ << " set as Static" << endl;
    }
}
//---------------------------------------------------------------------------//
void immersedBody::chceckBodyOp()
{
    if(bodyOperation_ != 5 || timesToSetStatic_ == -1)
        return;

}
