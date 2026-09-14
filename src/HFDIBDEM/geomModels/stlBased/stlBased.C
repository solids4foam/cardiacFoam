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
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "stlBased.H"
#include "vector2D.H"
#include "mathematicalConstants.H"
#include <cctype>
#include "IFstream.H"
#include "triangleFuncs.H"
#include <cmath>

using namespace Foam;

//---------------------------------------------------------------------------//
stlBased::stlBased
(
    const  fvMesh&   mesh,
    // const contactType cType,
    word      stlPath,
    scalar  thrSurf
)
:
geomModel(mesh,thrSurf),
bodySurfMesh_
(
    IOobject
    (
        stlPath,
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
stlPath_(stlPath)
{
    historyPoints_ = bodySurfMesh_.points();
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
vector stlBased::addModelReturnRandomPosition
(
    const bool allActiveCellsInMesh,
    const boundBox  cellZoneBounds,
    Random&          randGen
)
{
    vector ranVec(vector::zero);

    //meshSearch searchEng(mesh_);
    const pointField bSMeshPts(bodySurfMesh_.points());

    // get its center of mass
    vector CoM(vector::zero);
    forAll(bSMeshPts,point)
    {
        CoM += bSMeshPts[point];
    }
    CoM/= bSMeshPts.size();

    const vector validDirs = (geometricD + vector::one)/2;
    vector dirCorr(cmptMultiply((vector::one - validDirs),CoM));
    dirCorr += cmptMultiply((vector::one - validDirs),0.5*(mesh_.bounds().max() + mesh_.bounds().min()));

    boundBox bodySurfBounds(bSMeshPts);
    // compute the max scales to stay in active bounding box
    vector maxScales(cellZoneBounds.max() - bodySurfBounds.max());
    maxScales -= cellZoneBounds.min() - bodySurfBounds.min();
    maxScales *= 0.5*0.9;//0.Y is there just to be sure

    InfoH << addModel_Info << "-- addModelMessage-- "
        << "acceptable movements: " << maxScales << endl;

    scalar ranNum = 0;
    for (int i=0;i<3;i++)
    {
        ranNum = 2.0*maxScales[i]*randGen.sample01<scalar>() - 1.0*maxScales[i];
        ranVec[i] = ranNum;
    }

    ranVec = cmptMultiply(validDirs,ranVec);                            //translate only with respect to valid directions
    ranVec += dirCorr;

    return ranVec;
}
//---------------------------------------------------------------------------//
void stlBased::bodyMovePoints
(
    vector translVec
)
{
    pointField bodyPoints(bodySurfMesh_.points());
    bodyPoints += translVec;

    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stlBased::bodyScalePoints
(
    scalar scaleFac
)
{
    pointField bodyPoints(bodySurfMesh_.points());

    // get its center of mass
    vector CoM(vector::zero);
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }
    CoM/= bodyPoints.size();

    bodyPoints -= CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    bodySurfMesh_.scalePoints(scaleFac);
    bodyPoints = bodySurfMesh_.points();
    bodyPoints += CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stlBased::bodyRotatePoints
(
    scalar rotAngle,
    vector axisOfRot
)
{
    pointField bodyPoints(bodySurfMesh_.points());
    // get its center of mass
    vector CoM(vector::zero);
    forAll(bodyPoints,point)
    {
        CoM += bodyPoints[point];
    }
    CoM/= bodyPoints.size();

    tensor rotMatrix(Foam::cos(rotAngle)*tensor::I);

    rotMatrix += Foam::sin(rotAngle)*tensor(
            0.0,      -axisOfRot.z(),  axisOfRot.y(),
            axisOfRot.z(), 0.0,       -axisOfRot.x(),
        -axisOfRot.y(), axisOfRot.x(),  0.0
    );

    rotMatrix += (1.0-Foam::cos(rotAngle))*(axisOfRot * axisOfRot);

    bodyPoints -= CoM;
    bodyPoints = rotMatrix & bodyPoints;
    bodyPoints += CoM;
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
void stlBased::synchronPos(label owner)
{
    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    owner = (owner == -1) ? owner_ : owner;

    if (owner == Pstream::myProcNo())
    {
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream send(proci, pBufs);
            send << bodySurfMesh_.points();
        }
    }

    pBufs.finishedSends();
    // move body to points calculated by owner_
    UIPstream recv(owner, pBufs);
    pointField bodyPoints (recv);

    // move mesh
    bodySurfMesh_.movePoints(bodyPoints);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//--------------------------------------------------------------------------//
// --------------------- stlBased::readValveSliceAxis ------------------------
void stlBased::readValveSliceAxis_legacy(const dictionary& d)
{
    // --- Time & space parameters (from original immersedBody) ---
    period_        = d.lookupOrDefault<scalar>("period",       1.0);
    nearRamp_      = d.lookupOrDefault<scalar>("nearRamp",     0.004);
    activeLen_     = d.lookupOrDefault<scalar>("activeLen",    0.016);
    Fmax_          = d.lookupOrDefault<scalar>("Fmax",         0.995);
    timeLaw_       = d.lookupOrDefault<word>  ("timeLaw",      word("cos2"));
    leafThickness_ = d.lookupOrDefault<scalar>("leafThickness",0.0);
    spaceLaw_      = d.lookupOrDefault<word>  ("spaceLaw",     word("smootherstep"));

    // --- Fixed end (top / bottom) ---
    {
        word fixedEndWord = d.lookupOrDefault<word>("fixedEnd", word("top"));

        Foam::string fe = fixedEndWord;
        for (auto& ch : fe)
        {
            ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
        }
        fixedEndWord = word(fe);

        if (fixedEndWord == "bottom" || fixedEndWord == "cb")
        {
            mFixedEnd_ = FixedEnd::Bottom;
        }
        else
        {
            mFixedEnd_ = FixedEnd::Top;
        }
    }

    // --- Optional helpers ---
    flipNormalTowardValve_ =
        d.lookupOrDefault<Switch>("flipNormalTowardValve", false);
    valveRefFile_ =
        d.lookupOrDefault<fileName>("valveRefFile", fileName(""));

    // --- Timing windows (same semantics as old code) ---

    // defaults
    wClose1Beg_ = 0.625; wClose1End_ = 0.655;
    wOpenBeg_   = 0.940; wOpenEnd_   = 1.000;
    wClose2Beg_ = 0.0;   wClose2End_ = 0.0;

    if (timeLaw_ == "twoWindow")
    {
        if (d.found("closeWindow"))
        {
            List<scalar> cw(d.lookup("closeWindow"));
            if (cw.size() == 2) { wClose1Beg_ = cw[0]; wClose1End_ = cw[1]; }
        }
        if (d.found("openWindow"))
        {
            List<scalar> ow(d.lookup("openWindow"));
            if (ow.size() == 2) { wOpenBeg_ = ow[0]; wOpenEnd_ = ow[1]; }
        }
    }
    else if (timeLaw_ == "threeWindow")
    {
        if (d.found("close1Window"))
        {
            List<scalar> c1(d.lookup("close1Window"));
            if (c1.size() == 2) { wClose1Beg_ = c1[0]; wClose1End_ = c1[1]; }
        }
        if (d.found("openWindow"))
        {
            List<scalar> ow(d.lookup("openWindow"));
            if (ow.size() == 2) { wOpenBeg_ = ow[0]; wOpenEnd_ = ow[1]; }
        }
        if (d.found("close2Window"))
        {
            List<scalar> c2(d.lookup("close2Window"));
            if (c2.size() == 2) { wClose2Beg_ = c2[0]; wClose2End_ = c2[1]; }
        }
    }

    auto clamp01 = [](scalar v)->scalar
    {
        return (v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v));
    };

    wClose1Beg_ = clamp01(wClose1Beg_);
    wClose1End_ = clamp01(wClose1End_);
    wOpenBeg_   = clamp01(wOpenBeg_);
    wOpenEnd_   = clamp01(wOpenEnd_);
    wClose2Beg_ = clamp01(wClose2Beg_);
    wClose2End_ = clamp01(wClose2End_);

    // --- Axis: topCenter / bottomCenter, else fallback to bounds ---

    const bool haveCt = d.found("topCenter");
    const bool haveCb = d.found("bottomCenter");

    if (haveCt) mCt_ = vector(d.lookup("topCenter"));
    if (haveCb) mCb_ = vector(d.lookup("bottomCenter"));

    if (!haveCt || !haveCb)
    {
        const boundBox bb(bodySurfMesh_.points(), true);
        mCt_ = bb.max();
        mCb_ = bb.min();
    }

    // Original implementation: axis from Ct -> Cb
    vector axis = mCb_ - mCt_;
    mLen_ = mag(axis);

    if (mLen_ <= VSMALL)
    {
        FatalErrorInFunction
            << "Valve slice-axis: invalid Ct/Cb; cannot define axis for STL "
            << stlPath_ << exit(FatalError);
    }

    mAxisUnit_ = axis / mLen_;

    // --- Snapshot REST points (no drifting) ---
    historyPoints_ = bodySurfMesh_.points();

    mitralEnabled_ = true;

    InfoH << basic_Info
        << "Valve(slice-axis) enabled for STL: " << stlPath_
        << " Ct=" << mCt_
        << " Cb=" << mCb_
        << " L="  << mLen_
        << " period=" << period_
        << " timeLaw=" << timeLaw_
        << " spaceLaw=" << spaceLaw_
        << " Fmax=" << Fmax_
        << endl;
}
void stlBased::applyValveSliceAxis_legacy(const Time& runTime)
{
    if (!mitralEnabled_) return;
    if (historyPoints_.empty()) return;

    auto clamp01 = [](scalar x)->scalar
    {
        return (x < 0 ? 0 : (x > 1 ? 1 : x));
    };

    const scalar tNow = runTime.value();

    // Normalized phase φ ∈ [0,1)
    scalar phi = 0.0;
    if (period_ > VSMALL)
    {
        scalar tCycle = std::fmod(tNow, period_);
        if (tCycle < 0) tCycle += period_;
        phi = tCycle / period_;
    }

    // -------- Time gain Gt (faithful port of immersedBody::timeGain_) --------

    auto sstep = [&](scalar u)->scalar
    {
        u = clamp01(u);
        return (3.0*u*u - 2.0*u*u*u);
    };

    const scalar eps = SMALL;

    auto edgeUp = [&](scalar x, scalar a, scalar b)->scalar
    {
        if (a <= b)
        {
            if (x <= a) return 0.0;
            if (x >= b) return 1.0;
            return sstep((x - a)/max(eps, b - a));
        }
        // wrap [a,1) ∪ [0,b]
        if (x >= a) return sstep((x - a)/max(eps, 1.0 - a));
        if (x <= b) return sstep((x + (1.0 - a))/max(eps, 1.0 - a));
        return 0.0;
    };

    auto edgeDown = [&](scalar x, scalar a, scalar b)->scalar
    {
        if (a <= b)
        {
            if (x <= a) return 0.0;
            if (x >= b) return 1.0;
            return sstep((x - a)/max(eps, b - a));
        }
        // wrap
        if (x >= a) return sstep((x - a)/max(eps, 1.0 - a));
        if (x <= b) return sstep((x + (1.0 - a))/max(eps, 1.0 - a));
        return 1.0;
    };

    scalar Gt = 0.0;

    if (timeLaw_ == "twoWindow")
    {
        const scalar rise = edgeUp  (phi, wClose1Beg_, wClose1End_);
        const scalar fall = edgeDown(phi, wOpenBeg_,   wOpenEnd_);
        Gt = clamp01(rise * (1.0 - fall));
    }
    else if (timeLaw_ == "threeWindow")
    {
        const scalar r1 = edgeUp  (phi, wClose1Beg_,  wClose1End_);
        const scalar fo = edgeDown(phi, wOpenBeg_,    wOpenEnd_);
        const scalar r2 = edgeUp  (phi, wClose2Beg_,  wClose2End_);
        Gt = clamp01(r1*(1.0 - fo) + r2);
    }
    else if (timeLaw_ == "cos2")
    {
        Gt = 0.5*(1.0 - Foam::cos(2.0*Foam::constant::mathematical::pi*phi));
    }
    else if (timeLaw_ == "smoothstep")
    {
        Gt = sstep(phi);
    }
    else
    {
        const scalar s = Foam::sin(Foam::constant::mathematical::pi * phi);
        Gt = s*s;
    }

    // -------- Space gain along axis (ported from spaceGain_) --------

    auto spaceGain = [&](scalar xi)->scalar
    {
        if (xi < 0.0) return 0.0;

        if (xi <= nearRamp_)
        {
            const scalar s = clamp01(xi / max(VSMALL, nearRamp_));
            if (spaceLaw_ == "linear")     return s;
            if (spaceLaw_ == "smoothstep") return (3*s*s - 2*s*s*s);
            // smootherstep (default)
            return s*s*s*(10.0 - 15.0*s + 6.0*s*s);
        }

        if (xi <= activeLen_) return 1.0;

        return 0.0;
    };

    // -------- Rebuild from REST each step (no drift), like old code --------

    const pointField& X0 = historyPoints_;
    pointField newPts(X0.size());

    label moved = 0;
    scalar minF = GREAT, maxF = -GREAT;

    forAll(X0, i)
    {
        const vector d0 = vector(X0[i] - mCt_);

        // axial coordinate ξ from Ct along axis, clamped to [0, L]
        scalar xi = (d0 & mAxisUnit_);
        xi = (xi < 0.0 ? 0.0 : (xi > mLen_ ? mLen_ : xi));

        const vector Cxi = vector(mCt_) + xi*mAxisUnit_;

        // in-plane radial vector at rest
        vector r0 = vector(X0[i] - Cxi);
        r0 -= (r0 & mAxisUnit_) * mAxisUnit_;
        const scalar rmag0 = mag(r0);

        // distance from fixed end
        const scalar xiLocal =
            (mFixedEnd_ == FixedEnd::Bottom ? (mLen_ - xi) : xi);

        const scalar gS = spaceGain(xiLocal);
        const scalar F  = clamp01(Fmax_ * gS * Gt);

        if (rmag0 > VSMALL && F > SMALL)
        {
            const vector Xnew = Cxi + (1.0 - F) * r0;
            newPts[i] = point(Xnew);
            ++moved;
            minF = Foam::min(minF, F);
            maxF = Foam::max(maxF, F);
        }
        else
        {
            newPts[i] = X0[i]; // keep rest
        }
    }

    if (moved == 0)
    {
        minF = 0.0;
        maxF = 0.0;
    }

    if ((runTime.timeIndex() % 1) == 0)
    {
        InfoH << iB_Info
            << "valve(slice-axis): t=" << tNow
            << " phi=" << phi
            << " Gt="  << Gt
            << " movedPts=" << moved
            << " F[min,max]=[" << minF << "," << maxF << "]"
            << nl;
    }

    bodySurfMesh_.movePoints(newPts);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}

void stlBased::readValveXiField(const fileName& fName)
{
    valveXi0_.clear();
    valveXiLoaded_ = false;
    xiZeroTol_ = 0.0;
    xiMaxUsed_ = 0.0;

    if (fName.empty())
    {
        InfoH << basic_Info
            << "Valve xiFieldFile is empty; no xi field loaded." << endl;
        return;
    }

    IFstream xiFile(fName);
    if (!xiFile.good())
    {
        WarningInFunction
            << "Cannot open valve xi field file: " << fName << endl;
        return;
    }

    const pointField pts(bodySurfMesh_.points());
    const label nPts = pts.size();

    DynamicList<scalar> rawXi;
    rawXi.reserve(nPts);

    for (label i = 0; i < nPts; ++i)
    {
        scalar val;
        if (!(xiFile >> val))
        {
            FatalErrorInFunction
                << "Could not read xi value " << i
                << " from file " << fName
                << ". Expected exactly " << nPts
                << " scalar entries for STL " << stlPath_
                << exit(FatalError);
        }
        rawXi.append(val);
    }

    valveXi0_.setSize(nPts);

    scalar minVal = GREAT;
    scalar maxVal = -GREAT;
    label negCount = 0;

    forAll(rawXi, i)
    {
        minVal = min(minVal, rawXi[i]);
        maxVal = max(maxVal, rawXi[i]);
        if (rawXi[i] < 0.0) ++negCount;
    }

    const bool flipSign = (maxVal <= 0.0 && negCount > rawXi.size()/2);

    scalar minXiPos = GREAT;
    scalar maxXiPos = -GREAT;

    forAll(valveXi0_, i)
    {
        scalar x = rawXi[i];
        if (flipSign) x = -x;
        if (x < 0.0) x = 0.0;

        valveXi0_[i] = x;
        if (x > SMALL)
        {
            minXiPos = min(minXiPos, x);
            maxXiPos = max(maxXiPos, x);
        }
    }

    if (maxXiPos < SMALL)
    {
        WarningInFunction
            << "Xi field contains no positive values after processing for "
            << stlPath_ << endl;
        forAll(valveXi0_, i) valveXi0_[i] = 0.0;
        valveXiLoaded_ = true;
        return;
    }

    const scalar targetMax = (activeLen_ > SMALL ? activeLen_ : maxXiPos);
    const scalar scale = targetMax/maxXiPos;

    forAll(valveXi0_, i)
    {
        valveXi0_[i] *= scale;
    }

    scalar tolFromRamp = 0.25*nearRamp_;
    if (tolFromRamp <= SMALL) tolFromRamp = 1.0e-4;

    xiZeroTol_ = max(SMALL, min(tolFromRamp, 0.5*minXiPos));
    xiMaxUsed_ = targetMax;
    valveXiLoaded_ = true;

    InfoH << basic_Info
        << "Loaded valve xi field from " << fName
        << " for " << nPts << " STL points"
        << " raw[min,max]=[" << minVal << "," << maxVal << "]"
        << " scaledMax=" << targetMax
        << " xiZeroTol=" << xiZeroTol_
        << " flipSign=" << flipSign
        << endl;
}

void stlBased::buildAnnulusMask(const fileName& annulusFile)
{
    annulusMask_.clear();

    if (annulusFile.empty())
    {
        return;
    }

    triSurfaceMesh annSurf
    (
        IOobject
        (
            annulusFile,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const pointField annPts(annSurf.points());
    const pointField bodyPts(bodySurfMesh_.points());
    const label nPts = bodyPts.size();

    annulusMask_.setSize(nPts, false);

    if (annulusTol_ <= SMALL)
    {
        annulusTol_ = 1.0e-5;
    }

    const scalar tol2 = annulusTol_*annulusTol_;
    label nAnn = 0;

    forAll(bodyPts, i)
    {
        scalar minD2 = GREAT;
        forAll(annPts, j)
        {
            minD2 = min(minD2, magSqr(bodyPts[i] - annPts[j]));
        }

        if (minD2 < tol2)
        {
            annulusMask_[i] = true;
            ++nAnn;
        }
    }

    InfoH << basic_Info
        << "buildAnnulusMask: matched " << nAnn << " of " << nPts
        << " valve STL vertices to annulus STL " << annulusFile
        << " within tol=" << annulusTol_ << endl;

    if (nAnn == 0)
    {
        WarningInFunction
            << "No valve STL vertices matched annulus STL " << annulusFile
            << " within tol=" << annulusTol_ << endl;
    }
}

void stlBased::readValveSliceAxis_saddle(const dictionary& d)
{
    period_ = d.lookupOrDefault<scalar>("period", 1.0);
    nearRamp_ = d.lookupOrDefault<scalar>("nearRamp", 0.004);
    activeLen_ = d.lookupOrDefault<scalar>("activeLen", 0.016);
    Fmax_ = d.lookupOrDefault<scalar>("Fmax", 0.95);
    timeLaw_ = d.lookupOrDefault<word>("timeLaw", word("twoWindow"));
    leafThickness_ = d.lookupOrDefault<scalar>("leafThickness", 0.0);
    spaceLaw_ = d.lookupOrDefault<word>("spaceLaw", word("smootherstep"));

    word fixedEndWord = d.lookupOrDefault<word>("fixedEnd", word("bottom"));
    Foam::string fe = fixedEndWord;
    for (auto& ch : fe)
    {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    fixedEndWord = word(fe);
    mFixedEnd_ =
        (fixedEndWord == "bottom" || fixedEndWord == "cb")
      ? FixedEnd::Bottom
      : FixedEnd::Top;

    flipNormalTowardValve_ =
        d.lookupOrDefault<Switch>("flipNormalTowardValve", false);
    valveRefFile_ =
        d.lookupOrDefault<fileName>("valveRefFile", fileName(""));

    wClose1Beg_ = 0.625; wClose1End_ = 0.655;
    wOpenBeg_ = 0.940; wOpenEnd_ = 1.000;
    wClose2Beg_ = 0.0; wClose2End_ = 0.0;

    if (timeLaw_ == "twoWindow")
    {
        if (d.found("closeWindow"))
        {
            List<scalar> cw(d.lookup("closeWindow"));
            if (cw.size() == 2) { wClose1Beg_ = cw[0]; wClose1End_ = cw[1]; }
        }
        if (d.found("openWindow"))
        {
            List<scalar> ow(d.lookup("openWindow"));
            if (ow.size() == 2) { wOpenBeg_ = ow[0]; wOpenEnd_ = ow[1]; }
        }
    }
    else if (timeLaw_ == "threeWindow")
    {
        if (d.found("close1Window"))
        {
            List<scalar> c1(d.lookup("close1Window"));
            if (c1.size() == 2) { wClose1Beg_ = c1[0]; wClose1End_ = c1[1]; }
        }
        if (d.found("openWindow"))
        {
            List<scalar> ow(d.lookup("openWindow"));
            if (ow.size() == 2) { wOpenBeg_ = ow[0]; wOpenEnd_ = ow[1]; }
        }
        if (d.found("close2Window"))
        {
            List<scalar> c2(d.lookup("close2Window"));
            if (c2.size() == 2) { wClose2Beg_ = c2[0]; wClose2End_ = c2[1]; }
        }
    }

    wClose1Beg_ = valveClamp01_(wClose1Beg_);
    wClose1End_ = valveClamp01_(wClose1End_);
    wOpenBeg_ = valveClamp01_(wOpenBeg_);
    wOpenEnd_ = valveClamp01_(wOpenEnd_);
    wClose2Beg_ = valveClamp01_(wClose2Beg_);
    wClose2End_ = valveClamp01_(wClose2End_);

    const bool haveCt = d.found("topCenter");
    const bool haveCb = d.found("bottomCenter");

    if (haveCt) mCt_ = vector(d.lookup("topCenter"));
    if (haveCb) mCb_ = vector(d.lookup("bottomCenter"));

    if (!haveCt || !haveCb)
    {
        const boundBox bb(bodySurfMesh_.points(), true);
        mCb_ = bb.min();
        mCt_ = bb.max();
    }

    vector axis = mCt_ - mCb_;
    mLen_ = mag(axis);
    if (mLen_ <= VSMALL)
    {
        FatalErrorInFunction
            << "Valve slice-axis saddle mode: invalid Ct/Cb for STL "
            << stlPath_ << exit(FatalError);
    }
    mAxisUnit_ = axis/mLen_;

    xiFieldFile_ =
        d.lookupOrDefault<fileName>("xiFieldFile", fileName(""));
    if (!xiFieldFile_.empty())
    {
        readValveXiField(xiFieldFile_);
    }
    else
    {
        valveXiLoaded_ = false;
    }

    annulusTol_ =
        d.lookupOrDefault<scalar>("annulusTol", 0.25*leafThickness_);

    fileName annFile =
        d.lookupOrDefault<fileName>("annulusStlFile", fileName(""));

    if (!annFile.empty())
    {
        buildAnnulusMask(annFile);
    }
    else
    {
        annulusMask_.clear();
    }

    if (!valveXiLoaded_)
    {
        WarningInFunction
            << "Saddle valve mode requested, but xiFieldFile was not loaded. "
            << "Valve deformation disabled for " << stlPath_ << endl;
        mitralEnabled_ = false;
        mitralRefEnabled_ = false;
        return;
    }

    historyPoints_ = bodySurfMesh_.points();
    mitralEnabled_ = true;
    mitralRefEnabled_ = false;
    useLegacyValveAxis_ = false;

    InfoH << basic_Info
        << "Valve(saddle+xi) enabled for STL: " << stlPath_
        << " Cb=" << mCb_
        << " Ct=" << mCt_
        << " L=" << mLen_
        << " period=" << period_
        << " Fmax=" << Fmax_
        << " xiZeroTol=" << xiZeroTol_
        << " annulusMaskSize=" << annulusMask_.size()
        << endl;
}

void stlBased::applyValveSliceAxis_saddle(const Time& runTime)
{
    applyValveSliceAxis_saddle(runTime.value());
}

void stlBased::applyValveSliceAxis_saddle(const scalar tNow)
{
    if (!mitralEnabled_ && !mitralRefEnabled_) return;
    if (historyPoints_.empty() || !valveXiLoaded_) return;

    scalar phi = 0.0;
    if (period_ > VSMALL)
    {
        scalar tCycle = std::fmod(tNow, period_);
        if (tCycle < 0.0) tCycle += period_;
        phi = tCycle/period_;
    }
    const scalar Gt = valveTimeGain_(phi);

    const pointField& X0 = historyPoints_;
    const label nPts = X0.size();
    pointField newPts(nPts);

    vector axis = mCt_ - mCb_;
    const scalar Laxis = mag(axis);
    if (Laxis <= VSMALL) return;
    const vector u = axis/Laxis;

    const scalar xiAttach = min(xiZeroTol_, 0.05*nearRamp_);
    const bool useAnnulusMask = (annulusMask_.size() == nPts);

    label moved = 0;
    scalar minF = GREAT;
    scalar maxF = -GREAT;

    forAll(X0, i)
    {
        if ((useAnnulusMask && annulusMask_[i])
         || (!useAnnulusMask && valveXi0_[i] <= xiAttach))
        {
            newPts[i] = X0[i];
            continue;
        }

        vector d0 = vector(X0[i] - mCb_);
        scalar xiAxis = (d0 & u);
        xiAxis = max(scalar(0.0), min(Laxis, xiAxis));

        vector Cxi = vector(mCb_) + xiAxis*u;
        vector r0 = vector(X0[i] - Cxi);
        r0 -= (r0 & u)*u;

        const scalar rmag0 = mag(r0);
        if (rmag0 <= VSMALL)
        {
            newPts[i] = X0[i];
            continue;
        }

        const scalar xiLocal =
            (mFixedEnd_ == FixedEnd::Bottom ? xiAxis : (Laxis - xiAxis));

        const scalar gS = valveSpaceGain_(xiLocal);
        const scalar F = valveClamp01_(Fmax_*Gt*gS);

        if (F > SMALL)
        {
	  //const scalar dClose = min(F*rmag0, 0.95*rmag0);
	    const scalar dClose = F*rmag0;
            newPts[i] = point(vector(X0[i]) - dClose*(r0/rmag0));
            ++moved;
            minF = min(minF, F);
            maxF = max(maxF, F);
        }
        else
        {
            newPts[i] = X0[i];
        }
    }

    if (moved == 0)
    {
        minF = 0.0;
        maxF = 0.0;
    }

    InfoH << iB_Info
        << "valve(saddle+xi): t=" << tNow
        << " phi=" << phi
        << " Gt=" << Gt
        << " movedPts=" << moved
        << " F[min,max]=[" << minF << "," << maxF << "]"
        << " useAnnulusMask=" << useAnnulusMask
        << endl;

    bodySurfMesh_.movePoints(newPts);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}

void stlBased::readValveSliceAxis(const dictionary& d)
{
    const word mode =
        d.lookupOrDefault<word>("valveMode", word("auto"));
    const bool saddleRequested =
        d.found("xiFieldFile")
     || d.found("annulusStlFile")
     || d.lookupOrDefault<Switch>("saddleValve", false)
     || mode == "saddle"
     || mode == "saddleXi"
     || mode == "xi";

    useLegacyValveAxis_ = !saddleRequested;

    if (saddleRequested)
    {
        readValveSliceAxis_saddle(d);
    }
    else
    {
        readValveSliceAxis_legacy(d);
        mitralRefEnabled_ = false;
        useLegacyValveAxis_ = true;
    }
}

void stlBased::applyValveSliceAxis(const Time& runTime)
{
    if (useLegacyValveAxis_)
    {
        applyValveSliceAxis_legacy(runTime);
    }
    else
    {
        applyValveSliceAxis_saddle(runTime);
    }
}

//--------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// -------- Reference-based mitral slice-axis deformation/velocity -----------

void stlBased::readValveSliceAxisRef(const dictionary& d)
{
    const word mode =
        d.lookupOrDefault<word>("valveMode", word("auto"));
    const bool saddleRequested =
        d.found("xiFieldFile")
     || d.found("annulusStlFile")
     || d.lookupOrDefault<Switch>("saddleValve", false)
     || mode == "saddle"
     || mode == "saddleXi"
     || mode == "xi";

    if (saddleRequested)
    {
        readValveSliceAxis_saddle(d);
        mitralRefEnabled_ = mitralEnabled_;
        mitralEnabled_ = false;
        useLegacyValveAxis_ = false;

        InfoH << basic_Info
            << "Valve(saddle+xi, reference velocity) enabled for STL: "
            << stlPath_ << endl;
        return;
    }

    readValveSliceAxis_legacy(d);
    useLegacyValveAxis_ = true;

    //const pointField& X0 = bodySurfMesh_.points();

    const pointField X0(bodySurfMesh_.points());
	
    const label nPts = X0.size();

    valveRefXi_.setSize(nPts, 0.0);
    valveRefXiLocal_.setSize(nPts, 0.0);
    valveRefCxi_.setSize(nPts, vector::zero);
    valveRefR0_.setSize(nPts, vector::zero);

    forAll(X0, i)
    {
        const vector d0 = vector(X0[i] - mCt_);

        scalar xi = (d0 & mAxisUnit_);
        xi = max(scalar(0.0), min(xi, mLen_));

        const vector Cxi = vector(mCt_) + xi*mAxisUnit_;

        vector r0 = vector(X0[i] - Cxi);
        r0 -= (r0 & mAxisUnit_)*mAxisUnit_;

        const scalar xiLocal =
            (mFixedEnd_ == FixedEnd::Bottom ? (mLen_ - xi) : xi);

        valveRefXi_[i] = xi;
        valveRefXiLocal_[i] = xiLocal;
        valveRefCxi_[i] = Cxi;
        valveRefR0_[i] = r0;
    }

    mitralRefEnabled_ = true;
    mitralEnabled_ = false;   // avoid both models being active together

    InfoH << basic_Info
        << "Valve(slice-axis, reference-based) enabled for STL: " << stlPath_
        << " Ct=" << mCt_
        << " Cb=" << mCb_
        << " L="  << mLen_
        << " period=" << period_
        << " timeLaw=" << timeLaw_
        << " spaceLaw=" << spaceLaw_
        << " Fmax=" << Fmax_
        << endl;
}


void stlBased::applyValveSliceAxisRef(const scalar tNow)
{
    if (!mitralRefEnabled_) return;

    if (!useLegacyValveAxis_)
    {
        applyValveSliceAxis_saddle(tNow);
        return;
    }

    if (historyPoints_.empty()) return;
    if (historyPoints_.size() != valveRefXi_.size()) return;

    scalar phi = 0.0;
    if (period_ > VSMALL)
    {
        scalar tCycle = std::fmod(tNow, period_);
        if (tCycle < 0.0) tCycle += period_;
        phi = tCycle/period_;
    }

    const scalar Gt = valveTimeGain_(phi);

    const pointField& X0 = historyPoints_;
    pointField newPts(X0.size());

    label moved = 0;
    scalar minF = GREAT, maxF = -GREAT;

    forAll(X0, i)
    {
        const scalar xiLocal = valveRefXiLocal_[i];
        const scalar gS = valveSpaceGain_(xiLocal);
        const scalar F  = valveClamp01_(Fmax_ * gS * Gt);

        const vector Xnew = valveRefCxi_[i] + (1.0 - F)*valveRefR0_[i];
        newPts[i] = point(Xnew);

        if (mag(valveRefR0_[i]) > VSMALL && F > SMALL)
        {
            ++moved;
            minF = Foam::min(minF, F);
            maxF = Foam::max(maxF, F);
        }
    }

    if (moved == 0)
    {
        minF = 0.0;
        maxF = 0.0;
    }

    InfoH << iB_Info
        << "valve(slice-axis-ref): t=" << tNow
        << " phi=" << phi
        << " Gt=" << Gt
        << " movedPts=" << moved
        << " F[min,max]=[" << minF << "," << maxF << "]"
        << nl;

    bodySurfMesh_.movePoints(newPts);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}


vector stlBased::valveSliceAxisRefVelocityAtPoint
(
    const point& p,
    const scalar tNow
) const
{
    if (!mitralRefEnabled_ || period_ <= VSMALL)
    {
        return vector::zero;
    }

    const pointField curPts(bodySurfMesh_.points());

    if (!useLegacyValveAxis_)
    {
        if (curPts.empty() || historyPoints_.empty()
         || curPts.size() != historyPoints_.size())
        {
            return vector::zero;
        }

        label nearestI = 0;
        scalar minD2 = magSqr(vector(curPts[0] - p));
        for (label i = 1; i < curPts.size(); ++i)
        {
            const scalar d2 = magSqr(vector(curPts[i] - p));
            if (d2 < minD2)
            {
                minD2 = d2;
                nearestI = i;
            }
        }

        const label nPts = historyPoints_.size();
        const bool useAnnulusMask = (annulusMask_.size() == nPts);
        const scalar xiAttach = min(xiZeroTol_, 0.05*nearRamp_);

        if ((useAnnulusMask && annulusMask_[nearestI])
         || (!useAnnulusMask && valveXiLoaded_ && valveXi0_[nearestI] <= xiAttach))
        {
            return vector::zero;
        }

        vector axis = mCt_ - mCb_;
        const scalar Laxis = mag(axis);
        if (Laxis <= VSMALL)
        {
            return vector::zero;
        }
        const vector u = axis/Laxis;

        vector d0 = vector(historyPoints_[nearestI] - mCb_);
        scalar xiAxis = (d0 & u);
        xiAxis = max(scalar(0.0), min(Laxis, xiAxis));

        vector Cxi = vector(mCb_) + xiAxis*u;
        vector r0 = vector(historyPoints_[nearestI] - Cxi);
        r0 -= (r0 & u)*u;

        if (mag(r0) <= VSMALL)
        {
            return vector::zero;
        }

        const scalar xiLocal =
            (mFixedEnd_ == FixedEnd::Bottom ? xiAxis : (Laxis - xiAxis));
        const scalar gS = valveSpaceGain_(xiLocal);
        if (gS <= SMALL)
        {
            return vector::zero;
        }

        scalar phi = 0.0;
        scalar tCycle = std::fmod(tNow, period_);
        if (tCycle < 0.0) tCycle += period_;
        phi = tCycle/period_;

        const scalar F = valveClamp01_(Fmax_*valveTimeGain_(phi)*gS);
        // if (F >= 0.95)
        // {
        //     return vector::zero;
        // }

	 if (F > 1)
        {
            return vector::zero;
        }

        return -(Fmax_*gS*valveTimeGainDot_(tNow))*r0;
    }

    if
    (
        curPts.size() != valveRefR0_.size()
     || curPts.empty()
    )
    {
        WarningInFunction
            << "Reference-based mitral velocity data inconsistent for "
            << stlPath_ << ". Returning zero velocity."
            << endl;
        return vector::zero;
    }

    // Match current query point to current STL point index,
    // then use that SAME index to fetch frozen reference data.
    label nearestI = 0;
    scalar minD2 = magSqr(vector(curPts[0] - p));

    for (label i = 1; i < curPts.size(); ++i)
    {
        const scalar d2 = magSqr(vector(curPts[i] - p));
        if (d2 < minD2)
        {
            minD2 = d2;
            nearestI = i;
        }
    }

    const scalar xiLocal = valveRefXiLocal_[nearestI];
    const scalar gS = valveSpaceGain_(xiLocal);
    const scalar dGtDt = valveTimeGainDot_(tNow);

    return -(Fmax_ * gS * dGtDt) * valveRefR0_[nearestI];
}


//---------------------------------------------------------------------------//
// -------- Analytical helpers for deforming mitral slice-axis velocity ------

scalar stlBased::valveClamp01_(const scalar x) const
{
    return (x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x));
}

scalar stlBased::valveSmoothstep01_(const scalar u) const
{
    const scalar s = valveClamp01_(u);
    return 3.0*s*s - 2.0*s*s*s;
}

scalar stlBased::valveDSmoothstep01_(const scalar u) const
{
    if (u <= 0.0 || u >= 1.0) return 0.0;
    return 6.0*u - 6.0*u*u;
}

scalar stlBased::valveSmootherstep01_(const scalar u) const
{
    const scalar s = valveClamp01_(u);
    return s*s*s*(10.0 - 15.0*s + 6.0*s*s);
}

scalar stlBased::valveDSmootherstep01_(const scalar u) const
{
    if (u <= 0.0 || u >= 1.0) return 0.0;
    return 30.0*Foam::pow(u,4) - 60.0*Foam::pow(u,3) + 30.0*u*u;
}

// Rising edge from a -> b on periodic phase phi in [0,1)
scalar stlBased::valveEdgeUp_(const scalar phi, const scalar a, const scalar b) const
{
    const scalar aa = valveClamp01_(a);
    const scalar bb = valveClamp01_(b);

    scalar len = 0.0;
    scalar u   = -1.0;

    if (aa <= bb)
    {
        len = max(SMALL, bb - aa);
        if (phi <= aa) return 0.0;
        if (phi >= bb) return 1.0;
        u = (phi - aa)/len;
    }
    else
    {
        // wrapped interval [aa,1) U [0,bb]
        len = max(SMALL, (1.0 - aa) + bb);

        if (phi >= aa)
        {
            u = (phi - aa)/len;
        }
        else if (phi <= bb)
        {
            u = (phi + (1.0 - aa))/len;
        }
        else
        {
            return 0.0;
        }
    }

    return valveSmoothstep01_(u);
}

scalar stlBased::valveDEdgeUpDphi_(const scalar phi, const scalar a, const scalar b) const
{
    const scalar aa = valveClamp01_(a);
    const scalar bb = valveClamp01_(b);

    scalar len = 0.0;
    scalar u   = -1.0;

    if (aa <= bb)
    {
        len = max(SMALL, bb - aa);
        if (phi <= aa || phi >= bb) return 0.0;
        u = (phi - aa)/len;
    }
    else
    {
        len = max(SMALL, (1.0 - aa) + bb);

        if (phi >= aa)
        {
            u = (phi - aa)/len;
        }
        else if (phi <= bb)
        {
            u = (phi + (1.0 - aa))/len;
        }
        else
        {
            return 0.0;
        }
    }

    return valveDSmoothstep01_(u)/len;
}

// Exact match to edgeDown(...) used in applyValveSliceAxis(...)
scalar stlBased::valveEdgeDown_(const scalar phi, const scalar a, const scalar b) const
{
    const scalar aa = valveClamp01_(a);
    const scalar bb = valveClamp01_(b);

    if (aa <= bb)
    {
        const scalar len = max(SMALL, bb - aa);
        if (phi <= aa) return 0.0;
        if (phi >= bb) return 1.0;
        return valveSmoothstep01_((phi - aa)/len);
    }
    else
    {
        // wrapped version matches applyValveSliceAxis edgeDown exactly
        const scalar len = max(SMALL, 1.0 - aa);

        if (phi >= aa)
        {
            return valveSmoothstep01_((phi - aa)/len);
        }
        else if (phi <= bb)
        {
            return valveSmoothstep01_((phi + (1.0 - aa))/len);
        }
        else
        {
            return 1.0;
        }
    }
}

scalar stlBased::valveDEdgeDownDphi_(const scalar phi, const scalar a, const scalar b) const
{
    const scalar aa = valveClamp01_(a);
    const scalar bb = valveClamp01_(b);

    if (aa <= bb)
    {
        const scalar len = max(SMALL, bb - aa);
        if (phi <= aa || phi >= bb) return 0.0;
        const scalar u = (phi - aa)/len;
        return valveDSmoothstep01_(u)/len;
    }
    else
    {
        const scalar len = max(SMALL, 1.0 - aa);

        if (phi >= aa)
        {
            const scalar u = (phi - aa)/len;
            return valveDSmoothstep01_(u)/len;
        }
        else if (phi <= bb)
        {
            const scalar u = (phi + (1.0 - aa))/len;
            return valveDSmoothstep01_(u)/len;
        }
        else
        {
            return 0.0;
        }
    }
}

scalar stlBased::valveTimeGain_(const scalar phiIn) const
{
    scalar phi = phiIn;
    if (phi < 0.0) phi += 1.0;
    if (phi >= 1.0) phi = phi - std::floor(phi);

    scalar Gt = 0.0;

    if (timeLaw_ == "twoWindow")
    {
        const scalar rise = valveEdgeUp_(phi,   wClose1Beg_, wClose1End_);
        const scalar fall = valveEdgeDown_(phi, wOpenBeg_,   wOpenEnd_);
        Gt = valveClamp01_(rise * (1.0 - fall));
    }
    else if (timeLaw_ == "threeWindow")
    {
        const scalar r1 = valveEdgeUp_(phi,   wClose1Beg_, wClose1End_);
        const scalar fo = valveEdgeDown_(phi, wOpenBeg_,   wOpenEnd_);
        const scalar r2 = valveEdgeUp_(phi,   wClose2Beg_, wClose2End_);
        Gt = valveClamp01_(r1*(1.0 - fo) + r2);
    }
    else if (timeLaw_ == "cos2")
    {
        Gt = 0.5*(1.0 - Foam::cos(2.0*Foam::constant::mathematical::pi*phi));
    }
    else if (timeLaw_ == "smoothstep")
    {
        Gt = valveSmoothstep01_(phi);
    }
    else
    {
        const scalar s = Foam::sin(Foam::constant::mathematical::pi*phi);
        Gt = s*s;
    }

    return valveClamp01_(Gt);
}

scalar stlBased::valveTimeGainDot_(const scalar tNow) const
{
    if (period_ <= VSMALL) return 0.0;

    scalar tCycle = std::fmod(tNow, period_);
    if (tCycle < 0.0) tCycle += period_;

    const scalar phi    = tCycle/period_;
    const scalar dphidt = 1.0/period_;

    if (timeLaw_ == "twoWindow")
    {
        const scalar rise  = valveEdgeUp_(phi,   wClose1Beg_, wClose1End_);
        const scalar drise = valveDEdgeUpDphi_(phi,   wClose1Beg_, wClose1End_);

        const scalar fall  = valveEdgeDown_(phi, wOpenBeg_, wOpenEnd_);
        const scalar dfall = valveDEdgeDownDphi_(phi, wOpenBeg_, wOpenEnd_);

        const scalar dGtDphi = drise*(1.0 - fall) - rise*dfall;
        return dGtDphi*dphidt;
    }
    else if (timeLaw_ == "threeWindow")
    {
        const scalar r1   = valveEdgeUp_(phi,   wClose1Beg_, wClose1End_);
        const scalar dr1  = valveDEdgeUpDphi_(phi,   wClose1Beg_, wClose1End_);

        const scalar fo   = valveEdgeDown_(phi, wOpenBeg_, wOpenEnd_);
        const scalar dfo  = valveDEdgeDownDphi_(phi, wOpenBeg_, wOpenEnd_);

        const scalar dr2  = valveDEdgeUpDphi_(phi,   wClose2Beg_, wClose2End_);

        const scalar dGtDphi = dr1*(1.0 - fo) - r1*dfo + dr2;
        return dGtDphi*dphidt;
    }
    else if (timeLaw_ == "cos2")
    {
        return Foam::constant::mathematical::pi
             * Foam::sin(2.0*Foam::constant::mathematical::pi*phi)
             * dphidt;
    }
    else if (timeLaw_ == "smoothstep")
    {
        return valveDSmoothstep01_(phi)*dphidt;
    }
    else
    {
        // sin^2(pi*phi)
        return Foam::constant::mathematical::pi
             * Foam::sin(2.0*Foam::constant::mathematical::pi*phi)
             * dphidt;
    }
}

scalar stlBased::valveSpaceGain_(const scalar xiLocal) const
{
    if (xiLocal < 0.0) return 0.0;

    if (xiLocal <= nearRamp_)
    {
        const scalar s = valveClamp01_(xiLocal/max(VSMALL, nearRamp_));

        if (spaceLaw_ == "linear")     return s;
        if (spaceLaw_ == "smoothstep") return valveSmoothstep01_(s);

        return valveSmootherstep01_(s);
    }

    if (xiLocal <= activeLen_) return 1.0;

    return 0.0;
}

vector stlBased::valveSliceAxisVelocityAtPoint
(
    const point& p,
    const scalar tNow
) const
{
    if (!mitralEnabled_ || mLen_ <= VSMALL || period_ <= VSMALL)
    {
        return vector::zero;
    }

    if (!useLegacyValveAxis_)
    {
        const pointField curPts(bodySurfMesh_.points());

        if (curPts.empty() || historyPoints_.empty()
         || curPts.size() != historyPoints_.size())
        {
            return vector::zero;
        }

        label nearestI = 0;
        scalar minD2 = magSqr(vector(curPts[0] - p));
        for (label i = 1; i < curPts.size(); ++i)
        {
            const scalar d2 = magSqr(vector(curPts[i] - p));
            if (d2 < minD2)
            {
                minD2 = d2;
                nearestI = i;
            }
        }

        const label nPts = historyPoints_.size();
        const bool useAnnulusMask = (annulusMask_.size() == nPts);
        const scalar xiAttach = min(xiZeroTol_, 0.05*nearRamp_);

        if ((useAnnulusMask && annulusMask_[nearestI])
         || (!useAnnulusMask && valveXiLoaded_ && valveXi0_[nearestI] <= xiAttach))
        {
            return vector::zero;
        }

        vector axis = mCt_ - mCb_;
        const scalar Laxis = mag(axis);
        if (Laxis <= VSMALL)
        {
            return vector::zero;
        }
        const vector u = axis/Laxis;

        vector d0 = vector(historyPoints_[nearestI] - mCb_);
        scalar xiAxis = (d0 & u);
        xiAxis = max(scalar(0.0), min(Laxis, xiAxis));

        vector Cxi = vector(mCb_) + xiAxis*u;
        vector r0 = vector(historyPoints_[nearestI] - Cxi);
        r0 -= (r0 & u)*u;

        if (mag(r0) <= VSMALL)
        {
            return vector::zero;
        }

        const scalar xiLocal =
            (mFixedEnd_ == FixedEnd::Bottom ? xiAxis : (Laxis - xiAxis));
        const scalar gS = valveSpaceGain_(xiLocal);
        if (gS <= SMALL)
        {
            return vector::zero;
        }

        scalar tCycle = std::fmod(tNow, period_);
        if (tCycle < 0.0) tCycle += period_;
        const scalar phi = tCycle/period_;

        const scalar F = valveClamp01_(Fmax_*valveTimeGain_(phi)*gS);
        // if (F >= 0.95)
        // {
        //     return vector::zero;
        // }

	 if (F > 1)
        {
            return vector::zero;
        }

        return -(Fmax_*gS*valveTimeGainDot_(tNow))*r0;
    }

    // Current axial coordinate along Ct -> Cb axis
    vector d = vector(p - mCt_);
    scalar xi = (d & mAxisUnit_);
    xi = max(scalar(0.0), min(mLen_, xi));

    const vector Cxi = vector(mCt_) + xi*mAxisUnit_;

    // Current in-plane radial vector
    vector r = vector(p - Cxi);
    r -= (r & mAxisUnit_) * mAxisUnit_;

    const scalar xiLocal =
        (mFixedEnd_ == FixedEnd::Bottom ? (mLen_ - xi) : xi);

    const scalar gS = valveSpaceGain_(xiLocal);

    scalar phi = 0.0;
    if (period_ > VSMALL)
    {
        scalar tCycle = std::fmod(tNow, period_);
        if (tCycle < 0.0) tCycle += period_;
        phi = tCycle / period_;
    }

    const scalar Gt    = valveTimeGain_(phi);
    const scalar dGtDt = valveTimeGainDot_(tNow);

    const scalar F    = valveClamp01_(Fmax_ * gS * Gt);
    const scalar dFdt = Fmax_ * gS * dGtDt;

    // Exact Eulerian volumetric velocity:
    // x = C(xi) + (1 - F) r0  =>  r = (1 - F) r0
    // u = -dFdt * r0 = -(dFdt/(1-F)) * r
    const scalar oneMinusF = max(scalar(1e-6), scalar(1.0 - F));

    return -(dFdt / oneMinusF) * r;
}
//--------------------------------------------------------------------------//
void stlBased::readBreathingSphere(const dictionary& d)
{
    // Required entries:
    //  center      (x y z)
    //  baseRadius  r
    //  period      T

    breathingCenter_ = vector(d.lookup("center"));
    breathingBaseR_  = readScalar(d.lookup("baseRadius"));
    breathingPeriod_ = readScalar(d.lookup("period"));

    // Store the current STL points as reference (rest) shape
    breathingRefPts_ = bodySurfMesh_.points();

    breathingEnabled_ = true;

    Info<< "Breathing sphere motion enabled for " << stlPath_
        << " center=" << breathingCenter_
        << " baseRadius=" << breathingBaseR_
        << " period=" << breathingPeriod_
        << endl;
}
void stlBased::applyBreathingSphere(const Time& runTime)
{
    if (!breathingEnabled_) return;
    if (breathingPeriod_ <= SMALL) return;
    if (breathingRefPts_.empty()) return;

    const scalar T   = breathingPeriod_;
    const scalar t   = runTime.value();
    const scalar tau = std::fmod(t, T);  // [0, T)
    const scalar s   =
        1.0 + Foam::sin(constant::mathematical::pi * tau / T);
    // s: 1 -> 2 -> 1 over one period

    pointField newPts(breathingRefPts_.size());

    forAll(breathingRefPts_, i)
    {
        const vector d = breathingRefPts_[i] - breathingCenter_;
        newPts[i] = breathingCenter_ + s * d;
    }

    bodySurfMesh_.movePoints(newPts);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));

    // Optional debug:
    scalar maxDisp = 0.0;
    forAll(newPts, i)
        maxDisp = max(maxDisp, mag(newPts[i] - breathingRefPts_[i]));
    Info<< "Breathing@" << stlPath_
        << " t=" << t
        << " s=" << s
        << " max|dP|=" << maxDisp
        << endl;
}
//---------------------------------------------------------------------------//
//--------------------------------------------------------------------------//
void stlBased::readSinusoidalTranslation(const dictionary& d)
{
    // Dictionary entries:
    // amplitude  A;          [m]
    // period     T;          [s]  (default 1.0)
    // direction (dx dy dz);  [-]  (default (1 0 0))
    // phase      phi;        [rad] (default 0.0)

    sinusAmp_    = d.lookupOrDefault<scalar>("amplitude", 0.0);
    sinusPeriod_ = d.lookupOrDefault<scalar>("period", 1.0);
    sinusPhase_  = d.lookupOrDefault<scalar>("phase", 0.0);
    sinusDir_    = d.lookupOrDefault<vector>("direction", vector(1,0,0));
    sinusVerticalShift_ = d.lookupOrDefault<vector>("verticalShift", vector::zero);

    const scalar dirMag = mag(sinusDir_);
    if (dirMag > VSMALL)
    {
        sinusDir_ /= dirMag; // normalize
    }
    else
    {
        sinusDir_ = vector(1,0,0);
    }

    // Store current STL points as the reference (rest) shape
    sinusRefPts_ = bodySurfMesh_.points();

    // Enable only if amplitude is non-zero (but you can keep enabled even if 0)
    sinusoidalEnabled_ = true;

    Info<< "Sinusoidal translation enabled for " << stlPath_
        << " A=" << sinusAmp_
        << " T=" << sinusPeriod_
        << " phase=" << sinusPhase_
        << " dir=" << sinusDir_
        << " verticalShift=" << sinusVerticalShift_
        << endl;
}

vector stlBased::sinusoidalVelocity(const scalar tNow) const
{
    if (!sinusoidalEnabled_ || sinusPeriod_ <= SMALL || mag(sinusAmp_) <= VSMALL)
    {
        return vector::zero;
    }

    const scalar omega = 2.0*Foam::constant::mathematical::pi / sinusPeriod_;
    const scalar vel   = sinusAmp_ * omega * Foam::cos(omega*tNow + sinusPhase_);

    return vel * sinusDir_;
}

void stlBased::applySinusoidalTranslation(const scalar tNow)
{
    if (!sinusoidalEnabled_) return;
    if (sinusRefPts_.empty()) return;
    if (sinusPeriod_ <= SMALL) return;

    if (mag(sinusAmp_) <= VSMALL)
    {
        return;
    }

    const scalar omega = 2.0*Foam::constant::mathematical::pi / sinusPeriod_;
    const scalar disp  = sinusAmp_ * Foam::sin(omega*tNow + sinusPhase_);
    const vector dX    = disp * sinusDir_ + sinusVerticalShift_;

    pointField newPts(sinusRefPts_.size());
    forAll(sinusRefPts_, i)
    {
        newPts[i] = sinusRefPts_[i] + dX;
    }

    bodySurfMesh_.movePoints(newPts);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}

void stlBased::applySinusoidalTranslation(const Time& runTime)
{
    applySinusoidalTranslation(runTime.value());
}

void stlBased::readQuadraticBend(const dictionary& d)
{
    quadraticBendEnabled_ = true;

    quadraticBendAmp_    = readScalar(d.lookup("amplitude"));
    quadraticBendPeriod_ = readScalar(d.lookup("period"));
    quadraticBendPhase_  = d.found("phase") ? readScalar(d.lookup("phase")) : 0.0;

    quadraticBendYMin_   = readScalar(d.lookup("yMin"));
    quadraticBendYMax_   = readScalar(d.lookup("yMax"));
    quadraticBendHeight_ = quadraticBendYMax_ - quadraticBendYMin_;

    quadraticBendDir_ =
        d.found("direction")
      ? vector(d.lookup("direction"))
      : vector(1,0,0);

    const scalar dirMag = mag(quadraticBendDir_);
    if (dirMag > VSMALL)
    {
        quadraticBendDir_ /= dirMag;
    }
    else
    {
        quadraticBendDir_ = vector(1,0,0);
    }

    quadraticBendRefPts_ = bodySurfMesh_.points();
    quadraticBendXi_.setSize(quadraticBendRefPts_.size());

    forAll(quadraticBendRefPts_, i)
    {
        scalar xi =
            (quadraticBendRefPts_[i].y() - quadraticBendYMin_)
           /(quadraticBendHeight_ + SMALL);

        xi = max(scalar(0), min(scalar(1), xi));
        quadraticBendXi_[i] = xi;
    }

    Info << "readQuadraticBend: amplitude = " << quadraticBendAmp_
         << " period = " << quadraticBendPeriod_
         << " phase = " << quadraticBendPhase_
         << " yMin = " << quadraticBendYMin_
         << " yMax = " << quadraticBendYMax_
         << " direction = " << quadraticBendDir_
         << endl;
}

void stlBased::applyQuadraticBend(const scalar tNow)
{
    if (!quadraticBendEnabled_) return;
    if (quadraticBendRefPts_.empty()) return;
    if (quadraticBendPeriod_ <= SMALL) return;

    const scalar omega =
        2.0*Foam::constant::mathematical::pi/quadraticBendPeriod_;

    const scalar s = Foam::sin(omega*tNow + quadraticBendPhase_);

    pointField newPts(quadraticBendRefPts_.size());

    forAll(quadraticBendRefPts_, i)
    {
        const scalar xi = quadraticBendXi_[i];
        const scalar disp = quadraticBendAmp_ * sqr(xi) * s;

        newPts[i] = quadraticBendRefPts_[i] + disp*quadraticBendDir_;
    }

    bodySurfMesh_.movePoints(newPts);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}

void stlBased::applyQuadraticBend(const Time& t)
{
    applyQuadraticBend(t.value());
}

vector stlBased::quadraticBendVelocityAtPoint
(
    const point& p,
    const scalar tNow
) const
{
    if (!quadraticBendEnabled_ || quadraticBendPeriod_ <= SMALL)
    {
        return vector::zero;
    }

    scalar xi = (p.y() - quadraticBendYMin_)/(quadraticBendHeight_ + SMALL);
    xi = max(scalar(0), min(scalar(1), xi));

    const scalar omega =
        2.0*Foam::constant::mathematical::pi/quadraticBendPeriod_;

    const scalar velMag =
        quadraticBendAmp_ * sqr(xi) * omega
       *Foam::cos(omega*tNow + quadraticBendPhase_);

    return velMag * quadraticBendDir_;
}

//================== custom profile bend, 11 march 26 ======================

void stlBased::readCustomProfileBend(const dictionary& d)
{
    customProfileBendEnabled_ = true;

    customProfileAmpX_   = readScalar(d.lookup("amplitudeX"));
    customProfileAmpY_   = readScalar(d.lookup("amplitudeY"));
    customProfilePeriod_ = readScalar(d.lookup("period"));
    customProfilePhase_  = d.found("phase") ? readScalar(d.lookup("phase")) : 0.0;

    customProfileYMin_   = readScalar(d.lookup("yMin"));
    customProfileYMax_   = readScalar(d.lookup("yMax"));
    customProfileHeight_ = customProfileYMax_ - customProfileYMin_;

    customProfileXCenter_ = readScalar(d.lookup("xCenter"));

    customProfileRefPts_ = bodySurfMesh_.points();
    customProfileXi_.setSize(customProfileRefPts_.size());

    forAll(customProfileRefPts_, i)
    {
        scalar xi =
            (customProfileRefPts_[i].y() - customProfileYMin_)
           /(customProfileHeight_ + SMALL);

        xi = max(scalar(0), min(scalar(1), xi));
        customProfileXi_[i] = xi;
    }

    Info << "readCustomProfileBend:"
         << " amplitudeX = " << customProfileAmpX_
         << " amplitudeY = " << customProfileAmpY_
         << " period = " << customProfilePeriod_
         << " phase = " << customProfilePhase_
         << " yMin = " << customProfileYMin_
         << " yMax = " << customProfileYMax_
         << " xCenter = " << customProfileXCenter_
         << endl;
}


void stlBased::applyCustomProfileBend(const scalar tNow)
{
    if (!customProfileBendEnabled_) return;
    if (customProfileRefPts_.empty()) return;
    if (customProfilePeriod_ <= SMALL) return;

    const scalar omega =
        2.0*Foam::constant::mathematical::pi/customProfilePeriod_;

    const scalar sTime = Foam::sin(omega*tNow + customProfilePhase_);

    pointField newPts(customProfileRefPts_.size());

    forAll(customProfileRefPts_, i)
    {
        const point& X0 = customProfileRefPts_[i];
        const scalar xi = customProfileXi_[i];

	const scalar g =
	  6.0*Foam::pow(xi, 5)
	  - 15.0*Foam::pow(xi, 4)
	  + 10.0*Foam::pow(xi, 3);

	const scalar gp =
	  30.0*Foam::pow(xi, 4)
	  - 60.0*Foam::pow(xi, 3)
	  + 30.0*Foam::sqr(xi);

        const scalar xShift =
            customProfileAmpX_ * g * sTime;

        const scalar yShift =
           -(X0.x() - customProfileXCenter_)
            * (customProfileAmpY_/(customProfileHeight_ + SMALL))
            * gp * sTime;

        newPts[i] = point
        (
            X0.x() + xShift,
            X0.y() + yShift,
            X0.z()
        );
    }

    bodySurfMesh_.movePoints(newPts);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}


void stlBased::applyCustomProfileBend(const Time& t)
{
    applyCustomProfileBend(t.value());
}


vector stlBased::customProfileBendVelocityAtPoint
(
    const point& p,
    const scalar tNow
) const
{
    if (!customProfileBendEnabled_ || customProfilePeriod_ <= SMALL)
    {
        return vector::zero;
    }

    scalar xi =
        (p.y() - customProfileYMin_) / (customProfileHeight_ + SMALL);

    xi = max(scalar(0), min(scalar(1), xi));

    // SAME profile as applyCustomProfileBend(...)
    const scalar g =
      6.0*Foam::pow(xi, 5)
      - 15.0*Foam::pow(xi, 4)
      + 10.0*Foam::pow(xi, 3);

    const scalar gp =
      30.0*Foam::pow(xi, 4)
      - 60.0*Foam::pow(xi, 3)
      + 30.0*Foam::sqr(xi);

    const scalar omega =
        2.0*Foam::constant::mathematical::pi/customProfilePeriod_;

    const scalar cTime = Foam::cos(omega*tNow + customProfilePhase_);

    const scalar ux =
        customProfileAmpX_ * g * omega * cTime;

    const scalar uy =
       -(p.x() - customProfileXCenter_)
        * (customProfileAmpY_/(customProfileHeight_ + SMALL))
        * gp * omega * cTime;

    return vector(ux, uy, 0.0);
}


void stlBased::customProfileBendMaxVelocity
(
    const scalar tNow,
    scalar& maxUx,
    scalar& maxUy,
    scalar& maxUmag
) const
{
    maxUx = 0.0;
    maxUy = 0.0;
    maxUmag = 0.0;

    if (!customProfileBendEnabled_ || customProfileRefPts_.empty())
    {
        return;
    }

    forAll(customProfileRefPts_, i)
    {
        const vector v = customProfileBendVelocityAtPoint(customProfileRefPts_[i], tNow);

        maxUx = max(maxUx, mag(v.x()));
        maxUy = max(maxUy, mag(v.y()));
        maxUmag = max(maxUmag, mag(v));
    }
}

//---------------------------------------------------------------------------//
void stlBased::getIntersectionPoints
(
    const label index,
    const treeBoundBox& cubeBb,
    DynamicPointList& intersectionPoints
)
{
    const pointField& points = triSurf_->points();
    const typename triSurface::FaceType& f = (*triSurf_)[index];

    for (auto ind : f)
    {
        if (cubeBb.contains(points[ind]))
        {
            intersectionPoints.append(points[ind]);
        }
    }

    const point fc = f.centre(points);

    if (f.size() == 3)
    {
        return intersectBb
        (
            points[f[0]],
            points[f[1]],
            points[f[2]],
            cubeBb,
            intersectionPoints
        );
    }
    else
    {
        forAll(f, fp)
        {
            intersectBb
            (
                points[f[fp]],
                points[f[f.fcIndex(fp)]],
                fc,
                cubeBb,
                intersectionPoints
            );
        }
    }

    return;
}
//---------------------------------------------------------------------------//
void stlBased::intersectBb
(
    const point& p0,
    const point& p1,
    const point& p2,
    const treeBoundBox& cubeBb,
    DynamicPointList& intersectionPoints
)
{
    const vector p10 = p1 - p0;
    const vector p20 = p2 - p0;

    // cubeBb points; counted as if cell with vertex0 at cubeBb.min().
    const point& min = cubeBb.min();
    const point& max = cubeBb.max();

    const point& cube0 = min;
    const point  cube1(min.x(), min.y(), max.z());
    const point  cube2(max.x(), min.y(), max.z());
    const point  cube3(max.x(), min.y(), min.z());

    const point  cube4(min.x(), max.y(), min.z());
    const point  cube5(min.x(), max.y(), max.z());
    const point  cube7(max.x(), max.y(), min.z());

    //
    // Intersect all 12 edges of cube with triangle
    //

    point pInter;
    pointField origin(4);
    // edges in x direction
    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube5;
    origin[3] = cube4;

    scalar maxSx = max.x() - min.x();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 0, origin, maxSx, pInter))
    {
        intersectionPoints.append(pInter);
    }

    // edges in y direction
    origin[0] = cube0;
    origin[1] = cube1;
    origin[2] = cube2;
    origin[3] = cube3;

    scalar maxSy = max.y() - min.y();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 1, origin, maxSy, pInter))
    {
        intersectionPoints.append(pInter);
    }

    // edges in z direction
    origin[0] = cube0;
    origin[1] = cube3;
    origin[2] = cube7;
    origin[3] = cube4;

    scalar maxSz = max.z() - min.z();

    if (triangleFuncs::intersectAxesBundle(p0, p10, p20, 2, origin, maxSz, pInter))
    {
        intersectionPoints.append(pInter);
    }


    // Intersect triangle edges with bounding box
    if (cubeBb.intersects(p0, p1, pInter))
    {
        intersectionPoints.append(pInter);
    }
    if (cubeBb.intersects(p1, p2, pInter))
    {
        intersectionPoints.append(pInter);
    }
    if (cubeBb.intersects(p2, p0, pInter))
    {
        intersectionPoints.append(pInter);
    }
}
//---------------------------------------------------------------------------//
void stlBased::setBodyPosition(pointField pos)
{
    bodySurfMesh_.movePoints(pos);
    triSurf_.reset(new triSurface(bodySurfMesh_));
    triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
}
//---------------------------------------------------------------------------//
