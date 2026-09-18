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
    class for addition and removal of immersedBodies

SourceFiles
    stlBased.C

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef stlBased_H
#define stlBased_H

#include "geomModel.H"

namespace Foam
{

// Forward declaration
class Time;

class stlBased
:
    public geomModel
{
    private:

        typedef DynamicList<label>  DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point>  DynamicPointList;

        // --- Valve motion state --------------------------------------
        bool   mitralEnabled_{false};

        enum class FixedEnd { Top, Bottom };
        FixedEnd mFixedEnd_{ FixedEnd::Top };

        vector mCt_, mCb_;      // top & bottom centres
        vector mAxisUnit_;      // unit axis (Ct -> Cb or reverse)
        scalar mLen_{0.0};      // axis length

        // Timing
        scalar period_{1.0};
        word   timeLaw_{"twoWindow"};   // "twoWindow" | "threeWindow" | others
        scalar wClose1Beg_{0.0}, wClose1End_{0.0};
        scalar wOpenBeg_{0.0},   wOpenEnd_{0.0};
        scalar wClose2Beg_{0.0}, wClose2End_{0.0};   // for "threeWindow"

        // Spatial/gains
        scalar Fmax_{0.7};
        scalar nearRamp_{0.004};        // [m]
        scalar activeLen_{0.016};       // [m]
        word   spaceLaw_{"smootherstep"}; // "linear"|"smoothstep"|"smootherstep"

        // Misc
        scalar   leafThickness_{0.002};
        bool     flipNormalTowardValve_{true};
        fileName valveRefFile_;

        // Select the flat slice-axis or saddle-aware valve model.
        bool     useLegacyValveAxis_{false};

        // Per-vertex distance from the fixed annulus.
        scalarField valveXi0_;
        bool        valveXiLoaded_{false};

        fileName xiFieldFile_;

        // Saddle-mode xi/annulus attachment state
        scalar   xiZeroTol_{0.0};
        scalar   xiMaxUsed_{0.0};
        boolList annulusMask_;
        scalar   annulusTol_{0.0};

        // Annulus geometry for the saddle-aware model.
        autoPtr<triSurfaceMesh>    annulusSurfMesh_;
        autoPtr<triSurface>        annulusTriSurf_;
        autoPtr<triSurfaceSearch>  annulusSearch_;

        // Nearest annulus point for each leaflet vertex.
        pointField annulusRef_;

       // Reference-based mitral deformation/velocity mode.
       bool mitralRefEnabled_ = false;

       // Material/reference data for each STL point.
      scalarField valveRefXi_;
      scalarField valveRefXiLocal_;
      vectorField valveRefCxi_;
       vectorField valveRefR0_;
  
        // --- Breathing sphere state ----------------------------------
        bool   breathingEnabled_{false};
        vector breathingCenter_{vector::zero};
        scalar breathingBaseR_{0.0};
        scalar breathingPeriod_{1.0};
        pointField breathingRefPts_;   // reference (rest) STL points


        // --- Sinusoidal translation motion state ------------------------
        bool      sinusoidalEnabled_{false};
        scalar    sinusAmp_{0.0};          // [m]
        scalar    sinusPeriod_{1.0};       // [s]
        scalar    sinusPhase_{0.0};        // [rad]
        vector    sinusDir_{vector(1,0,0)}; // direction (will be normalized)
        vector    sinusVerticalShift_{vector::zero}; // optional constant offset [m]
        pointField sinusRefPts_;           // reference (rest) STL points

        //------------- leaflet bend -----------------------

        // Quadratic bending motion (reference-coordinate based)
        bool   quadraticBendEnabled_{false};

        scalar quadraticBendAmp_{0.0};
        scalar quadraticBendPeriod_{1.0};
        scalar quadraticBendPhase_{0.0};

        scalar quadraticBendYMin_{0.0};
        scalar quadraticBendYMax_{0.0};
        scalar quadraticBendHeight_{0.0};

        vector quadraticBendDir_{vector(1,0,0)};

        pointField quadraticBendRefPts_;
        scalarField quadraticBendXi_;

  //---------------------- custom bend, 11 march 26---------------

          // Custom profile bend (2D) with both Ux and Uy
        bool   customProfileBendEnabled_{false};

        scalar customProfileAmpX_{0.0};
        scalar customProfileAmpY_{0.0};
        scalar customProfilePeriod_{1.0};
        scalar customProfilePhase_{0.0};

        scalar customProfileYMin_{0.0};
        scalar customProfileYMax_{0.0};
        scalar customProfileHeight_{0.0};

        scalar customProfileXCenter_{0.0};

        pointField customProfileRefPts_;
        scalarField customProfileXi_;


  
    protected:

        triSurfaceMesh      bodySurfMesh_;
        word                stlPath_;
        autoPtr<triSurface> triSurf_;
        autoPtr<triSurfaceSearch> triSurfSearch_;

        pointField historyPoints_;

        void getIntersectionPoints
        (
            const label index,
            const treeBoundBox& cubeBb,
            DynamicPointList& intersectionPoints
        );

        void intersectBb
        (
            const point& p0,
            const point& p1,
            const point& p2,
            const treeBoundBox& cubeBb,
            DynamicPointList& intersectionPoints
        );

        // Flat slice-axis valve model.
        void readValveSliceAxis_legacy(const dictionary& d);
        void applyValveSliceAxis_legacy(const Time& t);

        // Saddle-aware valve model using xiFieldFile and annulus STL.
        void readValveSliceAxis_saddle(const dictionary& d);
        void applyValveSliceAxis_saddle(const Time& t);
        void applyValveSliceAxis_saddle(const scalar tNow);

        // Load the per-vertex annulus-distance field.
        void readValveXiField(const fileName& fName);

        // Build annulusMask_ from an annulus-only STL.
        void buildAnnulusMask(const fileName& annulusFile);

        // Pointwise analytical velocity helpers for deforming valves.
        scalar valveClamp01_(const scalar x) const;
        scalar valveSmoothstep01_(const scalar u) const;
        scalar valveDSmoothstep01_(const scalar u) const;
        scalar valveSmootherstep01_(const scalar u) const;
        scalar valveDSmootherstep01_(const scalar u) const;

        scalar valveEdgeUp_(const scalar phi, const scalar a, const scalar b) const;
        scalar valveDEdgeUpDphi_(const scalar phi, const scalar a, const scalar b) const;
        scalar valveEdgeDown_(const scalar phi, const scalar a, const scalar b) const;
        scalar valveDEdgeDownDphi_(const scalar phi, const scalar a, const scalar b) const;
        scalar valveTimeGain_(const scalar phi) const;

        scalar valveTimeGainDot_(const scalar tNow) const;
        scalar valveSpaceGain_(const scalar xiLocal) const;

    public:

        // Constructors ------------------------------------------------
        stlBased
        (
            const fvMesh&   mesh,
            // const contactType cType,
            word            stlPath,
            scalar          thrSurf
        );

        // Destructors -------------------------------------------------
        virtual ~stlBased() {}

        virtual void createImmersedBody
        (
            volScalarField& body,
            Field<label>&   octreeField,
            List<labelList>& cellPoints
        ) override
        {
            InfoH << iB_Info << "-- geomModelMessage-- "
                  << "notImplemented" << endl;
        }

        virtual std::shared_ptr<geomModel> getCopy() override = 0;

        vector addModelReturnRandomPosition
        (
            const bool      allActiveCellsInMesh,
            const boundBox  cellZoneBounds,
            Random&         randGen
        ) final override;

        void bodyMovePoints(vector translVec) final override;
        void bodyScalePoints(scalar scaleFac) final override;

        void bodyRotatePoints
        (
            scalar rotAngle,
            vector axisOfRot
        ) final override;

        // Valve (slice-axis) motion API -------------------------------
        bool mitralEnabled() const { return mitralEnabled_; }
        void readValveSliceAxis(const dictionary& d);   // parse sub-dict once
        void applyValveSliceAxis(const Time& t);        // deform STL for current time

  //--------------- new, 19 march 26

         //------------ valve definitions -----------------------
      // New reference-based mitral mode
  void readValveSliceAxisRef(const dictionary& d);
  void applyValveSliceAxisRef(const scalar tNow);

  bool mitralRefEnabled() const
  {
    return mitralRefEnabled_;
  }

  bool hasAnalyticalMitralRefVelocity() const
  {
    return mitralRefEnabled_;
  }

  vector valveSliceAxisRefVelocityAtPoint
  (
   const point& p,
   const scalar tNow
   ) const;

        // Breathing sphere API ----------------------------------------
        bool breathingEnabled() const { return breathingEnabled_; }
        void readBreathingSphere(const dictionary& d);
        void applyBreathingSphere(const Time& t);


       // Sinusoidal translation API -------------------------------------
       bool sinusoidalEnabled() const { return sinusoidalEnabled_; }
       void readSinusoidalTranslation(const dictionary& d);
       void applySinusoidalTranslation(const Time& t);

  // Evaluate sinusoidal translation at an explicit time.
  void applySinusoidalTranslation(const scalar tNow);

  // Analytical rigid translation velocity for the sinusoidal case.
  vector sinusoidalVelocity(const scalar tNow) const;

  bool hasAnalyticalSinusoidalVelocity() const { return sinusoidalEnabled_; }

       // Quadratic bending API ----------------------------------------
       bool quadraticBendEnabled() const { return quadraticBendEnabled_; }

       void readQuadraticBend(const dictionary& d);
       void applyQuadraticBend(const Time& t);
       void applyQuadraticBend(const scalar tNow);

       // Pointwise analytical velocity for deforming body
       vector quadraticBendVelocityAtPoint(const point& p, const scalar tNow) const;

       bool hasAnalyticalQuadraticBendVelocity() const
       {
           return quadraticBendEnabled_;
       }

        // Pointwise analytical velocity for deforming mitral slice-axis model.
        vector valveSliceAxisVelocityAtPoint(const point& p, const scalar tNow) const;

        bool hasAnalyticalMitralVelocity() const
        {
            return mitralEnabled_;
        }
       // Custom profile bend API ----------------------------------------
       bool customProfileBendEnabled() const { return customProfileBendEnabled_; }

       void readCustomProfileBend(const dictionary& d);
       void applyCustomProfileBend(const Time& t);
       void applyCustomProfileBend(const scalar tNow);

       vector customProfileBendVelocityAtPoint(const point& p, const scalar tNow) const;

       bool hasAnalyticalCustomProfileBendVelocity() const
       {
           return customProfileBendEnabled_;
       }

       void customProfileBendMaxVelocity
       (
           const scalar tNow,
           scalar& maxUx,
           scalar& maxUy,
           scalar& maxUmag
       ) const;
  
  //----------------------------------------------------------------
        // ---------------------- geomModel interface ------------------
        void resetBody() final override
        {
            triSurfaceMesh tmpSurf
            (
                IOobject
                (
                    stlPath_,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            bodySurfMesh_.movePoints(tmpSurf.points());
        }

        vector getCoM() final override
        {
            const pointField bSMeshPts(bodySurfMesh_.points());
            vector CoM(vector::zero);
            forAll(bSMeshPts, point)
            {
                CoM += bSMeshPts[point];
            }
            CoM /= bSMeshPts.size();
            return CoM;
        }

        boundBox getBounds() final override
        {
            boundBox bounds(bodySurfMesh_.points(), false);
            return bounds;
        }

        void synchronPos(label owner = -1) final override;

        void writeBodyInfo(label bodyId, IOdictionary& dict, word path) final override
        {
            triSurface triToRet(bodySurfMesh_);
            triToRet.write(path + "/stlFiles/" + name(bodyId) + ".stl", ".stl");
        }

        boolList pointInside(pointField pointI) final override
        {
            boolList returnList = triSurfSearch_().calcInside(pointI);
            return returnList;
        }

        bool pointInside(point pointI) final override
        {
            pointField points(1, pointI);
            boolList returnList = triSurfSearch_().calcInside(points);
            return returnList[0];
        }

        // volumeType getVolumeType(subVolume& sv, bool cIb) final override;

        void getClosestPointAndNormal
        (
            const point&  startPoint,
            const vector& span,
            point&        closestPoint,
            vector&       normal
        ) final override
        {
            pointIndexHit ibPointIndexHit = triSurfSearch_().nearest(startPoint, span);

            if (!ibPointIndexHit.hit())
            {
                // Moving/deforming STLs can leave the initial local span too
                // small for a few interface cells. Retry locally before using
                // the conservative CoM fallback.
                const vector retrySpan(4.0*span);
                ibPointIndexHit =
                    triSurfSearch_().nearest(startPoint, retrySpan);
            }

            if (ibPointIndexHit.hit())
            {
                List<pointIndexHit> ibPointIndexHitList(1, ibPointIndexHit);
                vectorField normalVectorField;

                const triSurfaceMesh& ibTempMesh(bodySurfMesh_);
                ibTempMesh.getNormal(ibPointIndexHitList, normalVectorField);

                normal       = normalVectorField[0];
                closestPoint = ibPointIndexHit.hitPoint();
            }
            else
            {
                InfoH << basic_Info << "Missing the closest point!" << endl;
                normal       = startPoint - getCoM();
                closestPoint = getCoM();
            }
        }

        // bool limitFinalSubVolume
        // (
        //     const subVolume& sv,
        //     bool cIb,
        //     boundBox& limBBox
        // ) final override;

        void initVars(pointField surfPoints, scalar thrSurf, scalar intSpan, bool sdBasedLambda)
        {
            bodySurfMesh_.movePoints(surfPoints);
            triSurf_.reset(new triSurface(bodySurfMesh_));
            triSurfSearch_.reset(new triSurfaceSearch(triSurf_()));
            historyPoints_ = surfPoints;
            thrSurf_       = thrSurf;
            intSpan_       = intSpan;
            sdBasedLambda_ = sdBasedLambda;
        }

        void setCoM() final override
        {
            const pointField bSMeshPts(bodySurfMesh_.points());
            CoM_ = vector::zero;
            forAll(bSMeshPts, point)
            {
                CoM_ += bSMeshPts[point];
            }
            CoM_ /= bSMeshPts.size();
        }

        pointField getBodyPoints() final override // Philip: this function should be changed to return a tmp!
        {
            return bodySurfMesh_.points()();
        }

        void setBodyPosition(pointField pos) final override;
};

} // End namespace Foam

#endif
