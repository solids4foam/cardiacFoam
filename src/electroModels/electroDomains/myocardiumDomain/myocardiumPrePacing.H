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
    Foam::myocardiumPrePacing

Description
    Single-cell pre-pacing of the myocardium's main tissues, driven by
    constant/prePacingProperties. A main tissue is the 'tissue' entry, or
    each ionicHeterogeneity region with its own constants and overrides.
    Blend-transition and gradientAxes-scaled cells keep their own constants
    and start from their dominant region's state.

    Registered on the myocardium mesh as "myocardiumPrePacing" so coupled
    models (active tension) can pace their own single cells against the
    same main tissues.

SourceFiles
    myocardiumPrePacing.C
\*---------------------------------------------------------------------------*/

#ifndef myocardiumPrePacing_H
#define myocardiumPrePacing_H

#include "regIOobject.H"
#include "fvMesh.H"
#include "volFields.H"
#include "ionicModel.H"
#include "prePacingIO.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class myocardiumPrePacing Declaration
\*---------------------------------------------------------------------------*/

class myocardiumPrePacing
:
    public regIOobject
{
    // Private Data

        const fvMesh& mesh_;

        //- Ionic model dictionary of the myocardium
        const dictionary electroProperties_;

        //- Case time step, the single-cell default [s]
        const scalar dt_;

        //- Main tissue names
        wordList regionNames_;

        //- Heterogeneity field value of each main tissue's single cell
        scalarList regionFieldValue_;

        //- Main tissue of every myocardium cell
        labelList cellRegion_;

        //- Cells in namedRegions blend transitions
        label nBlendedCells_;

        //- ionicHeterogeneity with transitions disabled, for single cells
        dictionary pureHeterogeneityDict_;
        bool hasHeterogeneity_;

        //- prePacingProperties settings of every main tissue, read once
        //  (a collective read in parallel)
        List<prePacingIO::PrePacingConfig> configs_;

        //- Main tissues paced so far and their converged ionic states
        boolList regionPaced_;
        List<scalarField> ionicStates_;


public:

    TypeName("myocardiumPrePacing");


    // Constructors

        //- Construct from the myocardium mesh and ionic dictionary.
        //  heterogeneityFieldValues holds, per myocardium cell, the
        //  ionicHeterogeneity field (region index for cellZoneRegions, t for
        //  namedRegions); nullptr without ionicHeterogeneity.
        myocardiumPrePacing
        (
            const fvMesh& mesh,
            const dictionary& electroProperties,
            const scalar dt,
            const scalarField* heterogeneityFieldValues,
            const label nCells
        );


    // Member Functions

        const wordList& regionNames() const
        {
            return regionNames_;
        }

        const labelList& cellRegion() const
        {
            return cellRegion_;
        }

        bool paced(const label regionI) const
        {
            return regionPaced_[regionI];
        }

        //- Single-cell time step of a main tissue [s]
        scalar singleCellDt(const label regionI) const;

        //- Processor that paces a main tissue
        static label owner(const label regionI)
        {
            return regionI % Pstream::nProcs();
        }

        //- prePacingProperties settings of a main tissue
        const prePacingIO::PrePacingConfig& config(const label regionI) const
        {
            return configs_[regionI];
        }

        //- A configured single cell of a main tissue with its pacing
        //  protocol; with fromPacedState, starting from its converged state
        autoPtr<ionicModel> newSingleCell
        (
            const label regionI,
            const bool fromPacedState
        ) const;

        //- Pace every enabled main tissue's single cell to convergence, each
        //  on its owner processor, and share the states with all processors
        void paceIonic();

        //- Seed tissueModel's states and Vm from the paced main tissues
        void seedIonic(ionicModel& tissueModel, volScalarField& Vm) const;

        virtual bool writeData(Ostream&) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
