#ifndef explicitLoopHandler_H
#define explicitLoopHandler_H

#include "fvCFD.H"
#include "ionicModel.H"

class explicitLoopHandler
{
private:
    const fvMesh& mesh_;
    ionicModel& ionicModel_;

    // Cached solver parameters
    scalar dx_;
    int dim_;
    scalar Def_max_;

    // Internal helper: compute dt stability limit
    scalar computeStableDt(const scalar CFL) const;

public:
    explicitLoopHandler(const fvMesh& mesh, ionicModel& model);

    // Precompute dx, dim, Def_max and update dt
    void initializeExplicit
    (
        scalar& dt,
        int& nsteps,
        const scalar chiVal,
        const scalar CmVal,
        const volTensorField& conductivity,
        const scalar CFL,
        const scalar dx,
        const int dim

    );

    // Perform one explicit time step
    void explicitLoop
    (
        const scalar t0,
        const scalar dt,
        volScalarField& Vm,
        volScalarField& Iion,
        Field<Field<scalar>>& states,
        volScalarField& externalStimulusCurrent,
        const labelList& stimulusCellIDs,
        const scalar stimulusIntensity,
        const scalar stimulusDuration,
        const dimensionedScalar& chi,
        const dimensionedScalar& Cm,
        const volTensorField& conductivity
    );
};

#endif
