#ifndef notsureLoopHandler_H
#define notsureLoopHandler_H

#include "fvCFD.H"
#include "ionicModel.H"
#include "pimpleControl.H"

class notsureLoopHandler
{
private:
    const fvMesh& mesh_;
    ionicModel& ionicModel_;

public:
    explicit notsureLoopHandler(const fvMesh& mesh, ionicModel& model);

    void notsureLoop
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
        const volTensorField& conductivity,
        pimpleControl& pimple
    );
};

#endif

