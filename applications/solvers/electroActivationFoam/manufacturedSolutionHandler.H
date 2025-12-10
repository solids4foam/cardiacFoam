#ifndef manufacturedSolutionHandler_H
#define manufacturedSolutionHandler_H

#include "fvCFD.H"
#include "ionicModel.H"
//#include "tmanufacturedFDA.H"

class manufacturedSolutionHandler
{
private:
    const fvMesh& mesh_;
    ionicModel& ionicModel_;

    scalar dx_;
    int dim_;
    int N_;

public:
    // Constructor
    manufacturedSolutionHandler(const fvMesh& mesh, ionicModel& model);

    void initializeManufactured
    (
        volScalarField& Vm,
        List<volScalarField*>& outFields,
        scalar dxstructured,
        int dim
    );


    // Manufactured-solution error post-processing
    void postProcess
    (
        const volScalarField& Vm,
        const List<volScalarField*>& outFields,
        const scalar dt,
        const int nsteps,
        const bool solveExplicit
    ) const;


};

#endif
