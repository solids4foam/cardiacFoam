#ifndef manufacturedSolutionHandler_H
#define manufacturedSolutionHandler_H

#include "fvCFD.H"
#include "ionicModel.H"
#include "tmanufacturedFDA.H"

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

    // Manufactured-solution analytic initialization
    void initializeManufactured
    (
        volScalarField& Vm,
        volScalarField& u1,
        volScalarField& u2,
        volScalarField& u3
    );

    // Non-manufactured initialization
    void initializeNonManufactured();

    // Manufactured-solution error post-processing
    void postProcess
    (
        const volScalarField& Vm,
        const volScalarField& u1,
        const volScalarField& u2,
        const volScalarField& u3,
        const scalar dt,
        const int nsteps,
        const bool solveExplicit
    ) const;

    scalar dx() const { return dx_; }
    int dim() const { return dim_; }
    int N() const { return N_; }
};

#endif
