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
    conductorLinearSolve

Description
    Linear solve with the lead-field conductor's operator A = M + C
    (forward) or its transpose M + C^T (lead fields), through solids4foam's
    PETSc layer (foamPetscSnesHelper).

    M is the implicit fvScalarMatrix of passiveConductor::assemble with its
    reference pin, C the explicit part stored by leadFieldAdjoint
    (probedOperator). The PETSc matrix is M, inserted by the helper with the
    diagonal of every patch, plus the probed entries of C (row and column
    swapped for the transpose). It is formed once; every right-hand side
    reuses it and its preconditioner.

    The solver entry of fvSolution (named at construction) needs
        solver   petsc;
        options  { ... PETSc options, snes_type ksponly ... }
    The first solve goes through the helper's SNES (one Newton step of the
    linear residual, ksponly), which forms A and sets up the KSP and the
    preconditioner from the options; later solves call KSPSolve on that KSP.
    A KSP that does not converge is fatal. The true residual of the real
    operator is the caller's check (leadFieldAdjoint::checkSolution, the
    lead-field verifier's direct solve).

    PETSc runs with floating-point traps masked: OpenFOAM's trap
    (FOAM_SIGFPE) stops PETSc's start-up and its preconditioners, which
    raise floating-point exceptions they expect to run untrapped. The traps
    are restored after each PETSc call.

    Cyclic patches are not supported (the helper inserts processor
    interfaces only); a heart in a torso has no periodic boundary.

    The class exists in builds with PETSc and the full solids4foam; any other
    build stops at construction.

SourceFiles
    conductorLinearSolve.C

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef conductorLinearSolve_H
#define conductorLinearSolve_H

#include "fvMatricesFwd.H"
#include "volFieldsFwd.H"
#include "dictionary.H"
#include "probedOperator.H"
#include <memory>

namespace Foam
{

class conductorLinearSolve
{
public:

    //- The system solved: A x = b or A^T x = b
    enum class direction { forward, transpose };

    //- How one solve went
    struct record
    {
        //- KSP iterations
        label iterations;

        //- KSP converged reason
        word reason;
    };


private:

    class impl;

    // The solver entry of fvSolution
    const dictionary entry_;

    std::unique_ptr<impl> implPtr_;


public:

    //- Solver for M + C (forward) or M + C^T (transpose). M is the
    //  implicit matrix with its reference pin and psi its unknown (whose
    //  interfaces M uses); C the probed explicit part on the same mesh.
    //  M, psi and C must outlive the solver.
    conductorLinearSolve
    (
        const word& solverEntry,
        fvScalarMatrix& M,
        const volScalarField& psi,
        const probedOperator& C,
        const direction dir
    );

    conductorLinearSolve(const conductorLinearSolve&) = delete;
    void operator=(const conductorLinearSolve&) = delete;

    ~conductorLinearSolve();


    //- x = A^-1 b (A^-T b for the transpose), from x = 0; fatal unless the
    //  KSP converged
    record solve(const scalarField& b, scalarField& x);

    //- The solver entry of fvSolution
    const dictionary& entry() const
    {
        return entry_;
    }
};

} // End namespace Foam

#endif

// ************************************************************************* //
