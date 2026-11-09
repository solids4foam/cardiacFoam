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

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#include "conductorLinearSolve.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "processorFvPatch.H"

#if defined(USE_PETSC) && !defined(USE_LIGHTWEIGHT_PHYSICSMODEL)

#include "foamPetscSnesHelper.H"
#include "petscErrorHandling.H"
#include "leadFieldAdjoint.H"
#include <cfenv>

namespace Foam
{

namespace
{

// Floating-point traps masked for the lifetime of the object, then the
// previous environment restored (flags raised meanwhile are discarded)
class untrappedFloatingPoint
{
    std::fenv_t saved_;

public:

    untrappedFloatingPoint()
    {
        std::feholdexcept(&saved_);
    }

    ~untrappedFloatingPoint()
    {
        std::feclearexcept(FE_ALL_EXCEPT);
        std::fesetenv(&saved_);
    }
};


// PETSc initialised with traps masked. A base of the solver listed before
// foamPetscSnesHelper, whose own PetscInitialize then finds PETSc running.
struct petscInitialisedUntrapped
{
    petscInitialisedUntrapped()
    {
        PetscBool initialised = PETSC_FALSE;
        PetscInitialized(&initialised);
        if (!initialised)
        {
            untrappedFloatingPoint guard;
            PetscInitialize(nullptr, nullptr, nullptr, nullptr);
        }
    }
};

} // End anonymous namespace


class conductorLinearSolve::impl
:
    private petscInitialisedUntrapped,
    public foamPetscSnesHelper
{
    fvScalarMatrix& M_;
    const volScalarField& psi_;
    const probedOperator& C_;
    const bool transpose_;

    // Right-hand side of the current solve
    scalarField b_;

    // Right-hand side vector of the solves after the first
    Vec bVec_;

    // The first solve has formed A and set up the KSP
    bool setUp_;

    label nCells() const
    {
        return M_.psi().mesh().nCells();
    }


public:

    impl
    (
        const word& solverEntry,
        fvScalarMatrix& M,
        const volScalarField& psi,
        const probedOperator& C,
        const direction dir
    )
    :
        petscInitialisedUntrapped(),
        foamPetscSnesHelper
        (
            solverEntry, fileName::null, M.psi().mesh(),
            foamPetscSnesHelper::CELLS, true, true
        ),
        M_(M),
        psi_(psi),
        C_(C),
        transpose_(dir == direction::transpose),
        b_(),
        bVec_(nullptr),
        setUp_(false)
    {
        // C's columns and the helper's rows must number cells alike
        if (C_.globalCells().offsets() != globalCells().offsets())
        {
            FatalErrorInFunction
                << "The probed operator and the PETSc matrix number the "
                << "cells differently." << exit(FatalError);
        }
    }

    ~impl()
    {
        if (bVec_)
        {
            VecDestroy(&bVec_);
        }
    }


    label initialiseJacobian(Mat& jac) override
    {
        const label n = nCells();
        const label start = globalCells().localStart();
        const label end = start + n;

        // Rows of C hold the cell and its face neighbours, so they bound
        // M's pattern; the transpose's rows can differ where the probing
        // pattern is not symmetric, which allocation on insertion absorbs
        List<PetscInt> dNnz(n, 1), oNnz(n, 0);
        const labelListList& columns = C_.columns();
        forAll(columns, cellI)
        {
            label d = 0, o = 0;
            for (const label col : columns[cellI])
            {
                (col >= start && col < end) ? ++d : ++o;
            }
            dNnz[cellI] = max(d, 1);
            oNnz[cellI] = o;
        }

        AssertPETSc(MatCreate(PETSC_COMM_WORLD, &jac));
        AssertPETSc(MatSetSizes(jac, n, n, PETSC_DETERMINE, PETSC_DETERMINE));
        AssertPETSc(MatSetType(jac, MATAIJ));
        AssertPETSc
        (
            MatXAIJSetPreallocation
            (
                jac, 1, dNnz.cdata(), oNnz.cdata(), nullptr, nullptr
            )
        );
        AssertPETSc
        (
            MatSetOption(jac, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE)
        );
        return 0;
    }

    label initialiseSolution(Vec& x) override
    {
        return foamPetscSnesHelper::initialiseSolution
        (
            x, M_.psi().mesh(), 1
        );
    }

    label formJacobian(Mat jac, const Vec) override
    {
        AssertPETSc(MatZeroEntries(jac));
        AssertPETSc(InsertFvMatrixIntoPETScMatrix(M_, jac, 0, 0, 1));

        const label start = globalCells().localStart();
        const labelListList& columns = C_.columns();
        const List<scalarList>& values = C_.values();
        forAll(columns, cellI)
        {
            const PetscInt row = start + cellI;
            forAll(columns[cellI], k)
            {
                const PetscInt col = columns[cellI][k];
                const PetscScalar v = values[cellI][k];
                if (transpose_)
                {
                    AssertPETSc(MatSetValue(jac, col, row, v, ADD_VALUES));
                }
                else
                {
                    AssertPETSc(MatSetValue(jac, row, col, v, ADD_VALUES));
                }
            }
        }

        AssertPETSc(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
        AssertPETSc(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
        return 0;
    }

    // f = A x - b with the stored operators: valid before A is formed
    label formResidual(Vec f, const Vec x) override
    {
        const label n = nCells();
        scalarField xf(n);
        {
            const PetscScalar* xa;
            AssertPETSc(VecGetArrayRead(x, &xa));
            forAll(xf, i)
            {
                xf[i] = xa[i];
            }
            AssertPETSc(VecRestoreArrayRead(x, &xa));
        }

        const scalarField r
        (
            leadFieldAdjoint::matrixProduct(M_, psi_, xf)
          + (transpose_ ? C_.applyTranspose(xf) : C_.apply(xf))
          - b_
        );

        PetscScalar* fa;
        AssertPETSc(VecGetArray(f, &fa));
        forAll(r, i)
        {
            fa[i] = r[i];
        }
        AssertPETSc(VecRestoreArray(f, &fa));
        return 0;
    }


    record solve(const scalarField& b, scalarField& x)
    {
        b_ = b;
        Vec xv = solution();
        AssertPETSc(VecSet(xv, 0.0));

        {
            untrappedFloatingPoint guard;
            if (!setUp_)
            {
                // Forms A, sets up the KSP and the preconditioner
                foamPetscSnesHelper::solve(false);
                AssertPETSc(VecDuplicate(xv, &bVec_));
                setUp_ = true;
            }
            else
            {
                PetscScalar* ba;
                AssertPETSc(VecGetArray(bVec_, &ba));
                forAll(b, i)
                {
                    ba[i] = b[i];
                }
                AssertPETSc(VecRestoreArray(bVec_, &ba));

                KSP ksp;
                AssertPETSc(SNESGetKSP(snes(), &ksp));
                AssertPETSc(KSPSolve(ksp, bVec_, xv));
            }
        }

        KSP ksp;
        AssertPETSc(SNESGetKSP(snes(), &ksp));
        PetscInt iterations = 0;
        AssertPETSc(KSPGetIterationNumber(ksp, &iterations));
        KSPConvergedReason reason;
        AssertPETSc(KSPGetConvergedReason(ksp, &reason));

        if (reason <= 0)
        {
            FatalErrorInFunction
                << "The conductor solve did not converge: KSP "
                << KSPConvergedReasons[reason] << " after " << iterations
                << " iterations." << exit(FatalError);
        }

        x.setSize(nCells());
        const PetscScalar* xa;
        AssertPETSc(VecGetArrayRead(xv, &xa));
        forAll(x, i)
        {
            x[i] = xa[i];
        }
        AssertPETSc(VecRestoreArrayRead(xv, &xa));

        return record{label(iterations), word(KSPConvergedReasons[reason])};
    }
};

} // End namespace Foam

#else

namespace Foam
{
class conductorLinearSolve::impl {};
}

#endif


namespace Foam
{

conductorLinearSolve::conductorLinearSolve
(
    const word& solverEntry,
    fvScalarMatrix& M,
    const volScalarField& psi,
    const probedOperator& C,
    const direction dir
)
:
    entry_(M.psi().mesh().solverDict(solverEntry)),
    implPtr_()
{
    const word solver(entry_.get<word>("solver"));
    if (solver != "petsc")
    {
        FatalIOErrorInFunction(entry_)
            << "The conductor solve needs 'solver petsc;' and PETSc "
            << "'options' in the fvSolution entry " << solverEntry
            << "; it has solver " << solver << '.' << exit(FatalIOError);
    }

    for (const fvPatch& patch : M.psi().mesh().boundary())
    {
        if (patch.coupled() && !isA<processorFvPatch>(patch))
        {
            FatalErrorInFunction
                << "The conductor solve does not support the coupled patch "
                << patch.name() << " (" << patch.type() << "): a heart in a "
                << "torso has no periodic boundary." << exit(FatalError);
        }
    }

#if defined(USE_PETSC) && !defined(USE_LIGHTWEIGHT_PHYSICSMODEL)
    implPtr_.reset(new impl(solverEntry, M, psi, C, dir));
#else
    FatalErrorInFunction
        << "The conductor solve needs cardiacFoam built with PETSc and the "
        << "full solids4foam (PETSC_DIR set, not the lightweight "
        << "physicsModel)." << exit(FatalError);
#endif
}


conductorLinearSolve::~conductorLinearSolve()
{}


conductorLinearSolve::record conductorLinearSolve::solve
(
    const scalarField& b,
    scalarField& x
)
{
#if defined(USE_PETSC) && !defined(USE_LIGHTWEIGHT_PHYSICSMODEL)
    return implPtr_->solve(b, x);
#else
    FatalErrorInFunction
        << "Built without PETSc." << exit(FatalError);
    return record{0, word::null};
#endif
}

} // End namespace Foam

// ************************************************************************* //
