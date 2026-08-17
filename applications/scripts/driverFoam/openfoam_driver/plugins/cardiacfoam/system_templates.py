#----------------------------------------------------------------------------#
# License
#     This file is part of cardiacFoam.
#
#     cardiacFoam is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     cardiacFoam is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with cardiacFoam.  If not, see <http://www.gnu.org/licenses/>.
#
# Module
#     plugins.cardiacfoam.system_templates
#
# Description
#     Baseline cardiacFoam system dictionaries, one set per myocardiumSolver.
#
# Author
#     Simao Nieto de Castro, UCD.
#----------------------------------------------------------------------------#

"""Baseline templates for system dictionaries (fvSolution, fvSchemes, controlDict).

These templates provide a clean starting point for generating cases from scratch.
Overrides can be applied to them using standard driverFoam mechanisms.

Every template here is cardiacFoam's: the controlDict names `cardiacFoam` as
its application, and the schemes/solution sets are keyed on the fields each
myocardiumSolver actually solves for (`Vm`, `psi`, `phiE`/`phiI`). That is why
this module lives in the plugin rather than in `specs/` -- there is no
solver-neutral fvSchemes to share.
"""

from __future__ import annotations


_FOAM_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1912                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      {object_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""


# -----------------------------------------------------------------------------
# controlDict templates
# -----------------------------------------------------------------------------

_CONTROL_DICT_BASE = """
application     cardiacFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time};
deltaT          {delta_t};

writeControl    runTime;
writeInterval   {write_interval};
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;

timeFormat      general;
timePrecision   6;
runTimeModifiable false;
"""

def build_control_dict(delta_t: float = 1e-4, end_time: float = 1.0, write_interval: float = 0.1) -> str:
    body = _CONTROL_DICT_BASE.format(
        delta_t=delta_t,
        end_time=end_time,
        write_interval=write_interval
    )
    return _FOAM_HEADER.format(object_name="controlDict") + body + "\n// ************************************************************************* //\n"


# -----------------------------------------------------------------------------
# fvSchemes templates
# -----------------------------------------------------------------------------

FV_SCHEMES_SINGLE_CELL = _FOAM_HEADER.format(object_name="fvSchemes") + """
ddtSchemes
{
    default         backward; // 2nd order
}

gradSchemes
{
    default         leastSquares;
}

divSchemes
{
    default         none;
    div(phiU,psi)   Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""

FV_SCHEMES_MONODOMAIN = _FOAM_HEADER.format(object_name="fvSchemes") + """
ddtSchemes
{
    default         none;
    ddt(Vm)         backward;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""

FV_SCHEMES_BIDOMAIN = _FOAM_HEADER.format(object_name="fvSchemes") + """
ddtSchemes
{
    default         none;
    ddt(Vm)         backward;
}

gradSchemes
{
    default         Gauss linear;
}

divSchemes
{
    default         none;
    div((conductivityIntracellular&grad(Vm))) Gauss linear;
    div((conductivityIntracellular&grad(phiE))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

// ************************************************************************* //
"""


# -----------------------------------------------------------------------------
# fvSolution templates
# -----------------------------------------------------------------------------

FV_SOLUTION_SINGLE_CELL = _FOAM_HEADER.format(object_name="fvSolution") + """
solvers
{
    "Vm|VmFinal|u|uFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-11;
        relTol          0.0;
    }

    "psi|psiFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-11;
        relTol          0.1;
    }
}

PIMPLE
{
    nOuterCorrectors    1000;

    residualControl
    {
        "Vm|psi"
        {
            relTol          0.0;
            tolerance       1e-6;
        }
    }
}

relaxationFactors
{
    equations
    {
        psi    0.9;
    }
}

// ************************************************************************* //
"""

FV_SOLUTION_MONODOMAIN = _FOAM_HEADER.format(object_name="fvSolution") + """
solvers
{
    "Vm|VmFinal|u|uFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-15;
        relTol          0.0;
    }
}

PIMPLE
{
    nOuterCorrectors     2;
}

// ************************************************************************* //
"""

FV_SOLUTION_BIDOMAIN = _FOAM_HEADER.format(object_name="fvSolution") + """
solvers
{
    "Vm|VmFinal|u|uFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-15;
        relTol          0.0;
    }

    "phiE|phiEFinal|phiI|phiIFinal"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-15;
        relTol          0.0;
    }
}

PIMPLE
{
    nOuterCorrectors     2;
}

// ************************************************************************* //
"""

def get_fv_schemes(solver_type: str) -> str:
    if "singleCell" in solver_type:
        return FV_SCHEMES_SINGLE_CELL
    elif "bidomain" in solver_type.lower():
        return FV_SCHEMES_BIDOMAIN
    else:
        return FV_SCHEMES_MONODOMAIN

def get_fv_solution(solver_type: str) -> str:
    if "singleCell" in solver_type:
        return FV_SOLUTION_SINGLE_CELL
    elif "bidomain" in solver_type.lower():
        return FV_SOLUTION_BIDOMAIN
    else:
        return FV_SOLUTION_MONODOMAIN
