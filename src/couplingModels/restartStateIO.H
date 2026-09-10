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

\*---------------------------------------------------------------------------*/

#ifndef restartStateIO_H
#define restartStateIO_H

#include "OSspecific.H"
#include "fvMesh.H"
#include "PtrList.H"
#include "scalarField.H"
#include "error.H"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>

namespace Foam
{
namespace restartStateIO
{

static const std::uint32_t magic = 0x43524653;
static const std::uint32_t version = 1;

inline fileName path(const fvMesh& mesh, const word& name)
{
    return
        mesh.time().path()/mesh.time().timeName()/mesh.dbDir()/name;
}

template<class T>
inline void writeValue(std::ofstream& os, const T& value)
{
    os.write(reinterpret_cast<const char*>(&value), sizeof(T));
}


template<class T>
inline void readValue(std::ifstream& is, T& value, const fileName& statePath)
{
    is.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!is)
    {
        FatalErrorInFunction
            << "Unexpected end of restart state file " << statePath
            << exit(FatalError);
    }
}


inline void validateHeader
(
    const word& model,
    const label stride,
    const label nRows,
    std::ifstream& is,
    const fileName& statePath
)
{
    std::uint32_t fileMagic = 0;
    std::uint32_t fileVersion = 0;
    std::uint32_t fileScalarSize = 0;
    std::uint32_t modelSize = 0;
    std::uint64_t fileStride = 0;
    std::uint64_t fileRows = 0;
    readValue(is, fileMagic, statePath);
    readValue(is, fileVersion, statePath);
    readValue(is, fileScalarSize, statePath);
    readValue(is, modelSize, statePath);
    readValue(is, fileStride, statePath);
    readValue(is, fileRows, statePath);

    if (modelSize == 0 || modelSize > 128)
    {
        FatalErrorInFunction
            << "Invalid restart state file " << statePath
            << exit(FatalError);
    }

    std::string fileModel(modelSize, '\0');
    is.read(&fileModel[0], modelSize);
    if (!is)
    {
        FatalErrorInFunction
            << "Unexpected end of restart state file " << statePath
            << exit(FatalError);
    }

    if
    (
        fileMagic != magic
     || fileVersion != version
     || fileModel != model
     || fileScalarSize != sizeof(scalar)
     || fileStride != static_cast<std::uint64_t>(stride)
     || fileRows != static_cast<std::uint64_t>(nRows)
    )
    {
        FatalErrorInFunction
            << "Invalid restart state file " << statePath
            << exit(FatalError);
    }

}

inline void writeHeader
(
    const word& model,
    const label stride,
    const label nRows,
    std::ofstream& os
)
{
    const std::uint32_t scalarSize = sizeof(scalar);
    const std::uint32_t modelSize = model.size();
    const std::uint64_t fileStride = stride;
    const std::uint64_t fileRows = nRows;
    writeValue(os, magic);
    writeValue(os, version);
    writeValue(os, scalarSize);
    writeValue(os, modelSize);
    writeValue(os, fileStride);
    writeValue(os, fileRows);
    os.write(model.c_str(), modelSize);
    if (!os)
    {
        FatalErrorInFunction
            << "Cannot write restart state file for " << model
            << exit(FatalError);
    }
}


inline scalar readScalar(std::ifstream& is, const fileName& statePath)
{
    scalar value = 0.0;
    readValue(is, value, statePath);
    return value;
}


inline void writeScalar
(
    std::ofstream& os,
    const scalar value,
    const fileName& statePath
)
{
    writeValue(os, value);
    if (!os)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath
            << exit(FatalError);
    }
}

inline void checkValue(const scalar value, const fileName& statePath)
{
    if (!std::isfinite(static_cast<double>(value)))
    {
        FatalErrorInFunction
            << "Non-finite value in restart state file " << statePath
            << exit(FatalError);
    }
}

inline bool readStates
(
    const fvMesh& mesh,
    const word& model,
    const label stride,
    PtrList<scalarField>& states
)
{
    const fileName statePath = path(mesh, model + "State");
    if (!isFile(statePath))
    {
        return false;
    }
    std::ifstream is(statePath.c_str(), std::ios::binary);
    if (!is)
    {
        FatalErrorInFunction
            << "Cannot read restart state file " << statePath
            << exit(FatalError);
    }
    validateHeader(model, stride, states.size(), is, statePath);
    forAll(states, integrationPtI)
    {
        forAll(states[integrationPtI], stateI)
        {
            states[integrationPtI][stateI] = readScalar(is, statePath);
            checkValue(states[integrationPtI][stateI], statePath);
        }
    }
    return true;
}

inline void writeStates
(
    const fvMesh& mesh,
    const word& model,
    const label stride,
    const PtrList<scalarField>& states
)
{
    const fileName statePath = path(mesh, model + "State");
    std::ofstream os(statePath.c_str(), std::ios::binary | std::ios::trunc);
    if (!os)
    {
        FatalErrorInFunction
            << "Cannot write restart state file " << statePath
            << exit(FatalError);
    }
    writeHeader(model, stride, states.size(), os);
    forAll(states, integrationPtI)
    {
        forAll(states[integrationPtI], stateI)
        {
            checkValue(states[integrationPtI][stateI], statePath);
            writeScalar(os, states[integrationPtI][stateI], statePath);
        }
    }
}

}
}

#endif
