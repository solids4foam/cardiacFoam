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

Application
    1DgraphToFoam

Description
    Convert a legacy ASCII VTK 1D graph stored as UNSTRUCTURED_GRID/POLYDATA
    lines into a Foam dictionary while preserving point and line metadata.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"
#include "vtkUnstructuredReader.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "vectorIOField.H"

using namespace Foam;

namespace
{

template<class Type>
void writeMappedFields
(
    Ostream& os,
    const objectRegistry& registry,
    const word& dictName,
    const labelList& map,
    const label expectedSize
)
{
    const UPtrList<const IOField<Type>> fields
    (
        registry.csorted<IOField<Type>>()
    );

    os  << dictName << nl
        << token::BEGIN_BLOCK << incrIndent << nl;

    for (const IOField<Type>& field : fields)
    {
        if (field.empty())
        {
            continue;
        }

        os  << indent << field.name() << nl
            << indent << token::BEGIN_LIST << incrIndent << nl;

        forAll(map, i)
        {
            const label sourceI = map[i];

            if (sourceI < 0 || sourceI >= field.size())
            {
                FatalErrorInFunction
                    << "Cannot map field '" << field.name()
                    << "': source index " << sourceI
                    << " is outside field size " << field.size() << "."
                    << exit(FatalError);
            }

            os << indent << field[sourceI] << nl;
        }

        os  << decrIndent << indent << token::END_LIST
            << token::END_STATEMENT << nl;
    }

    if (expectedSize >= 0)
    {
        for (const IOField<Type>& field : fields)
        {
            if (field.size() != expectedSize)
            {
                Info<< "Note: field '" << field.name() << "' has "
                    << field.size() << " values; expected " << expectedSize
                    << " before mapping." << nl;
            }
        }
    }

    os  << decrIndent << token::END_BLOCK << nl;
}


template<class Type>
void writeDirectFields
(
    Ostream& os,
    const objectRegistry& registry,
    const word& dictName,
    const label expectedSize
)
{
    const UPtrList<const IOField<Type>> fields
    (
        registry.csorted<IOField<Type>>()
    );

    os  << dictName << nl
        << token::BEGIN_BLOCK << incrIndent << nl;

    for (const IOField<Type>& field : fields)
    {
        if (field.size() != expectedSize)
        {
            FatalErrorInFunction
                << "Cannot write field '" << field.name()
                << "': expected " << expectedSize
                << " values but found " << field.size() << "."
                << exit(FatalError);
        }

        os  << indent << field.name() << nl
            << indent << field << token::END_STATEMENT << nl;
    }

    os  << decrIndent << token::END_BLOCK << nl;
}


labelList endpointNodes(const label nPoints, const labelListList& lines)
{
    labelList degree(nPoints, 0);

    forAll(lines, lineI)
    {
        const labelList& line = lines[lineI];

        forAll(line, i)
        {
            const label pointI = line[i];
            if (pointI < 0 || pointI >= nPoints)
            {
                FatalErrorInFunction
                    << "Line " << lineI << " uses point " << pointI
                    << " outside point count " << nPoints << "."
                    << exit(FatalError);
            }
        }

        for (label i = 1; i < line.size(); ++i)
        {
            ++degree[line[i - 1]];
            ++degree[line[i]];
        }
    }

    DynamicList<label> endpoints;
    forAll(degree, pointI)
    {
        if (degree[pointI] == 1)
        {
            endpoints.append(pointI);
        }
    }

    return endpoints.shrink();
}


scalarField lineLengths(const pointField& points, const labelListList& lines)
{
    scalarField lengths(lines.size(), 0.0);

    forAll(lines, lineI)
    {
        const labelList& line = lines[lineI];
        for (label i = 1; i < line.size(); ++i)
        {
            lengths[lineI] += mag(points[line[i]] - points[line[i - 1]]);
        }
    }

    return lengths;
}


template<class Type>
const IOField<Type>* findField
(
    const objectRegistry& registry,
    const wordList& names
)
{
    for (const word& name : names)
    {
        if (const IOField<Type>* field = registry.cfindObject<IOField<Type>>(name))
        {
            return field;
        }
    }

    return nullptr;
}


label roleAt
(
    const label nodeI,
    const IOField<label>* labelRoles,
    const IOField<scalar>* scalarRoles
)
{
    if (labelRoles)
    {
        return (*labelRoles)[nodeI];
    }

    if (scalarRoles)
    {
        return label(round((*scalarRoles)[nodeI]));
    }

    return 0;
}


label rootNodeFromFields
(
    const label nPoints,
    const IOField<label>* labelRoles,
    const IOField<scalar>* scalarRoles
)
{
    label rootNode = -1;

    for (label nodeI = 0; nodeI < nPoints; ++nodeI)
    {
        if (roleAt(nodeI, labelRoles, scalarRoles) == 1)
        {
            if (rootNode != -1)
            {
                FatalErrorInFunction
                    << "nodeRole marks more than one root node. Found "
                    << rootNode << " and " << nodeI << "."
                    << exit(FatalError);
            }
            rootNode = nodeI;
        }
    }

    if (rootNode == -1)
    {
        rootNode = 0;
        Info<< "No nodeRole root marker was found. Using node 0 as root."
            << nl;
    }

    return rootNode;
}


scalar edgeConductanceFromFields
(
    const label lineI,
    const objectRegistry& cellData,
    const labelList& lineMap
)
{
    static const wordList conductanceNames
    {
        "conductance",
        "diffusivity",
        "D",
        "sigma",
        "conductivity"
    };

    const IOField<scalar>* conductance =
        findField<scalar>(cellData, conductanceNames);

    if (!conductance)
    {
        return 1.0;
    }

    const label sourceI = lineMap[lineI];
    if (sourceI < 0 || sourceI >= conductance->size())
    {
        FatalErrorInFunction
            << "Cannot map edge conductance from field '"
            << conductance->name() << "': lineMap[" << lineI << "]="
            << sourceI << " is outside field size " << conductance->size()
            << "."
            << exit(FatalError);
    }

    return (*conductance)[sourceI];
}


pointField writeConductionSolverContract
(
    Ostream& os,
    const pointField& points,
    const labelListList& lines,
    const objectRegistry& pointData,
    const objectRegistry& cellData,
    const labelList& lineMap,
    const scalar graphStep
)
{
    if (lineMap.size() != lines.size())
    {
        FatalErrorInFunction
            << "lineMap has " << lineMap.size() << " entries but lines has "
            << lines.size() << "."
            << exit(FatalError);
    }

    const wordList nodeRoleNames{"nodeRole", "node_role", "role"};
    const IOField<label>* labelRoles =
        findField<label>(pointData, nodeRoleNames);
    const IOField<scalar>* scalarRoles =
        findField<scalar>(pointData, nodeRoleNames);

    const label rootNode =
        rootNodeFromFields(points.size(), labelRoles, scalarRoles);

    const labelList endpoints = endpointNodes(points.size(), lines);

    DynamicList<label> pvjNodes;
    if (labelRoles || scalarRoles)
    {
        forAll(points, nodeI)
        {
            if (roleAt(nodeI, labelRoles, scalarRoles) == 2)
            {
                pvjNodes.append(nodeI);
            }
        }
    }

    if (pvjNodes.empty())
    {
        forAll(endpoints, i)
        {
            if (endpoints[i] != rootNode)
            {
                pvjNodes.append(endpoints[i]);
            }
        }
    }

    pointField pvjLocations(pvjNodes.size());
    forAll(pvjNodes, i)
    {
        pvjLocations[i] = points[pvjNodes[i]];
    }

    // Expanded point list: original points + interpolated intermediate nodes.
    DynamicList<point> expandedPoints(points);
    label nextNodeId = points.size();

    DynamicList<scalarList> conductionEdges;
    forAll(lines, lineI)
    {
        const labelList& line = lines[lineI];

        if (line.size() < 2)
        {
            FatalErrorInFunction
                << "Line " << lineI << " has " << line.size()
                << " point(s). At least two are required."
                << exit(FatalError);
        }

        const scalar cond =
            edgeConductanceFromFields(lineI, cellData, lineMap);

        for (label i = 1; i < line.size(); ++i)
        {
            const point& pA = points[line[i - 1]];
            const point& pB = points[line[i]];
            const scalar len = mag(pB - pA);

            if (graphStep > 0 && len > graphStep)
            {
                const label nSeg =
                    max(1, label(std::ceil(len / graphStep)));
                const scalar subLen = len / nSeg;

                label prevNode = line[i - 1];
                for (label s = 1; s < nSeg; ++s)
                {
                    const scalar alpha = scalar(s) / scalar(nSeg);
                    expandedPoints.append((1 - alpha)*pA + alpha*pB);
                    const label curNode = nextNodeId++;
                    scalarList edge(4, 0.0);
                    edge[0] = prevNode;
                    edge[1] = curNode;
                    edge[2] = subLen;
                    edge[3] = cond;
                    conductionEdges.append(edge);
                    prevNode = curNode;
                }
                scalarList edge(4, 0.0);
                edge[0] = prevNode;
                edge[1] = line[i];
                edge[2] = subLen;
                edge[3] = cond;
                conductionEdges.append(edge);
            }
            else
            {
                scalarList edge(4, 0.0);
                edge[0] = line[i - 1];
                edge[1] = line[i];
                edge[2] = len;
                edge[3] = cond;
                conductionEdges.append(edge);
            }
        }
    }

    os  << "rootNode" << nl
        << rootNode << token::END_STATEMENT << nl << nl
        << "pvjNodes" << nl
        << labelList(pvjNodes.shrink()) << token::END_STATEMENT << nl << nl
        << "pvjLocations" << nl
        << pvjLocations << token::END_STATEMENT << nl << nl
        << "conductionEdges" << nl
        << List<scalarList>(conductionEdges.shrink())
        << token::END_STATEMENT << nl << nl;

    const label nSubdivided = nextNodeId - points.size();
    Info<< "1DgraphToFoam solver contract:" << nl
        << "  rootNode: " << rootNode << nl
        << "  pvjNodes: " << pvjNodes.size() << nl
        << "  conductionEdges: " << conductionEdges.size() << nl;
    if (nSubdivided > 0)
    {
        Info<< "  subdivided edges (graphStep=" << graphStep << " m): "
            << nSubdivided << " intermediate nodes inserted" << nl;
    }
    Info<< endl;

    return pointField(expandedPoints.shrink());
}

} // End anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert a legacy ASCII VTK 1D graph into constant/<name> while "
        "preserving POINT_DATA and line CELL_DATA."
    );

    argList::noParallel();
    argList::addArgument("vtk-file", "The input legacy ASCII VTK file");
    argList::addOption
    (
        "name",
        "word",
        "Name of the graph object written under constant/ (default: purkinjeGraph)"
    );
    argList::addOption
    (
        "graphStep",
        "scalar",
        "Maximum edge length after subdivision (m). "
        "Edges longer than this are split into equal sub-segments. "
        "Default 3e-4 m (0.3 mm). Set to 0 to disable."
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const fileName vtkFile(args.get<fileName>(1));
    const word graphName(args.getOrDefault<word>("name", "purkinjeGraph"));
    const scalar graphStep(args.getOrDefault<scalar>("graphStep", 3e-4));

    IFstream vtkStream(vtkFile);
    vtkUnstructuredReader reader(runTime, vtkStream);

    const pointField& points = reader.points();
    const labelListList& lines = reader.lines();

    if (points.empty())
    {
        FatalErrorInFunction
            << "No VTK POINTS were found in '" << vtkFile << "'."
            << exit(FatalError);
    }

    if (lines.empty())
    {
        FatalErrorInFunction
            << "No VTK line elements were found in '" << vtkFile
            << "'. Use VTK_LINE, VTK_POLY_LINE, or POLYDATA LINES."
            << exit(FatalError);
    }

    const fileName outputFile(runTime.path()/runTime.constant()/graphName);

    Info<< "Writing 1D graph to " << outputFile << nl
        << "  points: " << points.size() << nl
        << "  lines:  " << lines.size() << nl << endl;

    OFstream os(outputFile);

    os  << "FoamFile" << nl
        << token::BEGIN_BLOCK << incrIndent << nl
        << indent << "version     2.0;" << nl
        << indent << "format      ascii;" << nl
        << indent << "class       dictionary;" << nl
        << indent << "location    \"constant\";" << nl
        << indent << "object      " << graphName << token::END_STATEMENT << nl
        << decrIndent << token::END_BLOCK << nl << nl;

    const pointField expandedPoints
    (
        writeConductionSolverContract
        (
            os,
            points,
            lines,
            reader.pointData(),
            reader.cellData(),
            reader.lineMap(),
            graphStep
        )
    );

    os  << "points" << nl << expandedPoints << token::END_STATEMENT << nl << nl
        << "edges" << nl << lines << token::END_STATEMENT << nl << nl
        << "edgeVtkCellMap" << nl << reader.lineMap()
        << token::END_STATEMENT << nl << nl
        << "edgeLength" << nl << lineLengths(points, lines)
        << token::END_STATEMENT << nl << nl
        << "endpointNodes" << nl << endpointNodes(points.size(), lines)
        << token::END_STATEMENT << nl << nl;

    os  << "pointFields" << nl
        << token::BEGIN_BLOCK << incrIndent << nl;
    writeDirectFields<label>(os, reader.pointData(), "label", points.size());
    writeDirectFields<scalar>(os, reader.pointData(), "scalar", points.size());
    writeDirectFields<vector>(os, reader.pointData(), "vector", points.size());
    writeDirectFields<sphericalTensor>
    (
        os,
        reader.pointData(),
        "sphericalTensor",
        points.size()
    );
    writeDirectFields<symmTensor>
    (
        os,
        reader.pointData(),
        "symmTensor",
        points.size()
    );
    writeDirectFields<tensor>(os, reader.pointData(), "tensor", points.size());
    os  << decrIndent << token::END_BLOCK
        << token::END_STATEMENT << nl << nl;

    os  << "edgeFields" << nl
        << token::BEGIN_BLOCK << incrIndent << nl;
    writeMappedFields<label>
    (
        os,
        reader.cellData(),
        "label",
        reader.lineMap(),
        reader.cells().size() + reader.faces().size() + reader.lines().size()
    );
    writeMappedFields<scalar>
    (
        os,
        reader.cellData(),
        "scalar",
        reader.lineMap(),
        reader.cells().size() + reader.faces().size() + reader.lines().size()
    );
    writeMappedFields<vector>
    (
        os,
        reader.cellData(),
        "vector",
        reader.lineMap(),
        reader.cells().size() + reader.faces().size() + reader.lines().size()
    );
    writeMappedFields<sphericalTensor>
    (
        os,
        reader.cellData(),
        "sphericalTensor",
        reader.lineMap(),
        reader.cells().size() + reader.faces().size() + reader.lines().size()
    );
    writeMappedFields<symmTensor>
    (
        os,
        reader.cellData(),
        "symmTensor",
        reader.lineMap(),
        reader.cells().size() + reader.faces().size() + reader.lines().size()
    );
    writeMappedFields<tensor>
    (
        os,
        reader.cellData(),
        "tensor",
        reader.lineMap(),
        reader.cells().size() + reader.faces().size() + reader.lines().size()
    );
    os  << decrIndent << token::END_BLOCK
        << token::END_STATEMENT << nl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
