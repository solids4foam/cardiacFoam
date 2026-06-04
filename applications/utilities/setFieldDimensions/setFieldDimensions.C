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
    setFieldDimensions

Description
    Assign correct SI dimensions to OpenFOAM field files that were produced
    by newVtkUnstructuredToFoam, which writes every CELL_DATA field as
    dimensionless [0 0 0 0 0 0 0] because the VTK format carries no unit
    information.

    The utility patches only the header line that begins with "dimensions"
    and leaves the rest of the field file (including internalField and
    boundaryField) entirely unchanged. This makes it safe and fast even on
    large meshes with millions of cells.

    A built-in catalogue covers the known cardiacFoam fields that need
    correction after VTK import. Additional fields can be targeted with
    -field and -dim.

Usage
    setFieldDimensions [OPTIONS]

Options
    -field <name>
        Process a single named field. If omitted, all fields in the built-in
        catalogue that exist in the time directory are processed.

    -dim "[M L T K mol A cd]"
        Dimension set to apply. Required when -field names a field not in the
        built-in catalogue. Must be quoted and space-separated, e.g.
            -dim "[-1 -3 3 0 0 2 0]"

    -time <t>
        Time directory to process. Defaults to "0".

    -dryRun
        Print proposed changes without writing any file.

Examples
    # Fix all known fields in 0/ (typical post-VTK-import step)
    setFieldDimensions -case ./myCase

    # Fix a specific field
    setFieldDimensions -field Diffusivity -case ./myCase

    # Fix a custom field with explicit dimensions
    setFieldDimensions -field myPressure -dim "[1 -1 -2 0 0 0 0]" -case ./myCase

    # Preview without writing
    setFieldDimensions -dryRun -case ./myCase

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dimensionSet.H"
#include "OStringStream.H"
#include "IStringStream.H"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace Foam;

// ─────────────────────────────────────────────────────────────────────────────
// Built-in catalogue: field name → correct SI dimensionSet
// Only fields that need correction after VTK import are listed here.
// Genuinely dimensionless fields (fiber, sheet, uvc_transmural, tags) are
// intentionally absent — they are already correct.
// ─────────────────────────────────────────────────────────────────────────────

struct FieldDimEntry
{
    const char* name;
    dimensionSet dims;
    const char* note;
};

static const FieldDimEntry catalogue[] =
{
    // Conductivity tensor / scalar  [S/m] = kg^-1 m^-3 s^3 A^2
    { "Diffusivity",               dimensionSet(-1, -3, 3, 0, 0, 2, 0),
      "conductivity tensor from VTK import (S/m)" },
    { "conductivity",              dimensionSet(-1, -3, 3, 0, 0, 2, 0),
      "conductivity tensor from VTK import (S/m)" },
    { "bodyAndOrgansConductivity", dimensionSet(-1, -3, 3, 0, 0, 2, 0),
      "torso/organ conductivity scalar for ECG bath domain (S/m)" },
};
static const int catalogueSize = sizeof(catalogue) / sizeof(catalogue[0]);


// ─────────────────────────────────────────────────────────────────────────────
// Patch the "dimensions" line in an OpenFOAM field file.
// Returns true if the file was changed (or would be changed in dry-run).
// ─────────────────────────────────────────────────────────────────────────────

bool patchDimensions
(
    const fileName& filePath,
    const dimensionSet& newDims,
    bool dryRun
)
{
    // Read entire file into memory — we must not parse the (potentially huge)
    // internalField token list; we only need to touch one header line.
    std::ifstream in(filePath.c_str(), std::ios::binary);
    if (!in.good())
    {
        WarningInFunction
            << "Cannot open " << filePath << " for reading" << endl;
        return false;
    }

    std::string content
    (
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>()
    );
    in.close();

    // Find the "dimensions" keyword line (not inside a comment).
    // We search for "\ndimensions" (or "^dimensions" at file start) followed
    // by whitespace and the bracket-delimited token.
    const std::string keyword = "dimensions";
    std::size_t pos = std::string::npos;

    // Search from the top — dimensions is always near the file header.
    // A valid match must be at the start of a line (preceded by \n or at
    // position 0). Lines like "// dimensions" are skipped automatically
    // because the "//" means the preceding character is not \n.
    std::size_t search = 0;
    while (search < content.size())
    {
        std::size_t found = content.find(keyword, search);
        if (found == std::string::npos) break;

        // Must be at line start: position 0 or directly after a newline.
        if (found > 0 && content[found - 1] != '\n')
        {
            search = found + 1;
            continue;
        }

        pos = found;
        break;
    }

    if (pos == std::string::npos)
    {
        WarningInFunction
            << "No 'dimensions' entry found in " << filePath << endl;
        return false;
    }

    // Find the end of this logical entry (the semicolon)
    std::size_t semi = content.find(';', pos);
    if (semi == std::string::npos)
    {
        WarningInFunction
            << "Malformed dimensions entry in " << filePath << endl;
        return false;
    }

    // Build the current dimension string from the file for reporting
    std::string oldEntry = content.substr(pos, semi - pos + 1);

    // Build the replacement line
    OStringStream oss;
    oss << newDims;
    const std::string newEntry =
        "dimensions      " + std::string(oss.str().c_str()) + ";";

    // Check if already correct (avoid unnecessary writes).
    // oss.str() is the formatted new dimension (e.g. "[ -1 -3  3  0  0  2  0 ]").
    // If oldEntry already contains that string, the file is already correct.
    if (oldEntry.find(oss.str().c_str()) != std::string::npos)
    {
        Info << "  ok    " << filePath.name()
             << " (dimensions already correct)" << nl;
        return false;
    }

    // Report
    Info << "  " << (dryRun ? "would fix" : "fixing") << "  "
         << filePath.name() << nl
         << "          old: " << oldEntry << nl
         << "          new: " << newEntry << nl;

    if (dryRun) return true;

    // Splice the replacement into the content
    const std::string patched =
        content.substr(0, pos)
      + newEntry
      + content.substr(semi + 1);

    std::ofstream out(filePath.c_str(), std::ios::binary);
    if (!out.good())
    {
        WarningInFunction
            << "Cannot open " << filePath << " for writing" << endl;
        return false;
    }
    out.write(patched.c_str(), static_cast<std::streamsize>(patched.size()));
    return true;
}


// ─────────────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    argList::noParallel();
    argList::addNote
    (
        "Assign correct SI dimensions to dimensionless fields produced by\n"
        "newVtkUnstructuredToFoam. Patches only the 'dimensions' header line;\n"
        "internalField and boundaryField data are never parsed or rewritten."
    );

    argList::addOption
    (
        "field",
        "word",
        "process only this named field (default: all fields in the built-in catalogue)"
    );
    argList::addOption
    (
        "dim",
        "string",
        "dimension set to apply, e.g. \"[-1 -3 3 0 0 2 0]\" — required when\n"
        "-field names a field not in the built-in catalogue"
    );
    argList::addOption
    (
        "time",
        "word",
        "time directory to process (default: 0)"
    );
    argList::addBoolOption
    (
        "dryRun",
        "print proposed changes without writing any file"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool dryRun = args.found("dryRun");
    const word timeName = args.opt<word>("time", "0");

    Info<< "\n========== setFieldDimensions ==========\n" << endl;
    if (dryRun) Info<< "  (dry-run mode — no files will be written)\n" << endl;

    const fileName timeDir = runTime.path() / timeName;

    if (!isDir(timeDir))
    {
        FatalErrorInFunction
            << "Time directory not found: " << timeDir << nl
            << "Use -time to specify a different directory."
            << exit(FatalError);
    }

    if (args.found("field"))
    {
        // ── Single-field mode ─────────────────────────────────────────────
        const word fieldName = args["field"];

        // Resolve dimensions: user-supplied or catalogue look-up
        dimensionSet newDims(dimless);
        bool found = false;

        if (args.found("dim"))
        {
            const string dimStr = args["dim"];
            IStringStream is(dimStr);
            is >> newDims;
            found = true;
        }
        else
        {
            for (int i = 0; i < catalogueSize; ++i)
            {
                if (fieldName == catalogue[i].name)
                {
                    newDims = catalogue[i].dims;
                    found = true;
                    break;
                }
            }
        }

        if (!found)
        {
            FatalErrorInFunction
                << "Field '" << fieldName << "' is not in the built-in catalogue.\n"
                << "Provide -dim \"[M L T K mol A cd]\" to specify dimensions explicitly."
                << exit(FatalError);
        }

        const fileName fieldPath = timeDir / fieldName;
        patchDimensions(fieldPath, newDims, dryRun);
    }
    else
    {
        // ── Catalogue mode: process all known fields that exist ───────────
        if (args.found("dim"))
        {
            FatalErrorInFunction
                << "-dim requires -field to be specified."
                << exit(FatalError);
        }

        label nFixed = 0;
        label nAlready = 0;
        label nMissing = 0;

        Info<< "Processing built-in catalogue in " << timeDir << ":\n" << endl;

        for (int i = 0; i < catalogueSize; ++i)
        {
            const fileName fieldPath = timeDir / catalogue[i].name;

            if (!isFile(fieldPath))
            {
                Info<< "  skip  " << catalogue[i].name
                    << " (not present in " << timeName << "/)" << nl;
                ++nMissing;
                continue;
            }

            Info<< "  [" << catalogue[i].note << "]" << nl;
            if (patchDimensions(fieldPath, catalogue[i].dims, dryRun))
                ++nFixed;
            else
                ++nAlready;
        }

        Info<< "\nSummary: "
            << nFixed   << " fixed, "
            << nAlready << " already correct, "
            << nMissing << " not present."
            << nl << endl;
    }

    Info<< "========================================\n" << endl;

    return 0;
}

// ************************************************************************* //
