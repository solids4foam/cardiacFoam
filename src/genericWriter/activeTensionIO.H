#ifndef activeTensionIO_H
#define activeTensionIO_H

#include "OFstream.H"
#include "scalarField.H"
#include "wordList.H"
#include "volFields.H"
#include "error.H"

namespace Foam
{

class activeTensionIO
{
public:
    // Return user list when provided, otherwise fall back to defaults
    static wordList exportedFieldNames
    (
        const wordList& userList,
        const wordList& defaultNames
    );

    static void writeHeader
    (
        OFstream& os,
        const wordList& names
    );

    static void write
    (
        const scalar t,
        OFstream& os,
        const scalarField& values
    );

    static void exportStateFields
    (
        const scalarField& Ta,
        const wordList& exportedNames,
        PtrList<volScalarField>& outFields
    );
};

} // End namespace Foam

#endif
