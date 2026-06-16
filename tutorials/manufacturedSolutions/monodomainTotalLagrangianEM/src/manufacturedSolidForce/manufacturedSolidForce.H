#ifndef manufacturedSolidForce_H
#define manufacturedSolidForce_H

#include "fvOption.H"

namespace Foam
{
namespace fv
{

class manufacturedSolidForce
:
    public option
{
    vector amplitude_;
    scalar Tmax_;
    scalar V0_;
    scalar gamma_;
        //- Passive mechanical parameters
        scalar E_;
        scalar nu_;

        //- Density
        scalar rho0_;

public:

    TypeName("manufacturedSolidForce");

    manufacturedSolidForce
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~manufacturedSolidForce();

    virtual bool isActive() { return true; }

    virtual void addSup
    (
        fvMatrix<vector>& eqn,
        const label fieldi
    );

    virtual bool read(const dictionary& dict);
};

}
}

#endif
