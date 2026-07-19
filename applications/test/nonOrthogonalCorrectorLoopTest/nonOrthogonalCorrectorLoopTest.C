#include "nonOrthogonalCorrectorLoop.H"
#include "IOstreams.H"
#include <cassert>

using namespace Foam;

int main()
{
    // count overload: N correctors => N+1 assemble-and-solve passes
    for (const label n : {label(0), label(1), label(2), label(5)})
    {
        label calls = 0;
        correctNonOrthogonalLoop(n, [&]() { ++calls; });
        assert(calls == n + 1);
    }

    // negative count clamps to a single pass
    {
        label calls = 0;
        correctNonOrthogonalLoop(label(-3), [&]() { ++calls; });
        assert(calls == 1);
    }

    Info<< "nonOrthogonalCorrectorLoop count semantics OK" << endl;
    return 0;
}
