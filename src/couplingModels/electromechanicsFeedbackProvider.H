#ifndef CouplingSignalProvider_H
#define CouplingSignalProvider_H

#include "scalar.H"
#include "label.H"

namespace Foam
{

enum class CouplingSignal // this can be expanded for the necessary eletromechanic feedback. For mechanobiology, see later.
{
    Vm,        // membrane potential mV
    Act,       // activation signal [0..1]
    Cai,       // intracellular Ca  mM
};

class CouplingSignalProvider
{
public:
    virtual ~CouplingSignalProvider() = default;

    // Does this provider support the requested signal?
    virtual bool hasSignal(const CouplingSignal s) const = 0;
    virtual scalar signal(const label i, const CouplingSignal s) const = 0;
};

} // End namespace Foam

#endif

