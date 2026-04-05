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

#ifndef electrophysicsSystem_H
#define electrophysicsSystem_H

#include "autoPtr.H"
#include "PtrList.H"
#include "electroDomainInterface.H"
#include "conductionSystemDomain.H"
#include "ecgDomain.H"
#include "electroDomainCoupler.H"
#include "myocardiumDomain.H"
#include "schemes/electrophysicsAdvanceScheme.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class electrophysicsSystem
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~

Description
    Orchestrates the execution of a multi-domain electrophysiology system
    using a domain separation architecture.

    Manages:
    - Primary domain:        myocardium (3D reaction-diffusion)
    - Pre-primary domains:   conduction networks (e.g., Purkinje 1D)
    - Post-primary domains:  secondary computations (e.g., ECG)
    - Coupling models:       bidirectional/unidirectional interactions
    - Advance scheme:        execution strategy (e.g., staggered coupling)

Architecture
    See ARCHITECTURE.md for detailed explanation of:
    - Domain separation philosophy
    - electroStateProvider interface and one-way coupling pattern
    - Dependency inversion between solvers and domains
    - Why this design improves extensibility and testability

\*---------------------------------------------------------------------------*/

class electrophysicsSystem
{
    autoPtr<MyocardiumDomain> myocardium_;
    autoPtr<electrophysicsAdvanceScheme> advanceScheme_;
    PtrList<electroDomainInterface> upstreamDomains_;
    PtrList<ElectroDomainCoupler> upstreamCouplingModels_;
    PtrList<electroDomainInterface> downstreamDomains_;
    PtrList<ElectroDomainCoupler> downstreamCouplingModels_;

public:

    electrophysicsSystem() = default;

    void setMyocardium(MyocardiumDomain* myocardium)
    {
        myocardium_.reset(myocardium);
    }

    void setAdvanceScheme(electrophysicsAdvanceScheme* advanceScheme)
    {
        advanceScheme_.reset(advanceScheme);
    }

    void appendUpstreamDomain(electroDomainInterface* domain)
    {
        if (!domain)
        {
            return;
        }

        const label i = upstreamDomains_.size();
        upstreamDomains_.setSize(i + 1);
        upstreamDomains_.set(i, domain);
    }

    void setUpstreamDomain(electroDomainInterface* domain)
    {
        clearUpstreamDomain();
        appendUpstreamDomain(domain);
    }

    void clearUpstreamDomain()
    {
        upstreamDomains_.clear();
    }

    void setConductionSystemDomain
    (
        ConductionSystemDomain* conductionSystemDomain
    )
    {
        setUpstreamDomain(conductionSystemDomain);
    }

    void clearConductionSystemDomain()
    {
        clearUpstreamDomain();
    }

    void setUpstreamCouplingModel(ElectroDomainCoupler* couplingModel)
    {
        clearUpstreamCouplingModel();
        appendUpstreamCouplingModel(couplingModel);
    }

    void appendUpstreamCouplingModel(ElectroDomainCoupler* couplingModel)
    {
        if (!couplingModel)
        {
            return;
        }

        const label i = upstreamCouplingModels_.size();
        upstreamCouplingModels_.setSize(i + 1);
        upstreamCouplingModels_.set(i, couplingModel);
    }

    void clearUpstreamCouplingModel()
    {
        upstreamCouplingModels_.clear();
    }

    void appendDownstreamDomain(electroDomainInterface* domain)
    {
        if (!domain)
        {
            return;
        }

        const label i = downstreamDomains_.size();
        downstreamDomains_.setSize(i + 1);
        downstreamDomains_.set(i, domain);
    }

    void setDownstreamDomain(electroDomainInterface* domain)
    {
        clearDownstreamDomain();
        appendDownstreamDomain(domain);
    }

    void clearDownstreamDomain()
    {
        downstreamDomains_.clear();
    }

    void setECGDomain(ECGDomain* ecgDomain)
    {
        setDownstreamDomain(ecgDomain);
    }

    void clearECGDomain()
    {
        clearDownstreamDomain();
    }

    void setDownstreamCouplingModel
    (
        ElectroDomainCoupler* couplingModel
    )
    {
        clearDownstreamCouplingModel();
        appendDownstreamCouplingModel(couplingModel);
    }

    void appendDownstreamCouplingModel
    (
        ElectroDomainCoupler* couplingModel
    )
    {
        if (!couplingModel)
        {
            return;
        }

        const label i = downstreamCouplingModels_.size();
        downstreamCouplingModels_.setSize(i + 1);
        downstreamCouplingModels_.set(i, couplingModel);
    }

    void clearDownstreamCouplingModel()
    {
        downstreamCouplingModels_.clear();
    }

    bool hasUpstreamDomain() const
    {
        return !upstreamDomains_.empty();
    }

    bool hasConductionSystemDomain() const
    {
        return hasUpstreamDomain();
    }

    bool hasUpstreamCouplingModel() const
    {
        return !upstreamCouplingModels_.empty();
    }

    bool hasDownstreamDomain() const
    {
        return !downstreamDomains_.empty();
    }

    bool hasECGDomain() const
    {
        return hasDownstreamDomain();
    }

    bool hasDownstreamCouplingModel() const
    {
        return !downstreamCouplingModels_.empty();
    }

    MyocardiumDomain& myocardium()
    {
        return myocardium_();
    }

    const MyocardiumDomain& myocardium() const
    {
        return myocardium_();
    }

    bool hasMyocardium() const
    {
        return myocardium_.valid();
    }

    const PtrList<electroDomainInterface>& upstreamDomains() const
    {
        return upstreamDomains_;
    }

    const PtrList<ElectroDomainCoupler>& upstreamCouplingModels() const
    {
        return upstreamCouplingModels_;
    }

    const PtrList<electroDomainInterface>& downstreamDomains() const
    {
        return downstreamDomains_;
    }

    const PtrList<ElectroDomainCoupler>& downstreamCouplingModels() const
    {
        return downstreamCouplingModels_;
    }

    bool advance
    (
        scalar t0,
        scalar dt,
        pimpleControl* pimplePtr,
        electrophysicsAdvanceTimings& timings
    )
    {
        return advanceScheme_->advance
        (
            t0,
            dt,
            *this,
            pimplePtr,
            timings
        );
    }

    void prepareUpstreamCouplings(scalar t0, scalar dt);

    void advanceUpstreamDomains(scalar t0, scalar dt);

    void preparePrimaryCouplings(scalar t0, scalar dt);

    void prepareDownstreamCouplings(scalar t0, scalar dt);

    void advanceDownstreamDomains(scalar t0, scalar dt);

    void writeUpstreamDomain();

    void writeConductionSystemDomain()
    {
        writeUpstreamDomain();
    }

    void writeUpstreamCouplingModel();

    void endConductionSystemDomain()
    {
        endUpstreamDomain();
    }

    void endUpstreamDomain();

    void writeDownstreamDomain();

    void writeECGDomain()
    {
        writeDownstreamDomain();
    }

    void endUpstreamCouplingModel();

    void endECGDomain()
    {
        endDownstreamDomain();
    }

    void endDownstreamDomain();

    void endDownstreamCouplingModel();
};

} // End namespace Foam

#endif

// ************************************************************************* //
