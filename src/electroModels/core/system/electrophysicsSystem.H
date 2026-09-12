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

Class
    electrophysicsSystem

Description
    Orchestrates the execution of a multi-domain electrophysiology system
    using a domain separation architecture.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#ifndef electrophysicsSystem_H
#define electrophysicsSystem_H

#include "autoPtr.H"
#include "PtrList.H"
#include "electroDomainInterface.H"
#include "electroStateDomain.H"
#include "electroDomainCoupler.H"
#include "myocardiumDomainInterface.H"
#include "electrophysicsAdvanceScheme.H"

namespace Foam
{

class electrophysicsSystem
{
    autoPtr<myocardiumDomainInterface> myocardium_;
    autoPtr<electroStateDomain> potentialDomain_;
    autoPtr<electrophysicsAdvanceScheme> advanceScheme_;
    PtrList<electroDomainInterface> conductionDomains_;
    PtrList<electroDomainCoupler> conductionCouplingModels_;
    PtrList<electroDomainInterface> ecgDomains_;
    PtrList<electroDomainCoupler> ecgCouplingModels_;

public:

    electrophysicsSystem() = default;

    void setMyocardium(myocardiumDomainInterface* myocardium)
    {
        myocardium_.reset(myocardium);
    }

    void setPotentialDomain(electroStateDomain* domain)
    {
        potentialDomain_.reset(domain);
    }

    void clearPotentialDomain()
    {
        if (potentialDomain_.valid() && myocardium_.valid())
        {
            myocardium_().unbindExternalPhiE();
        }
        potentialDomain_.clear();
    }

    void setAdvanceScheme(electrophysicsAdvanceScheme* advanceScheme)
    {
        advanceScheme_.reset(advanceScheme);
    }

    void appendConductionDomain(electroDomainInterface* domain)
    {
        if (!domain)
        {
            return;
        }

        const label i = conductionDomains_.size();
        conductionDomains_.setSize(i + 1);
        conductionDomains_.set(i, domain);
    }

    void clearConductionDomains()
    {
        conductionDomains_.clear();
    }

    void appendConductionCoupling(electroDomainCoupler* couplingModel)
    {
        if (!couplingModel)
        {
            return;
        }

        const label i = conductionCouplingModels_.size();
        conductionCouplingModels_.setSize(i + 1);
        conductionCouplingModels_.set(i, couplingModel);
    }

    void clearConductionCouplings()
    {
        conductionCouplingModels_.clear();
    }

    void appendECGDomain(electroDomainInterface* domain)
    {
        if (!domain)
        {
            return;
        }

        const label i = ecgDomains_.size();
        ecgDomains_.setSize(i + 1);
        ecgDomains_.set(i, domain);
    }

    void clearECGDomains()
    {
        ecgDomains_.clear();
    }

    void appendECGCoupling(electroDomainCoupler* couplingModel)
    {
        if (!couplingModel)
        {
            return;
        }

        const label i = ecgCouplingModels_.size();
        ecgCouplingModels_.setSize(i + 1);
        ecgCouplingModels_.set(i, couplingModel);
    }

    void clearECGCouplings()
    {
        ecgCouplingModels_.clear();
    }

    bool hasConductionDomains() const
    {
        return !conductionDomains_.empty();
    }

    bool hasConductionCouplings() const
    {
        return !conductionCouplingModels_.empty();
    }

    bool hasECGDomains() const
    {
        return !ecgDomains_.empty();
    }

    bool hasECGCouplings() const
    {
        return !ecgCouplingModels_.empty();
    }

    myocardiumDomainInterface& myocardium()
    {
        return myocardium_();
    }

    const myocardiumDomainInterface& myocardium() const
    {
        return myocardium_();
    }

    bool hasMyocardium() const
    {
        return myocardium_.valid();
    }

    bool hasPotentialDomain() const
    {
        return potentialDomain_.valid();
    }

    const electroStateProvider& potentialStateProvider() const
    {
        return potentialDomain_();
    }

    const PtrList<electroDomainInterface>& conductionDomains() const
    {
        return conductionDomains_;
    }

    const PtrList<electroDomainCoupler>& conductionCouplings() const
    {
        return conductionCouplingModels_;
    }

    const PtrList<electroDomainInterface>& ecgDomains() const
    {
        return ecgDomains_;
    }

    const PtrList<electroDomainCoupler>& ecgCouplings() const
    {
        return ecgCouplingModels_;
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

    void prepareConductionCouplings(scalar t0, scalar dt);

    void advanceConductionDomains(scalar t0, scalar dt);

    void prepareMyocardiumCouplings(scalar t0, scalar dt);

    void prepareECGCouplings(scalar t0, scalar dt);

    void preparePotentialDomain(scalar t0, scalar dt);

    void advancePotentialDomain(scalar t0, scalar dt);

    void advanceECGDomains(scalar t0, scalar dt);

    void writeConductionDomains();

    void writeConductionCouplings();

    void writePotentialDomain();

    void endConductionDomains();

    void writeECGDomains();

    void endConductionCouplings();

    void endECGDomains();

    void endECGCouplings();
};

} // End namespace Foam

#endif

// ************************************************************************* //
