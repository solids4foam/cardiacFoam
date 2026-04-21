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
#include "electroDomainCoupler.H"
#include "myocardiumDomainInterface.H"
#include "electrophysicsAdvanceScheme.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class electrophysicsSystem
                    ~~~~~~~~~~~~~~~~~~~~~~~~~~

Description
    Orchestrates the execution of a multi-domain electrophysiology system
    using a domain separation architecture.

    Terminology used by this container:
    - primary domain        = myocardium
    - conduction domains    = pre-myocardium graph / Purkinje domains
    - ECG domains           = post-myocardium ECG domains advanced after myocardium

Architecture
    See ARCHITECTURE.md for detailed explanation of:
    - Domain separation philosophy
    - electroStateProvider interface and one-way coupling pattern
    - Dependency inversion between solvers and domains
    - Why this design improves extensibility and testability

\*---------------------------------------------------------------------------*/

class electrophysicsSystem
{
    autoPtr<myocardiumDomainInterface> myocardium_;
    autoPtr<electrophysicsAdvanceScheme> advanceScheme_;
    PtrList<electroDomainInterface> conductionDomains_;
    PtrList<ElectroDomainCoupler> conductionCouplingModels_;
    PtrList<electroDomainInterface> ecgDomains_;
    PtrList<ElectroDomainCoupler> ecgCouplingModels_;

public:

    electrophysicsSystem() = default;

    void setMyocardium(myocardiumDomainInterface* myocardium)
    {
        myocardium_.reset(myocardium);
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

    void appendConductionCoupling(ElectroDomainCoupler* couplingModel)
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

    void appendECGCoupling(ElectroDomainCoupler* couplingModel)
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

    const PtrList<electroDomainInterface>& conductionDomains() const
    {
        return conductionDomains_;
    }

    const PtrList<ElectroDomainCoupler>& conductionCouplings() const
    {
        return conductionCouplingModels_;
    }

    const PtrList<electroDomainInterface>& ecgDomains() const
    {
        return ecgDomains_;
    }

    const PtrList<ElectroDomainCoupler>& ecgCouplings() const
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

    void advanceECGDomains(scalar t0, scalar dt);

    void writeConductionDomains();

    void writeConductionCouplings();

    void endConductionDomains();

    void writeECGDomains();

    void endConductionCouplings();

    void endECGDomains();

    void endECGCouplings();
};

} // End namespace Foam

#endif

// ************************************************************************* //
