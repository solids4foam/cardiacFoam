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
    ionicHeterogeneityProbe

Description
    Probe a BuenoOrovio transmural heterogeneity setup with synthetic
    transmural-distance samples. Writes Vm traces, action-potential metrics,
    and adjacent-sample smoothness checks.

Author
    Simao Nieto de Castro. All rights reserved.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "ionicModel.H"
#include "ionicHeterogeneity.H"
#include "Switch.H"

#include <cmath>

using namespace Foam;

namespace
{

dictionary electroModelDict(const IOdictionary& electroDict)
{
    if (electroDict.found("electroModel"))
    {
        word electroModelName;
        electroDict.lookup("electroModel") >> electroModelName;

        const word coeffsName(electroModelName + "Coeffs");

        if (!electroDict.found(coeffsName))
        {
            FatalErrorInFunction
                << "Expected sub-dictionary '" << coeffsName
                << "' in electroProperties for electroModel '"
                << electroModelName << "'."
                << exit(FatalError);
        }

        return electroDict.subDict(coeffsName);
    }

    if (electroDict.found("myocardiumSolver"))
    {
        word myocardiumSolverName;
        electroDict.lookup("myocardiumSolver") >> myocardiumSolverName;

        const word coeffsName(myocardiumSolverName + "Coeffs");

        if (electroDict.found(coeffsName))
        {
            return electroDict.subDict(coeffsName);
        }
    }

    return dictionary(electroDict);
}


struct ProbeMetrics
{
    scalar t;
    scalar rmp;
    scalar peak;
    scalar amplitude;
    scalar maxDvdt;
    scalar apd30;
    scalar apd50;
    scalar apd70;
    scalar apd90;
    bool repolarized;
    bool valid;
    word flags;
};


scalar crossingTime
(
    const scalarField& time,
    const scalarField& vm,
    const label startI,
    const scalar level
)
{
    for (label i = startI + 1; i < vm.size(); ++i)
    {
        if (vm[i - 1] > level && vm[i] <= level)
        {
            const scalar denom = vm[i] - vm[i - 1];
            if (mag(denom) <= VSMALL)
            {
                return time[i];
            }

            const scalar w = (level - vm[i - 1])/denom;
            return time[i - 1] + w*(time[i] - time[i - 1]);
        }
    }

    return -1.0;
}


ProbeMetrics computeMetrics
(
    const scalar t,
    const scalarField& time,
    const scalarField& vm,
    const scalar minAmplitude,
    const scalar minPeak,
    const scalar maxRecoveryRise,
    const bool checkSecondaryRise
)
{
    ProbeMetrics m;
    m.t = t;
    m.rmp = vm[0];
    m.peak = vm[0];
    m.amplitude = 0.0;
    m.maxDvdt = 0.0;
    m.apd30 = -1.0;
    m.apd50 = -1.0;
    m.apd70 = -1.0;
    m.apd90 = -1.0;
    m.repolarized = false;
    m.valid = true;
    m.flags = "ok";

    label peakI = 0;
    label activationI = 0;

    forAll(vm, i)
    {
        if (vm[i] > m.peak)
        {
            m.peak = vm[i];
            peakI = i;
        }

        if (i > 0)
        {
            const scalar dvdt = (vm[i] - vm[i - 1])/(time[i] - time[i - 1]);
            if (dvdt > m.maxDvdt)
            {
                m.maxDvdt = dvdt;
                activationI = i;
            }
        }
    }

    m.amplitude = m.peak - m.rmp;
    const scalar activationTime = time[activationI];

    const scalar level30 = m.rmp + 0.70*m.amplitude;
    const scalar level50 = m.rmp + 0.50*m.amplitude;
    const scalar level70 = m.rmp + 0.30*m.amplitude;
    const scalar level90 = m.rmp + 0.10*m.amplitude;

    const scalar t30 = crossingTime(time, vm, peakI, level30);
    const scalar t50 = crossingTime(time, vm, peakI, level50);
    const scalar t70 = crossingTime(time, vm, peakI, level70);
    const scalar t90 = crossingTime(time, vm, peakI, level90);

    if (t30 >= 0.0) { m.apd30 = t30 - activationTime; }
    if (t50 >= 0.0) { m.apd50 = t50 - activationTime; }
    if (t70 >= 0.0) { m.apd70 = t70 - activationTime; }
    if (t90 >= 0.0) { m.apd90 = t90 - activationTime; }

    m.repolarized = t90 >= 0.0;

    bool secondaryRise = false;
    scalar postPeakMin = vm[peakI];
    for (label i = peakI + 1; i < vm.size(); ++i)
    {
        postPeakMin = min(postPeakMin, vm[i]);
        if (vm[i] - postPeakMin > maxRecoveryRise)
        {
            secondaryRise = true;
            break;
        }
    }

    DynamicList<word> flags;

    if (m.amplitude < minAmplitude)
    {
        flags.append("lowAmplitude");
    }
    if (m.peak < minPeak)
    {
        flags.append("lowPeak");
    }
    if (!m.repolarized)
    {
        flags.append("noAPD90");
    }
    if
    (
        m.apd30 < 0.0
     || m.apd50 < 0.0
     || m.apd70 < 0.0
     || m.apd90 < 0.0
     || !(m.apd30 <= m.apd50 && m.apd50 <= m.apd70 && m.apd70 <= m.apd90)
    )
    {
        flags.append("badAPDOrdering");
    }
    if (checkSecondaryRise && secondaryRise)
    {
        flags.append("secondaryRise");
    }

    if (!flags.empty())
    {
        m.valid = false;
        string s(flags[0]);
        for (label i = 1; i < flags.size(); ++i)
        {
            s += "_";
            s += flags[i];
        }
        m.flags = word(s);
    }

    return m;
}


scalar waveformRMS(const scalarField& a, const scalarField& b)
{
    scalar sum = 0.0;

    forAll(a, i)
    {
        sum += sqr(a[i] - b[i]);
    }

    return std::sqrt(sum/scalar(a.size()));
}


bool bounded(const scalar value, const scalar a, const scalar b, const scalar tol)
{
    return value >= min(a, b) - tol && value <= max(a, b) + tol;
}


void writeMetricsHeader(OFstream& os)
{
    os  << "t,valid,RMP,peak,amplitude,max_dVdt,"
        << "APD30,APD50,APD70,APD90,repolarized,shapeFlags" << nl;
}


void writeMetrics(OFstream& os, const ProbeMetrics& m)
{
    os  << m.t << "," << m.valid << "," << m.rmp << "," << m.peak << ","
        << m.amplitude << "," << m.maxDvdt << "," << m.apd30 << ","
        << m.apd50 << "," << m.apd70 << "," << m.apd90 << ","
        << m.repolarized << "," << m.flags << nl;
}

} // End unnamed namespace


int main(int argc, char *argv[])
{
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "\n========== ionicHeterogeneityProbe ==========\n\n";

    IOdictionary electroDict
    (
        IOobject
        (
            "electroProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary probeDict
    (
        IOobject
        (
            "ionicHeterogeneityProbe",
            runTime.constant(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    dictionary modelDict(electroModelDict(electroDict));

    if (probeDict.found("singleCellStimulus"))
    {
        modelDict.add
        (
            "singleCellStimulus",
            probeDict.subDict("singleCellStimulus"),
            true
        );
    }

    if (!modelDict.found("ionicHeterogeneity"))
    {
        FatalErrorInFunction
            << "No ionicHeterogeneity dictionary found in the selected "
            << "electroProperties model coefficients."
            << exit(FatalError);
    }

    const word ionicModelName(modelDict.lookup("ionicModel"));

    const label nSamples =
        probeDict.lookupOrDefault<label>("nSamples", 101);
    const scalar durationMs =
        probeDict.lookupOrDefault<scalar>("duration", 1000.0);
    const scalar dtMs =
        probeDict.lookupOrDefault<scalar>("dt", 0.1);
    const label writeEvery =
        max(label(1), probeDict.lookupOrDefault<label>("writeEvery", 1));

    if (nSamples < 2 || durationMs <= 0.0 || dtMs <= 0.0)
    {
        FatalErrorInFunction
            << "Invalid probe controls: nSamples=" << nSamples
            << ", duration=" << durationMs
            << ", dt=" << dtMs << "."
            << exit(FatalError);
    }

    const scalar minAmplitude =
        probeDict.lookupOrDefault<scalar>("minAmplitude", 50.0);
    const scalar minPeak =
        probeDict.lookupOrDefault<scalar>("minPeak", 0.0);
    const scalar maxRecoveryRise =
        probeDict.lookupOrDefault<scalar>("maxRecoveryRise", 5.0);
    const Switch checkSecondaryRise =
        probeDict.lookupOrDefault<Switch>("checkSecondaryRise", false);
    const scalar maxAPDJump =
        probeDict.lookupOrDefault<scalar>("maxAPDJump", 20.0);
    const scalar maxWaveformRMS =
        probeDict.lookupOrDefault<scalar>("maxWaveformRMS", 10.0);
    const Switch checkAPDEnvelope =
        probeDict.lookupOrDefault<Switch>("checkAPDEnvelope", true);
    const Switch failOnAPDEnvelope =
        probeDict.lookupOrDefault<Switch>("failOnAPDEnvelope", true);
    const scalar maxAPDBoundTolerance =
        probeDict.lookupOrDefault<scalar>("maxAPDBoundTolerance", 5.0);

    const dictionary& heterogeneityDict =
        modelDict.subDict("ionicHeterogeneity");
    const word heterogeneityMode =
        heterogeneityDict.lookupOrDefault<word>("mode", "transmuralBands");

    if (heterogeneityMode == "cellZoneRegions")
    {
        FatalErrorInFunction
            << "ionicHeterogeneityProbe does not support mode "
            << "cellZoneRegions: there is no continuous distance field to "
            << "sweep. Use transmuralBands or namedRegions."
            << exit(FatalError);
    }

    scalar endoMInterface = 0.3;
    scalar mEpiInterface = 0.7;

    if (heterogeneityMode == "namedRegions")
    {
        // Best-effort defaults for the envelope check: use the boundary
        // between the first and second, and second and third, regions
        // once sorted by range. Users can always override via
        // endoReferenceT/mCellReferenceT/epiReferenceT in the probe dict.
        const List<ionicHeterogeneity::NamedFieldRegion> regions =
            ionicHeterogeneity::parseNamedFieldRegions
            (
                heterogeneityDict.subDict("regions")
            );

        if (regions.size() >= 3)
        {
            endoMInterface = regions[0].rangeMax;
            mEpiInterface = regions[regions.size() - 2].rangeMax;
        }
    }
    else
    {
        endoMInterface =
            heterogeneityDict.lookupOrDefault<scalar>("endoMInterface", 0.3);
        mEpiInterface =
            heterogeneityDict.lookupOrDefault<scalar>("mEpiInterface", 0.7);
    }
    const scalar endoReferenceT =
        probeDict.lookupOrDefault<scalar>("endoReferenceT", 0.0);
    const scalar mCellReferenceT =
        probeDict.lookupOrDefault<scalar>
        (
            "mCellReferenceT",
            0.5*(endoMInterface + mEpiInterface)
        );
    const scalar epiReferenceT =
        probeDict.lookupOrDefault<scalar>("epiReferenceT", 1.0);

    scalarField tSamples(nSamples, 0.0);
    forAll(tSamples, i)
    {
        tSamples[i] = scalar(i)/scalar(nSamples - 1);
    }

    autoPtr<ionicModel> modelPtr =
        ionicModel::New(modelDict, nSamples, dtMs, true);

    ionicModel& model = modelPtr();
    model.configureIonicHeterogeneity
    (
        tSamples,
        modelDict.subDict("ionicHeterogeneity")
    );

    const label nSteps = label(ceil(durationMs/dtMs));
    scalarField times(nSteps + 1, 0.0);
    PtrList<scalarField> vmTraces(nSamples);

    forAll(vmTraces, sampleI)
    {
        vmTraces.set(sampleI, new scalarField(nSteps + 1, 0.0));
        vmTraces[sampleI][0] = model.signal(sampleI, CouplingSignal::VM);
    }

    scalarField Vm(nSamples, 0.0);
    scalarField Im(nSamples, 0.0);

    for (label stepI = 1; stepI <= nSteps; ++stepI)
    {
        const scalar t0 = (stepI - 1)*dtMs/1000.0;
        model.solveODE(t0, dtMs/1000.0, Vm, Im);
        times[stepI] = stepI*dtMs;

        forAll(vmTraces, sampleI)
        {
            vmTraces[sampleI][stepI] =
                model.signal(sampleI, CouplingSignal::VM);
        }
    }

    const fileName outputDir
    (
        runTime.path()/"postProcessing"/"ionicHeterogeneityProbe"
    );
    mkDir(outputDir);

    OFstream traces(outputDir/"Vm_traces.csv");
    traces << "time,t,Vm" << nl;

    forAll(vmTraces, sampleI)
    {
        forAll(times, stepI)
        {
            if (stepI % writeEvery == 0 || stepI == times.size() - 1)
            {
                traces
                    << times[stepI] << ","
                    << tSamples[sampleI] << ","
                    << vmTraces[sampleI][stepI] << nl;
            }
        }
    }

    List<ProbeMetrics> metrics(nSamples);

    OFstream metricsFile(outputDir/"AP_metrics.csv");
    writeMetricsHeader(metricsFile);

    forAll(metrics, sampleI)
    {
        metrics[sampleI] =
            computeMetrics
            (
                tSamples[sampleI],
                times,
                vmTraces[sampleI],
                minAmplitude,
                minPeak,
                maxRecoveryRise,
                checkSecondaryRise
            );
        writeMetrics(metricsFile, metrics[sampleI]);
    }

    OFstream smoothnessFile(outputDir/"smoothness_report.csv");
    smoothnessFile
        << "t0,t1,dAPD30,dAPD50,dAPD70,dAPD90,waveformRMS,"
        << "validTransition" << nl;

    label invalidSamples = 0;
    label invalidTransitions = 0;
    label invalidAPDEnvelopes = 0;

    forAll(metrics, sampleI)
    {
        if (!metrics[sampleI].valid)
        {
            ++invalidSamples;
        }

        if (sampleI == 0)
        {
            continue;
        }

        const ProbeMetrics& a = metrics[sampleI - 1];
        const ProbeMetrics& b = metrics[sampleI];
        const scalar dAPD30 = mag(b.apd30 - a.apd30);
        const scalar dAPD50 = mag(b.apd50 - a.apd50);
        const scalar dAPD70 = mag(b.apd70 - a.apd70);
        const scalar dAPD90 = mag(b.apd90 - a.apd90);
        const scalar rms = waveformRMS(vmTraces[sampleI - 1], vmTraces[sampleI]);

        const bool validTransition =
            a.valid
         && b.valid
         && dAPD30 <= maxAPDJump
         && dAPD50 <= maxAPDJump
         && dAPD70 <= maxAPDJump
         && dAPD90 <= maxAPDJump
         && rms <= maxWaveformRMS;

        if (!validTransition)
        {
            ++invalidTransitions;
        }

        smoothnessFile
            << a.t << "," << b.t << ","
            << dAPD30 << "," << dAPD50 << ","
            << dAPD70 << "," << dAPD90 << ","
            << rms << "," << validTransition << nl;
    }

    if (checkAPDEnvelope)
    {
        const label endoRefI =
            min
            (
                label(nSamples - 1),
                max(label(0), label(round(endoReferenceT*(nSamples - 1))))
            );
        const label mCellRefI =
            min
            (
                label(nSamples - 1),
                max(label(0), label(round(mCellReferenceT*(nSamples - 1))))
            );
        const label epiRefI =
            min
            (
                label(nSamples - 1),
                max(label(0), label(round(epiReferenceT*(nSamples - 1))))
            );

        const ProbeMetrics& endoRef = metrics[endoRefI];
        const ProbeMetrics& mCellRef = metrics[mCellRefI];
        const ProbeMetrics& epiRef = metrics[epiRefI];

        OFstream envelopeFile(outputDir/"APD_envelope_report.csv");
        envelopeFile
            << "t,zone,referenceA,referenceB,"
            << "APD30Bounded,APD50Bounded,APD70Bounded,APD90Bounded,"
            << "validEnvelope" << nl;

        forAll(metrics, sampleI)
        {
            const ProbeMetrics& m = metrics[sampleI];
            const bool endoToM = m.t <= mCellReferenceT;
            const ProbeMetrics& a = endoToM ? endoRef : mCellRef;
            const ProbeMetrics& b = endoToM ? mCellRef : epiRef;

            const bool apd30Bounded =
                bounded(m.apd30, a.apd30, b.apd30, maxAPDBoundTolerance);
            const bool apd50Bounded =
                bounded(m.apd50, a.apd50, b.apd50, maxAPDBoundTolerance);
            const bool apd70Bounded =
                bounded(m.apd70, a.apd70, b.apd70, maxAPDBoundTolerance);
            const bool apd90Bounded =
                bounded(m.apd90, a.apd90, b.apd90, maxAPDBoundTolerance);

            const bool validEnvelope =
                m.valid
             && a.valid
             && b.valid
             && apd30Bounded
             && apd50Bounded
             && apd70Bounded
             && apd90Bounded;

            if (!validEnvelope)
            {
                ++invalidAPDEnvelopes;
            }

            envelopeFile
                << m.t << ","
                << (endoToM ? "endoM" : "mEpi") << ","
                << a.t << "," << b.t << ","
                << apd30Bounded << ","
                << apd50Bounded << ","
                << apd70Bounded << ","
                << apd90Bounded << ","
                << validEnvelope << nl;
        }
    }

    Info<< "Wrote " << outputDir/"Vm_traces.csv" << nl
        << "Wrote " << outputDir/"AP_metrics.csv" << nl
        << "Wrote " << outputDir/"smoothness_report.csv" << nl
        << "Invalid samples: " << invalidSamples << "/" << nSamples << nl
        << "Invalid adjacent transitions: " << invalidTransitions << "/"
        << (nSamples - 1) << nl
        << "Invalid APD envelopes: " << invalidAPDEnvelopes << "/"
        << nSamples << nl
        << "\nCompleted.\n";

    return
        invalidSamples
     || invalidTransitions
     || (failOnAPDEnvelope && invalidAPDEnvelopes)
      ? 1
      : 0;
}

// ************************************************************************* //
