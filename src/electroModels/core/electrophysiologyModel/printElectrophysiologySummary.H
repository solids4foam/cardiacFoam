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

Description
    Code snippet to print a formatted summary of the electrophysiology solver
    configuration, including myocardium, conduction networks, and ECG domains.

Author
    Simao Nieto de Castro, UCD.
\*---------------------------------------------------------------------------*/

{
    const word myoSolver(this->lookupOrDefault<word>("myocardiumSolver", "none"));
    if (myoSolver != "none" && this->found(myoSolver + "Coeffs"))
    {
        const dictionary& coeffs = this->subDict(myoSolver + "Coeffs");

        const word advanceScheme(this->lookupOrDefault<word>("electrophysicsAdvanceScheme", "staggeredElectrophysicsAdvanceScheme"));

        const label meshDim = this->mesh().nGeometricD();
        string dimStr = Foam::name(meshDim) + "D";
        if (coeffs.found("conductionNetworkDomains"))
        {
            dimStr += "-1D";
        }

        string solverHeader = dimStr + " solver";
        string solverUnderline(solverHeader.size(), '-');

        Info<< nl
            << "=========================================================" << nl
            << "         _   _                                           " << nl
            << "        / \\_/ \\                                          " << nl
            << "        \\     /          Cardiac Electrophysiology       " << nl
            << "         \\   /                                           " << nl
            << "          \\_/                                            " << nl
            << "=========================================================" << nl
            << "  Advance Scheme      : " << advanceScheme << nl << nl
            << "  " << solverHeader << nl
            << "  " << solverUnderline << nl;

        // 1. Myocardium
        const word myoIonic = coeffs.lookupOrDefault<word>("ionicModel", "none");
        Info<< "  myocardium          : " << myoSolver << nl
            << "    Ionic Model       : " << myoIonic << nl;

        if (coeffs.found("ionicHeterogeneity"))
        {
            const dictionary& het = coeffs.subDict("ionicHeterogeneity");
            const word mode = het.lookupOrDefault<word>("mode", "none");
            if (mode == "transmuralBands")
            {
                Info<< "    Heterogeneity     : endo |"
                    << het.lookupOrDefault<scalar>("endoMInterface", 0.0) << "| mcells |"
                    << het.lookupOrDefault<scalar>("mEpiInterface", 0.0) << "| epi" << nl;
            }
            else
            {
                Info<< "    Heterogeneity     : " << mode << nl;
            }
        }
        else
        {
            const word tissue = coeffs.lookupOrDefault<word>("tissue", "none");
            if (tissue != "none")
            {
                Info<< "    Tissue Region     : " << tissue << nl;
            }
        }

        // 2. Conduction / Purkinje
        if (coeffs.found("conductionNetworkDomains"))
        {
            const dictionary& condDicts = coeffs.subDict("conductionNetworkDomains");
            forAllConstIters(condDicts, iter)
            {
                if (iter().isDict())
                {
                    const word condName = iter().keyword();
                    const dictionary& condDict = iter().dict();
                    const word condSolverWrapper = condDict.lookupOrDefault<word>("conductionSystemDomain", "none");

                    if (condSolverWrapper != "none" && condDict.found(condSolverWrapper + "Coeffs"))
                    {
                        const dictionary& condCoeffs = condDict.subDict(condSolverWrapper + "Coeffs");
                        const word innerSolver = condCoeffs.lookupOrDefault<word>("conductionSystemSolver", condSolverWrapper);
                        const word cIonic = condCoeffs.lookupOrDefault<word>("ionicModel", "none");

                        string padCond("  " + condName);
                        while(padCond.size() < 22) padCond += " ";

                        Info<< nl << padCond.c_str() << ": " << innerSolver << nl;
                        if (cIonic != "none")
                        {
                            Info<< "    Ionic Model       : " << cIonic << nl;
                        }

                        if (condCoeffs.found("ionicHeterogeneity"))
                        {
                            const dictionary& het = condCoeffs.subDict("ionicHeterogeneity");
                            const word mode = het.lookupOrDefault<word>("mode", "none");
                            if (mode == "transmuralBands")
                            {
                                Info<< "    Heterogeneity     : endo |"
                                    << het.lookupOrDefault<scalar>("endoMInterface", 0.0) << "| mcells |"
                                    << het.lookupOrDefault<scalar>("mEpiInterface", 0.0) << "| epi" << nl;
                            }
                            else
                            {
                                Info<< "    Heterogeneity     : " << mode << nl;
                            }
                        }
                        else
                        {
                            const word cTissue = condCoeffs.lookupOrDefault<word>("tissue", "none");
                            if (cTissue != "none")
                            {
                                Info<< "    Tissue Region     : " << cTissue << nl;
                            }
                        }
                    }
                    else
                    {
                        string padCond("  " + condName);
                        while(padCond.size() < 22) padCond += " ";
                        Info<< nl << padCond.c_str() << ": " << condSolverWrapper << nl;
                    }
                }
            }
        }

        // 3. ECG Domains
        if (coeffs.found("ecgDomains"))
        {
            const dictionary& ecgDicts = coeffs.subDict("ecgDomains");
            forAllConstIters(ecgDicts, iter)
            {
                if (iter().isDict())
                {
                    const word ecgName = iter().keyword();
                    const dictionary& ecgDict = iter().dict();
                    const word ecgSolver = ecgDict.lookupOrDefault<word>("ecgSolver", "none");

                    string padEcg("  " + ecgName);
                    while(padEcg.size() < 22) padEcg += " ";

                    Info<< nl << padEcg.c_str() << ": " << ecgSolver << nl;

                    if (ecgDict.found("electrodePositions"))
                    {
                        const dictionary& leadsDict = ecgDict.subDict("electrodePositions");
                        wordList leads = leadsDict.toc();
                        Info<< "    Leads             : " << leads.size() << " leads" << nl;
                    }
                }
            }
        }

        Info<< "=========================================================" << nl << endl;
    }
}
