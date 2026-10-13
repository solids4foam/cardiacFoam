# Lead-field ECG file inventory

Paths are relative to the cardiacFoam repository. This inventory is the
feature diff against PR #59, not the combined diff of PR #59 against main.
A means added; M means an existing file was extended. Generated cases,
meshes, fields, logs and local historical evidence are excluded. The inventory
includes the homogeneous-template eikonal extension.

| Change | File | Responsibility |
|---|---|---|
| M | [.gitignore](../.gitignore) | Generated/local artifact rules |
| M | [applications/utilities/README.md](../applications/utilities/README.md) | Documentation |
| A | [applications/utilities/buildLeadFields/Make/files](../applications/utilities/buildLeadFields/Make/files) | Build registration |
| A | [applications/utilities/buildLeadFields/Make/options](../applications/utilities/buildLeadFields/Make/options) | Build registration |
| A | [applications/utilities/buildLeadFields/README.md](../applications/utilities/buildLeadFields/README.md) | Documentation |
| A | [applications/utilities/buildLeadFields/buildLeadFields.C](../applications/utilities/buildLeadFields/buildLeadFields.C) | Static construction or verification application |
| A | [applications/utilities/buildLeadFields/examples/wilsonLeads](../applications/utilities/buildLeadFields/examples/wilsonLeads) | Static construction or verification application |
| A | [applications/utilities/buildLeadFields/utility.manifest.toml](../applications/utilities/buildLeadFields/utility.manifest.toml) | Utility discovery |
| A | [applications/utilities/checkLeadFieldRecovery/Make/files](../applications/utilities/checkLeadFieldRecovery/Make/files) | Build registration |
| A | [applications/utilities/checkLeadFieldRecovery/Make/options](../applications/utilities/checkLeadFieldRecovery/Make/options) | Build registration |
| A | [applications/utilities/checkLeadFieldRecovery/README.md](../applications/utilities/checkLeadFieldRecovery/README.md) | Documentation |
| A | [applications/utilities/checkLeadFieldRecovery/checkLeadFieldRecovery.C](../applications/utilities/checkLeadFieldRecovery/checkLeadFieldRecovery.C) | Static construction or verification application |
| A | [applications/utilities/checkLeadFieldRecovery/utility.manifest.toml](../applications/utilities/checkLeadFieldRecovery/utility.manifest.toml) | Utility discovery |
| A | [applications/utilities/setManufacturedLeadFields/Make/files](../applications/utilities/setManufacturedLeadFields/Make/files) | Build registration |
| A | [applications/utilities/setManufacturedLeadFields/Make/options](../applications/utilities/setManufacturedLeadFields/Make/options) | Build registration |
| A | [applications/utilities/setManufacturedLeadFields/README.md](../applications/utilities/setManufacturedLeadFields/README.md) | Documentation |
| A | [applications/utilities/setManufacturedLeadFields/setManufacturedLeadFields.C](../applications/utilities/setManufacturedLeadFields/setManufacturedLeadFields.C) | Static construction or verification application |
| A | [applications/utilities/setManufacturedLeadFields/utility.manifest.toml](../applications/utilities/setManufacturedLeadFields/utility.manifest.toml) | Utility discovery |
| A | [docs/DUAL_ECG_VERIFICATION.md](DUAL_ECG_VERIFICATION.md) | Documentation |
| A | [docs/EIKONAL_LEAD_FIELD_ECG.md](EIKONAL_LEAD_FIELD_ECG.md) | Documentation |
| A | [docs/LEAD_FIELD_ECG_FILES.md](LEAD_FIELD_ECG_FILES.md) | Documentation |
| A | [docs/LEAD_FIELD_ECG_IMPLEMENTATION.md](LEAD_FIELD_ECG_IMPLEMENTATION.md) | Documentation |
| A | [docs/LEAD_FIELD_ECG_NEXT_STEPS.md](LEAD_FIELD_ECG_NEXT_STEPS.md) | Documentation |
| A | [docs/LEAD_FIELD_ECG_PR_VALIDATION.md](LEAD_FIELD_ECG_PR_VALIDATION.md) | Documentation |
| M | [src/electroModels/Make/files](../src/electroModels/Make/files) | Build registration |
| M | [src/electroModels/Make/options](../src/electroModels/Make/options) | Build registration |
| M | [src/electroModels/core/ARCHITECTURE.md](../src/electroModels/core/ARCHITECTURE.md) | Documentation |
| M | [src/electroModels/core/system/electrophysicsSystemBuilder.C](../src/electroModels/core/system/electrophysicsSystemBuilder.C) | Production interface/implementation |
| M | [src/electroModels/ecgModels/README.md](../src/electroModels/ecgModels/README.md) | Documentation |
| M | [src/electroModels/ecgModels/eikonalECG/README.md](../src/electroModels/ecgModels/eikonalECG/README.md) | Documentation |
| M | [src/electroModels/ecgModels/eikonalECG/eikonalECG.C](../src/electroModels/ecgModels/eikonalECG/eikonalECG.C) | Production interface/implementation |
| M | [src/electroModels/ecgModels/eikonalECG/eikonalECG.H](../src/electroModels/ecgModels/eikonalECG/eikonalECG.H) | Production interface/implementation |
| A | [src/electroModels/ecgModels/eikonalECG/eikonalTemplateSource.C](../src/electroModels/ecgModels/eikonalECG/eikonalTemplateSource.C) | Production interface/implementation |
| A | [src/electroModels/ecgModels/eikonalECG/eikonalTemplateSource.H](../src/electroModels/ecgModels/eikonalECG/eikonalTemplateSource.H) | Production interface/implementation |
| A | [src/electroModels/ecgModels/eikonalLeadFieldECG/README.md](../src/electroModels/ecgModels/eikonalLeadFieldECG/README.md) | Documentation |
| A | [src/electroModels/ecgModels/eikonalLeadFieldECG/eikonalLeadFieldECG.C](../src/electroModels/ecgModels/eikonalLeadFieldECG/eikonalLeadFieldECG.C) | Production interface/implementation |
| A | [src/electroModels/ecgModels/eikonalLeadFieldECG/eikonalLeadFieldECG.H](../src/electroModels/ecgModels/eikonalLeadFieldECG/eikonalLeadFieldECG.H) | Production interface/implementation |
| A | [src/electroModels/ecgModels/leadFieldECG/README.md](../src/electroModels/ecgModels/leadFieldECG/README.md) | Documentation |
| A | [src/electroModels/ecgModels/leadFieldECG/leadFieldData.C](../src/electroModels/ecgModels/leadFieldECG/leadFieldData.C) | Production interface/implementation |
| A | [src/electroModels/ecgModels/leadFieldECG/leadFieldData.H](../src/electroModels/ecgModels/leadFieldECG/leadFieldData.H) | Production interface/implementation |
| A | [src/electroModels/ecgModels/leadFieldECG/leadFieldECG.C](../src/electroModels/ecgModels/leadFieldECG/leadFieldECG.C) | Production interface/implementation |
| A | [src/electroModels/ecgModels/leadFieldECG/leadFieldECG.H](../src/electroModels/ecgModels/leadFieldECG/leadFieldECG.H) | Production interface/implementation |
| M | [src/electroModels/ecgModels/torsoECG/torsoECG.C](../src/electroModels/ecgModels/torsoECG/torsoECG.C) | Production interface/implementation |
| M | [src/electroModels/electroDomains/ecgDomain/ecgDomain.C](../src/electroModels/electroDomains/ecgDomain/ecgDomain.C) | Production interface/implementation |
| M | [src/electroModels/electroDomains/ecgDomain/ecgDomain.H](../src/electroModels/electroDomains/ecgDomain/ecgDomain.H) | Production interface/implementation |
| M | [src/electroModels/electroDomains/ecgDomain/ecgSolver.C](../src/electroModels/electroDomains/ecgDomain/ecgSolver.C) | Production interface/implementation |
| M | [src/electroModels/electroDomains/ecgDomain/ecgSolver.H](../src/electroModels/electroDomains/ecgDomain/ecgSolver.H) | Production interface/implementation |
| M | [src/genericWriter/conductivityFieldIO.C](../src/genericWriter/conductivityFieldIO.C) | Production interface/implementation |
| M | [src/genericWriter/conductivityFieldIO.H](../src/genericWriter/conductivityFieldIO.H) | Production interface/implementation |
| M | [src/genericWriter/ecgModelIO.C](../src/genericWriter/ecgModelIO.C) | Production interface/implementation |
| M | [src/genericWriter/ecgModelIO.H](../src/genericWriter/ecgModelIO.H) | Production interface/implementation |
| M | [src/verificationModels/Make/files](../src/verificationModels/Make/files) | Build registration |
| M | [src/verificationModels/README.md](../src/verificationModels/README.md) | Documentation |
| A | [src/verificationModels/ecgVerification/manufacturedEikonalLeadFieldECGVerifier.C](../src/verificationModels/ecgVerification/manufacturedEikonalLeadFieldECGVerifier.C) | Independent verification interface/implementation |
| A | [src/verificationModels/ecgVerification/manufacturedEikonalLeadFieldECGVerifier.H](../src/verificationModels/ecgVerification/manufacturedEikonalLeadFieldECGVerifier.H) | Independent verification interface/implementation |
| A | [src/verificationModels/ecgVerification/manufacturedLeadFieldECGVerifier.C](../src/verificationModels/ecgVerification/manufacturedLeadFieldECGVerifier.C) | Independent verification |
| A | [src/verificationModels/ecgVerification/manufacturedLeadFieldECGVerifier.H](../src/verificationModels/ecgVerification/manufacturedLeadFieldECGVerifier.H) | Independent verification |
| M | [tutorials/Alltest-regression](../tutorials/Alltest-regression) | Tutorial inputs and verification controls |
| M | [tutorials/README.md](../tutorials/README.md) | Documentation |
| M | [tutorials/configurationChecks/README.md](../tutorials/configurationChecks/README.md) | Documentation |
| A | [tutorials/configurationChecks/base/leadFieldColumn/Allprepare](../tutorials/configurationChecks/base/leadFieldColumn/Allprepare) | Maintained case runner/cleanup |
| A | [tutorials/configurationChecks/base/leadFieldColumn/checkLeadFields.py](../tutorials/configurationChecks/base/leadFieldColumn/checkLeadFields.py) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/constant/physicsProperties](../tutorials/configurationChecks/base/leadFieldColumn/constant/physicsProperties) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/blockMeshDict](../tutorials/configurationChecks/base/leadFieldColumn/system/blockMeshDict) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/controlDict](../tutorials/configurationChecks/base/leadFieldColumn/system/controlDict) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/createBafflesDict](../tutorials/configurationChecks/base/leadFieldColumn/system/createBafflesDict) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/fvSchemes](../tutorials/configurationChecks/base/leadFieldColumn/system/fvSchemes) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/fvSolution](../tutorials/configurationChecks/base/leadFieldColumn/system/fvSolution) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/setTorsoOrganConductivityFieldDict](../tutorials/configurationChecks/base/leadFieldColumn/system/setTorsoOrganConductivityFieldDict) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/topoSetBaffleDict](../tutorials/configurationChecks/base/leadFieldColumn/system/topoSetBaffleDict) | Case input/check registration |
| A | [tutorials/configurationChecks/base/leadFieldColumn/system/topoSetDict](../tutorials/configurationChecks/base/leadFieldColumn/system/topoSetDict) | Case input/check registration |
| M | [tutorials/configurationChecks/checks](../tutorials/configurationChecks/checks) | Case input/check registration |
| M | [tutorials/configurationChecks/comparisons](../tutorials/configurationChecks/comparisons) | Case input/check registration |
| A | [tutorials/configurationChecks/post/leadFieldGuards](../tutorials/configurationChecks/post/leadFieldGuards) | Maintained configuration control |
| A | [tutorials/configurationChecks/post/leadFieldHeartOnly](../tutorials/configurationChecks/post/leadFieldHeartOnly) | Maintained configuration control |
| A | [tutorials/configurationChecks/post/leadFieldParallelBuild](../tutorials/configurationChecks/post/leadFieldParallelBuild) | Maintained configuration control |
| A | [tutorials/configurationChecks/post/leadFieldRecovery](../tutorials/configurationChecks/post/leadFieldRecovery) | Maintained configuration control |
| A | [tutorials/configurationChecks/post/leadFieldRestartRounding](../tutorials/configurationChecks/post/leadFieldRestartRounding) | Maintained configuration control |
| A | [tutorials/configurationChecks/regression/manifests/bidomainLeadFieldBaseline](../tutorials/configurationChecks/regression/manifests/bidomainLeadFieldBaseline) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/bidomainLeadFieldECG](../tutorials/configurationChecks/regression/manifests/bidomainLeadFieldECG) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/bidomainLeadFieldECGParallel](../tutorials/configurationChecks/regression/manifests/bidomainLeadFieldECGParallel) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/bidomainLeadFieldRecovery](../tutorials/configurationChecks/regression/manifests/bidomainLeadFieldRecovery) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/bidomainLeadFieldRecoveryParallel](../tutorials/configurationChecks/regression/manifests/bidomainLeadFieldRecoveryParallel) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldBaseline](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldBaseline) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECG](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECG) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECGField](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECGField) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECGLegacy](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECGLegacy) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECGParallel](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldECGParallel) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldGuards](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldGuards) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldHeartOnly](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldHeartOnly) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldParallelBuild](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldParallelBuild) | Expected configuration output layout |
| A | [tutorials/configurationChecks/regression/manifests/monodomainLeadFieldRestartRounding](../tutorials/configurationChecks/regression/manifests/monodomainLeadFieldRestartRounding) | Expected configuration output layout |
| A | [tutorials/configurationChecks/variants/bidomainLeadFieldECG/constant/electroProperties](../tutorials/configurationChecks/variants/bidomainLeadFieldECG/constant/electroProperties) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/bidomainLeadFieldECG/constant/leadFieldProperties](../tutorials/configurationChecks/variants/bidomainLeadFieldECG/constant/leadFieldProperties) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECG/constant/electroProperties](../tutorials/configurationChecks/variants/monodomainLeadFieldECG/constant/electroProperties) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECG/constant/leadFieldProperties](../tutorials/configurationChecks/variants/monodomainLeadFieldECG/constant/leadFieldProperties) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGField/0/ConductivityExtracellular](../tutorials/configurationChecks/variants/monodomainLeadFieldECGField/0/ConductivityExtracellular) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGField/0/ConductivityIntracellular](../tutorials/configurationChecks/variants/monodomainLeadFieldECGField/0/ConductivityIntracellular) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGField/constant/electroProperties](../tutorials/configurationChecks/variants/monodomainLeadFieldECGField/constant/electroProperties) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGField/constant/leadFieldProperties](../tutorials/configurationChecks/variants/monodomainLeadFieldECGField/constant/leadFieldProperties) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/0/conductivityExtracellular](../tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/0/conductivityExtracellular) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/0/conductivityIntracellular](../tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/0/conductivityIntracellular) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/constant/electroProperties](../tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/constant/electroProperties) | Configuration variant input |
| A | [tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/constant/leadFieldProperties](../tutorials/configurationChecks/variants/monodomainLeadFieldECGLegacy/constant/leadFieldProperties) | Configuration variant input |
| M | [tutorials/manufacturedSolutions/README.md](../tutorials/manufacturedSolutions/README.md) | Documentation |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/Allclean](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/Allclean) | Maintained case runner/cleanup |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/Allrun](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/Allrun) | Maintained case runner/cleanup |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/README.md](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/README.md) | Documentation |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/constant/electroProperties](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/constant/electroProperties) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/constant/leadFieldProperties](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/constant/leadFieldProperties) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/constant/physicsProperties](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/constant/physicsProperties) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/regression/leadFieldECG.reference](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/regression/leadFieldECG.reference) | Maintained numerical regression |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/regression/regressionTest.sh](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/regression/regressionTest.sh) | Maintained numerical regression |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/regression/source.reference](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/regression/source.reference) | Maintained numerical regression |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/setup/pseudoECGObserver](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/setup/pseudoECGObserver) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/blockMeshDict.3D](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/blockMeshDict.3D) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/controlDict](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/controlDict) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/decomposeParDict](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/decomposeParDict) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/fvSchemes](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/fvSchemes) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/fvSolution](../tutorials/manufacturedSolutions/bidomainLeadFieldECG/system/fvSolution) | Case input/check registration |
| M | [tutorials/manufacturedSolutions/eikonalECG/Allclean](../tutorials/manufacturedSolutions/eikonalECG/Allclean) | Maintained tutorial runner |
| M | [tutorials/manufacturedSolutions/eikonalECG/Allrun](../tutorials/manufacturedSolutions/eikonalECG/Allrun) | Maintained tutorial runner |
| M | [tutorials/manufacturedSolutions/eikonalECG/README.md](../tutorials/manufacturedSolutions/eikonalECG/README.md) | Documentation |
| M | [tutorials/manufacturedSolutions/eikonalECG/constant/electroProperties](../tutorials/manufacturedSolutions/eikonalECG/constant/electroProperties) | Tutorial inputs and verification controls |
| A | [tutorials/manufacturedSolutions/eikonalECG/constant/leadFieldProperties](../tutorials/manufacturedSolutions/eikonalECG/constant/leadFieldProperties) | Tutorial inputs and verification controls |
| A | [tutorials/manufacturedSolutions/eikonalECG/setup/leadFieldObserver](../tutorials/manufacturedSolutions/eikonalECG/setup/leadFieldObserver) | Tutorial inputs and verification controls |
| A | [tutorials/manufacturedSolutions/eikonalLeadFieldECG/Allclean](../tutorials/manufacturedSolutions/eikonalLeadFieldECG/Allclean) | Maintained tutorial runner |
| A | [tutorials/manufacturedSolutions/eikonalLeadFieldECG/Allrun](../tutorials/manufacturedSolutions/eikonalLeadFieldECG/Allrun) | Maintained tutorial runner |
| A | [tutorials/manufacturedSolutions/eikonalLeadFieldECG/README.md](../tutorials/manufacturedSolutions/eikonalLeadFieldECG/README.md) | Documentation |
| A | [tutorials/manufacturedSolutions/eikonalLeadFieldECG/regression/leadFieldECG.reference](../tutorials/manufacturedSolutions/eikonalLeadFieldECG/regression/leadFieldECG.reference) | Numerical regression bounds |
| A | [tutorials/manufacturedSolutions/eikonalLeadFieldECG/regression/regressionTest.sh](../tutorials/manufacturedSolutions/eikonalLeadFieldECG/regression/regressionTest.sh) | Tutorial inputs and verification controls |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/Allclean](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/Allclean) | Maintained case runner/cleanup |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/Allrun](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/Allrun) | Maintained case runner/cleanup |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/README.md](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/README.md) | Documentation |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/constant/electroProperties](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/constant/electroProperties) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/constant/leadFieldProperties](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/constant/leadFieldProperties) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/constant/physicsProperties](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/constant/physicsProperties) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/regression/leadFieldECG.reference](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/regression/leadFieldECG.reference) | Maintained numerical regression |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/regression/regressionTest.sh](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/regression/regressionTest.sh) | Maintained numerical regression |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/regression/source.reference](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/regression/source.reference) | Maintained numerical regression |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/setup/pseudoECGObserver](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/setup/pseudoECGObserver) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/blockMeshDict.3D](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/blockMeshDict.3D) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/controlDict](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/controlDict) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/decomposeParDict](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/decomposeParDict) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/fvSchemes](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/fvSchemes) | Case input/check registration |
| A | [tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/fvSolution](../tutorials/manufacturedSolutions/monodomainLeadFieldECG/system/fvSolution) | Case input/check registration |
| A | [tutorials/regression/checkDualECG.py](../tutorials/regression/checkDualECG.py) | Maintained numerical regression |
| A | [tutorials/regression/checkEikonalDualECG.py](../tutorials/regression/checkEikonalDualECG.py) | Tutorial inputs and verification controls |
