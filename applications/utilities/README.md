# cardiacFoam utilities

`applications/Allwmake` builds every utility directory below in both full and
lightweight modes. The executable names come from each `Make/files`; each linked
README is maintained beside its utility implementation.

| Purpose | Executable/runtime name | Build mode/backend | Source boundary | Known equivalence limitation | Documentation |
|---|---|---|---|---|---|
| Convert a VTK line graph to an OpenFOAM edge mesh | `1DgraphToFoam` | full + lightweight; serial | maintained C++ | Conversion preserves supported data only; validate imported topology | [README](1DgraphToFoam/README.md) |
| Detect or explicitly rescale mesh coordinates | `checkMeshGeometry` | full + lightweight; serial | maintained C++ | Detection is read-only unless writing is explicitly requested | [README](checkMeshGeometry/README.md) |
| Inspect ionic heterogeneity weights and assignments | `ionicHeterogeneityProbe` | full + lightweight; serial | maintained C++ | Diagnostic output does not replace a solver regression | [README](ionicHeterogeneityProbe/README.md) |
| List variables exposed by the selected ionic model | `listCellModelsVariables` | full + lightweight; serial | maintained C++ and generated model metadata consumer | Reports metadata; it does not prove runtime export availability | [README](listCellModelsVariables/README.md) |
| Convert an unstructured VTK mesh to OpenFOAM | `newVtkUnstructuredToFoam` | full + lightweight; serial | maintained C++ reader/converter | Conversion supports the documented VTK subset only | [README](newVtkUnstructuredToFoam/README.md) |
| Recompute pseudo-ECG output from stored fields | `recomputePseudoECG` | full + lightweight; serial | maintained C++ | Equivalence requires identical stored fields, electrodes, and settings | [README](recomputePseudoECG/README.md) |
| Advance one configured Purkinje graph domain | `runPurkinjeGraph` | full + lightweight; serial | maintained C++ | Standalone graph evolution omits myocardium coupling currents | [README](runPurkinjeGraph/README.md) |
| Create myocardial fibre/sheet fields | `setFibreField` | full + lightweight; utility-defined parallel behavior | maintained C++ | Field equivalence depends on mesh, patches, and method settings | [README](setFibreField/README.md) |
| Set field dimensions in an existing field file | `setFieldDimensions` | full + lightweight; serial | maintained C++ | Changes dimensions metadata, not field values | [README](setFieldDimensions/README.md) |
| Assign torso-organ conductivity fields | `setTorsoOrganConductivityField` | full + lightweight; serial | maintained C++ | Result depends on input labels and conductivity dictionary | [README](setTorsoOrganConductivityField/README.md) |
| Diagnose manufactured bath-bidomain interface fluxes | `bathBidomainInterfaceMetrics` | full + lightweight; serial/reconstructed | maintained C++ | One-sided fluxes depend on the configured gradient scheme | [README](bathBidomainInterfaceMetrics/README.md) |
| Sweep configured ionic-model currents | `sweepCurrents` | full + lightweight; serial | maintained C++ and generated model metadata consumer | A sweep is not equivalent to a spatial solver trajectory | [README](sweepCurrents/README.md) |

Utilities are hand-maintained applications. They may consume generated ionic
Names/equation metadata, but they do not own it. A directory name is not the
build contract: when adding or renaming a utility, keep its `EXE` entry, this
inventory, its local README, and any documented command examples synchronized.
