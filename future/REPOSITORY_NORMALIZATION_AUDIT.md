# Repository Normalization Characteristics and Audit

This document records the implementation characteristics that should be
preserved when the cardiacFoam structure is transposed to another solver or
repository. It is a normalization contract, not a second tutorial manual.
Detailed mathematical justification, experiment-specific decisions, and
failure analysis belong in the README.md of the relevant tutorial or study.

## 1. Main C++ objectives

The C++ implementation has one primary architectural objective: keep the
solver libraries generic and reusable while making tutorial behavior explicit
through OpenFOAM dictionaries.

The corresponding rules are:

- Runtime selection is the public extension boundary. Models, solvers,
  couplers, ECG solvers, and verification models are selected by dictionary
  type/selector entries and registered in the appropriate runtime table.
- C++ libraries own algorithms, interfaces, field ownership, and dictionary
  parsing. A tutorial owns its manufactured formula, mesh-specific setup,
  expected values, and experiment orchestration.
- Tutorial-specific C++ is permitted only when it is genuinely part of the
  tutorial experiment. The electromechanics MMS case is the precedent: its
  local boundary condition, body-force option, symbolic expression, and
  verification harness live under that case's src/ and verification/.
- Dictionary reads must correspond to real, documented input leaves. Avoid
  hidden tutorial names, duplicated constants, and solver branches that can
  only be understood by reading one particular case.
- Shared verification models should expose a stable interface and be selected
  from the case dictionary. The case README explains what is verified and why;
  the verifier implements the measurement and output contract.
- Generated ionic equations and other generated source must be changed through
  their generator/template pipeline, not normalized by hand in generated
  headers.

The runtime flow to preserve is:

~~~
top-level physics selector
  -> electro/physics model selector
  -> solver or domain selector
  -> optional coupler / ECG / verification model selector
  -> numerical kernel
~~~

The main public surfaces are documented in
[src/ARCHITECTURE.md](../src/ARCHITECTURE.md) and
[src/electroModels/ARCHITECTURE.md](../src/electroModels/ARCHITECTURE.md).

## 2. Dictionary normalization

### 2.1 Comment policy

OpenFOAM dictionaries are executable input, not design documents.

- Keep the standard FoamFile header when required by the repository's
  OpenFOAM convention.
- Do not place multi-line explanations, derivations, historical context, or
  troubleshooting notes in a dictionary.
- If a semantic comment is necessary, keep it to one short line adjacent to
  the affected entry.
- Move all justification to the case or study README.md.
- Prefer self-describing key names, units, and stable block structure over
  explanatory comments.

Good:

~~~
dimension "3D";  // Selects the 3-D manufactured solution.
~~~

Not canonical:

~~~
/* Several paragraphs explaining the history of this parameter,
   alternative formulations, and why this value was chosen. */
dimension "3D";
~~~

The standard OpenFOAM banner is an allowed format header; it is not a
substitute for case documentation.

### 2.2 Stable dictionary surfaces

For an electro-only tutorial, keep the input surface predictable:

~~~
constant/
├── physicsProperties       # top-level physics family
└── electroProperties       # electro selector and solver coefficients
~~~

The normal selection pattern is:

~~~
// constant/physicsProperties
type electroModel;

// constant/electroProperties
myocardiumSolver monodomainSolver;
monodomainSolverCoeffs
{
    ...
}
~~~

The selected solver's settings belong in the matching <selector>Coeffs block.
Nested blocks represent configuration groups; they are not separate top-level
configuration files unless the runtime requires a separate region or OpenFOAM
object.

For a coupled electromechanics tutorial, preserve the explicit region split:

~~~
constant/
├── physicsProperties
├── electroMechanicalProperties
├── electro/
│   └── electroProperties
└── solid/
    ├── mechanicalProperties
    └── solidProperties

system/
├── controlDict
├── electro/
│   ├── fvSchemes
│   └── fvSolution
└── solid/
    ├── fvSchemes
    └── fvSolution
~~~

This region split is a deliberate variant, not a reason to flatten the
electro-only cases.

### 2.3 Generalization rules for inputs

- Use the same key names and block nesting for the same runtime concept in
  every tutorial.
- Put case values in dictionaries; put the reason for those values in the
  README.
- Keep verification selectors, manufactured-solution parameters, and output
  selections inside the solver coefficient block that consumes them.
- Keep mesh variants in system/ with explicit names such as
  blockMeshDict.1D, blockMeshDict.2D, and blockMeshDict.3D.
- Treat generated defaults, meshes, processor directories, logs, and
  post-processing output as derived artifacts. They must be ignored or
  removed by Allclean; the tracked dictionary is the source of truth.
- When a parameter is shared by C++ components, read it from one dictionary
  location. Do not repeat a physics value in a verifier, a script, and a
  tutorial dictionary.

## 3. Canonical tutorial architecture

The maintained tutorial index is
[tutorials/README.md](../tutorials/README.md). A normal runnable case follows
this shape:

~~~
<tutorial>/
├── README.md                 # purpose, inputs, rationale, outputs, usage
├── Allrun                    # canonical serial/parallel entry point
├── Allclean                  # removes derived case output
├── 0/                        # initial fields, when required
├── constant/                 # runtime dictionaries and tracked static data
├── system/                   # control, discretization, mesh, sampling
├── regression/               # case-level reference and regression checker
├── reference/                # frozen convergence/reference data
├── setup/                    # automation, post-processing, studies
└── postProcessing/           # generated output; never the source of truth
~~~

Not every case needs every directory. The invariant is that the case root is
self-contained and its README is the entry point.

### Manufactured-solution cases

The manufactured-solution family is grouped under
tutorials/manufacturedSolutions/ by the verified physical scope:

| Case | Main purpose | Structural characteristic |
| --- | --- | --- |
| monodomainPseudoECG | Monodomain field and pseudo-ECG verification | Monodomain coefficients plus ECG domain |
| bidomain | Bidomain field convergence | Intracellular/extracellular coefficients |
| bathBidomain | Bidomain with bath and ECG ownership | Bath domain, interface, and conductivity inputs |
| eikonalECG | Activation-time and ECG verification | Eikonal coefficients plus ECG domain |
| monodomain1D3D | 1-D Purkinje / 3-D monodomain coupling | Graph input and domainCouplings block |
| monodomainTotalLagrangianEM | Coupled electromechanics MMS | Electro and solid regions plus case-local C++ |

The first four are the common electro-only pattern. The last two are
specializations that retain the same case-level documentation and automation
contract while adding domain-specific inputs.

### Study directories

Each convergence or sensitivity study belongs below the owning case:

~~~
setup/studies/<studyName>/
├── README.md                 # study purpose and interpretation
├── sweep_*.json              # machine-readable sweep definition
├── run/ or scripts           # optional orchestration helpers
├── aggregate*.py             # optional result aggregation
├── results/                  # generated and ignored
└── .gitignore
~~~

The study README states the question being tested, command to run it, output
location, and whether the result is canonical or exploratory. It does not move
long explanations into JSON or OpenFOAM dictionaries.

## 4. Documentation placement

Use the nearest README that owns the explanation:

| Information | Location |
| --- | --- |
| Repository-wide architecture | root README.md |
| Library/runtime ownership | src/*/README.md or ARCHITECTURE.md |
| Case purpose and scientific rationale | tutorial README.md |
| Sweep definition and interpretation | setup/studies/<name>/README.md |
| Script mechanics | script comments/docstrings, kept concise |
| Input values and runtime selection | OpenFOAM dictionaries |
| Frozen numerical expectations | reference/ and regression files |

The rule is simple: dictionaries say what to run; READMEs say why it is
organized that way.

## 5. Audit of the current repository

This audit was performed against the current working tree and is intended to
be kept with the normalization contract.

| Check | Result | Evidence / action |
| --- | --- | --- |
| Common electro-only dictionary surface | Pass for maintained manufactured cases | physicsProperties selects electroModel; electroProperties selects the solver and owns the matching Coeffs block. |
| Regioned electromechanics input surface | Pass | monodomainTotalLagrangianEM separates constant/electro and constant/solid, with matching system regions. |
| Case-local mathematical justification | Pass in representative manufactured cases | Rationale is in case and study READMEs, not in the maintained manufactured dictionaries. |
| Semantic dictionary comments are one-line | Mostly pass, with a known exception | Manufactured case dictionaries contain only the OpenFOAM banner/separators; tutorials/template/constant/electroProperties contains extensive explanatory comments and is not canonical under this rule. |
| Reusable C++ versus tutorial-local code | Pass with deliberate MMS extension | Shared verifiers live under src/verificationModels; electromechanics-only MMS support lives under the tutorial's local src/. |
| Generated artifacts separated from inputs | Pass by convention | Allclean, .gitignore, reference/, and setup/studies/*/results/ establish the source/derived boundary. |
| Manufactured-solution index completeness | Needs update | tutorials/manufacturedSolutions/README.md lists four cases but the tree also contains monodomain1D3D and monodomainTotalLagrangianEM. |
| Top-level tutorial index completeness | Needs review | tutorials/README.md documents the canonical suite but should be checked whenever a manufactured case is added or renamed. |

The two Needs update entries are documentation drift, not a reason to alter
the input architecture. The template comment volume is normalization debt;
its explanatory content should move to a template README or a separate input
reference when the template is next revised.

## 6. Transposition checklist

When applying this contract to another repository, verify the following in
order:

1. Identify the C++ runtime-selection layers and write down the dictionary
   entry that selects each layer.
2. Define one stable input surface for the common case and explicit region
   variants only where the solver requires them.
3. Create one self-contained tutorial directory with README, run/clean entry
   points, inputs, system files, and references before adding study variants.
4. Keep dictionary comments to one line and move all rationale to Markdown.
5. Put every sweep under its owning tutorial, with a README and ignored result
   directory.
6. Separate reusable C++ algorithms from tutorial-specific manufactured
   formulas, boundary conditions, and expected-value logic.
7. Audit the tutorial index and README index whenever a case is added.
8. Run the case's regression check and confirm that cleaning the case removes
   only generated data.

The resulting repository should be understandable from its READMEs and runtime
selectors without requiring historical context from dictionary comments.
