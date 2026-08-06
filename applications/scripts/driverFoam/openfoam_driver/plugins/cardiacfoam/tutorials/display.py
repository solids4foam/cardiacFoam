from openfoam_driver.tutorials_display import TutorialDisplay

TUTORIALS: tuple[TutorialDisplay, ...] = (
    TutorialDisplay(
        id="singleCell",
        title="Single-cell action potential",
        summary=(
            "Run a single-cell sweep over an ionic model to inspect AP "
            "morphology. Useful for drug-effect studies."
        ),
        thumbnail="/tutorials/single-cell.png",
        tags=("single-cell", "ionic-model", "AP"),
        preset={
            "anatomy.mesh": "single-cell",
            "physics.ionic_model": "TenTusscher",
        },
    ),
    TutorialDisplay(
        id="niederer2012",
        title="Niederer 2012 verification benchmark",
        summary=(
            "The Niederer et al. 2012 N-version benchmark for cardiac "
            "tissue electrophysiology. Validates monodomain solvers."
        ),
        thumbnail="/tutorials/niederer-2012.png",
        tags=("benchmark", "verification", "monodomain"),
        preset={
            "anatomy.mesh": "niederer-slab",
            "physics.ionic_model": "TenTusscher",
        },
    ),
    TutorialDisplay(
        id="manufacturedFDA",
        title="Manufactured solution (monodomain)",
        summary=(
            "Method of manufactured solutions on the monodomain "
            "equation. Used to verify spatial/temporal convergence."
        ),
        thumbnail="/tutorials/manufactured-fda.png",
        tags=("manufactured-solution", "verification", "monodomain"),
        preset={
            "anatomy.mesh": "fda-cuboid",
            "physics.ionic_model": "FentonKarma",
        },
    ),
    TutorialDisplay(
        id="manufacturedFDABidomain",
        title="Manufactured solution (bidomain)",
        summary=(
            "Same MMS verification at bidomain resolution. Pairs with "
            "the monodomain variant for cross-formulation comparison."
        ),
        thumbnail="/tutorials/manufactured-fda-bidomain.png",
        tags=("manufactured-solution", "verification", "bidomain"),
        preset={
            "anatomy.mesh": "fda-cuboid",
            "physics.ionic_model": "FentonKarma",
        },
    ),
    TutorialDisplay(
        id="manufacturedFDABathBidomain",
        title="Manufactured solution (bath bidomain)",
        summary=(
            "FDA bidomain-with-bath manufactured solution with a grounded "
            "bath electrode and bath ECG potential verification."
        ),
        thumbnail="/tutorials/manufactured-fda-bath-bidomain.png",
        tags=("manufactured-solution", "verification", "bidomain", "bath-ecg"),
        preset={
            "anatomy.mesh": "fda-bath-cuboid",
            "physics.ionic_model": "FentonKarma",
        },
    ),
    TutorialDisplay(
        id="manufacturedEikonalECG",
        title="Manufactured solution (eikonal ECG)",
        summary=(
            "Manufactured eikonal activation-time verification with template "
            "surrogate ECG and quadrature ECG reference."
        ),
        thumbnail="/tutorials/manufactured-eikonal-ecg.png",
        tags=("manufactured-solution", "verification", "eikonal", "ecg"),
        preset={
            "anatomy.mesh": "unit-domain",
            "physics.ionic_model": "none",
        },
    ),
    TutorialDisplay(
        id="manufacturedMonodomainTotalLagrangianEM",
        title="Manufactured electromechanics (MMS)",
        summary=(
            "Electromechanics verification on a fully coupled manufactured field. "
            "Vm, D, lambda, and Ta are rigorous MMS targets."
        ),
        thumbnail="/tutorials/manufactured-electromechanics-bc.png",
        tags=("manufactured-solution", "verification", "electromechanics"),
        preset={
            "anatomy.mesh": "unit-domain",
            "physics.ionic_model": "monodomainFDAManufactured",
        },
    ),
    TutorialDisplay(
        id="manufacturedMonodomain1D3D",
        title="Manufactured Purkinje-myocardium coupling (MMS)",
        summary=(
            "Coupled 1D Purkinje graph / 3D monodomain manufactured-solution "
            "convergence, across decoupled, unidirectional, and bidirectional "
            "PVJ transfer regimes."
        ),
        thumbnail="/tutorials/manufactured-monodomain-1d3d.png",
        tags=("manufactured-solution", "verification", "purkinje", "coupling"),
        preset={
            "anatomy.mesh": "unit-domain",
            "physics.ionic_model": "monodomainFDAManufactured",
        },
    ),
    TutorialDisplay(
        id="manufacturedPurkinjeGraph",
        title="Manufactured solution (Purkinje graph)",
        summary=(
            "Manufactured monodomain solution on a 1D Purkinje graph coupled "
            "to a 3D domain. Traces mesh refinement convergence on the "
            "Hines-ordered network."
        ),
        thumbnail="/tutorials/manufactured-purkinje-graph.png",
        tags=("manufactured-solution", "verification", "purkinje", "1D-3D"),
        preset={
            "anatomy.mesh": "purkinje-graph",
            "physics.ionic_model": "monodomainFDAManufactured",
        },
    ),
    TutorialDisplay(
        id="heartSolverComparison",
        title="Heart solver comparison (eikonal / monodomain / bidomain)",
        summary=(
            "Compares whole solver stacks -- eikonal, monodomain, bidomain, "
            "and a mixed monodomain-tissue/eikonal-Purkinje variant -- over "
            "one shared real heart anatomy (mesh + Purkinje graph)."
        ),
        thumbnail="/tutorials/heart-solver-comparison.png",
        tags=("real-anatomy", "purkinje", "eikonal", "monodomain", "bidomain", "solver-comparison"),
        preset={
            "anatomy.mesh": "heart-purkinje-graph",
        },
    ),

    TutorialDisplay(
        id="restitutionCurves",
        title="Restitution curves (S1–S2 protocol)",
        summary=(
            "S1–S2 pacing protocol that traces APD restitution. Useful "
            "for arrhythmia substrate studies."
        ),
        thumbnail="/tutorials/restitution-curves.png",
        tags=("single-cell", "S1-S2", "restitution"),
        preset={
            "anatomy.mesh": "single-cell",
            "physics.ionic_model": "TenTusscher",
            "stimulus.protocol": "s1s2",
        },
    ),
    TutorialDisplay(
        id="monodomainAndEikonal1DCableCVConvergence",
        title="1D Cable CV Convergence (Monodomain & Eikonal)",
        summary=(
            "1D cable verification protocol to extract continuous conduction "
            "velocity profiles and perform mesh resolution convergence sweeps."
        ),
        thumbnail="/tutorials/cable-cv-convergence.png",
        tags=("cable", "cv", "convergence", "monodomain", "eikonal"),
        preset={
            "anatomy.mesh": "cable-1d",
            "physics.ionic_model": "BuenoOrovio",
        },
    ),
)
