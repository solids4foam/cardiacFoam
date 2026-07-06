from openfoam_driver.dashboard.case_scan import scan_cases


def test_scan_finds_cases_and_groups_by_family(tutorial_tree):
    cases = {c.case_id: c for c in scan_cases(tutorial_tree)}
    assert "manufacturedSolutions/monodomainPseudoECG" in cases
    assert "PATHOS/RBBB" in cases
    assert cases["PATHOS/RBBB"].family == "PATHOS"


def test_goal_taken_from_readme_body(tutorial_tree):
    cases = {c.case_id: c for c in scan_cases(tutorial_tree)}
    goal = cases["manufacturedSolutions/monodomainPseudoECG"].goal
    assert goal.startswith("Verifies the monodomain pseudo-ECG")


def test_solver_stack_parsed(tutorial_tree):
    cases = {c.case_id: c for c in scan_cases(tutorial_tree)}
    c = cases["manufacturedSolutions/monodomainPseudoECG"]
    assert c.solver == "monodomainSolver"
    assert c.ionic_model == "TNNP"


def test_processor_and_setup_dirs_excluded(tutorial_tree, tmp_path):
    (tutorial_tree / "manufacturedSolutions" / "monodomainPseudoECG" / "processor0"
     / "constant").mkdir(parents=True)
    ids = {c.case_id for c in scan_cases(tutorial_tree)}
    assert not any("processor0" in i for i in ids)


def test_garbled_dict_does_not_raise(tutorial_tree):
    bad = tutorial_tree / "PATHOS" / "RBBB" / "constant" / "electroProperties"
    bad.write_text("\x00\x00 not a dict")
    cases = {c.case_id: c for c in scan_cases(tutorial_tree)}
    assert cases["PATHOS/RBBB"].solver is None
