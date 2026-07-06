from openfoam_driver.dashboard.models import CaseCard, CaseView


def test_casecard_defaults_and_serialization():
    card = CaseCard(case_id="PATHOS/RBBB", family="PATHOS", name="RBBB")
    assert card.status == "not_run"
    assert card.regression == "none"
    assert card.metrics == {}
    assert card.plots == []
    assert card.is_3d is False
    d = card.to_dict()
    assert d["case_id"] == "PATHOS/RBBB"
    assert d["metrics"] == {}


def test_caseview_wraps_card():
    card = CaseCard(case_id="x/y", family="x", name="y")
    view = CaseView(card=card, user_notes="hi", tags=["final"], captions={"p.png": "cap"})
    assert view.tags == ["final"]
    assert view.captions["p.png"] == "cap"
