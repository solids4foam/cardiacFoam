import pytest

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from openfoam_driver.dashboard.app import create_app  # noqa: E402


@pytest.fixture
def client(tutorial_tree, store_db):
    app = create_app(root=tutorial_tree, store_path=store_db, enable_3d=False)
    return TestClient(app)


def test_inventory_lists_cases(client):
    r = client.get("/")
    assert r.status_code == 200
    assert "monodomainPseudoECG" in r.text
    assert "RBBB" in r.text


def test_case_detail_shows_goal_and_metrics(client):
    r = client.get("/case/manufacturedSolutions/monodomainPseudoECG")
    assert r.status_code == 200
    assert "pseudo-ECG" in r.text
    assert "ecg_peak_to_peak" in r.text


def test_3d_case_has_model_viewer(tutorial_tree, store_db):
    case = tutorial_tree / "PATHOS" / "RBBB"
    (case / "heart.vtu").write_text("x")
    app = create_app(root=tutorial_tree, store_path=store_db, enable_3d=False)
    client = TestClient(app)
    r = client.get("/case/PATHOS/RBBB")
    assert "<model-viewer" in r.text


def test_annotation_persists(client):
    r = client.post(
        "/case/PATHOS/RBBB/annotation",
        data={"notes": "final run", "tags": "ischemia paper, final"},
        follow_redirects=False,
    )
    assert r.status_code in (302, 303)
    detail = client.get("/case/PATHOS/RBBB")
    assert "final run" in detail.text
    assert "ischemia-paper" in detail.text


def test_caption_persists(client):
    client.post(
        "/case/manufacturedSolutions/monodomainPseudoECG/caption",
        data={"figure_id": "ecg.png", "caption": "Lead trace"},
        follow_redirects=False,
    )
    caps = client.get("/api/catalog.json").json()
    assert caps["cases"]  # smoke: catalog serializes


def test_rescan_picks_up_new_case(client, tutorial_tree):
    new = tutorial_tree / "electrophysiologyProtocols" / "singleCell"
    (new / "constant").mkdir(parents=True)
    (new / "system").mkdir(parents=True)
    (new / "system" / "controlDict").write_text("application singleCellSolver;\n")
    (new / "README.md").write_text("# singleCell\n\nSingle cell ODE.\n")
    client.post("/rescan", follow_redirects=False)
    r = client.get("/")
    assert "singleCell" in r.text
