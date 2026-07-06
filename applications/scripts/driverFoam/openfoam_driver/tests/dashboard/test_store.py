from openfoam_driver.dashboard.store import AnnotationStore


def test_annotation_roundtrip(store_db):
    s = AnnotationStore(store_db)
    s.upsert_annotation("PATHOS/RBBB", notes="looks good", tags=["Final", "ischemia paper"])
    a = s.get_annotation("PATHOS/RBBB")
    assert a.notes == "looks good"
    assert a.tags == ["final", "ischemia-paper"]  # normalized


def test_empty_annotation_default(store_db):
    s = AnnotationStore(store_db)
    a = s.get_annotation("nope/x")
    assert a.notes == ""
    assert a.tags == []


def test_caption_roundtrip(store_db):
    s = AnnotationStore(store_db)
    s.upsert_caption("PATHOS/RBBB", "ecg_plots.png", "Simulated V1-V6")
    caps = s.get_captions("PATHOS/RBBB")
    assert caps["ecg_plots.png"] == "Simulated V1-V6"


def test_all_tags_union(store_db):
    s = AnnotationStore(store_db)
    s.upsert_annotation("a/b", notes="", tags=["final"])
    s.upsert_annotation("c/d", notes="", tags=["draft", "final"])
    assert s.all_tags() == ["draft", "final"]
