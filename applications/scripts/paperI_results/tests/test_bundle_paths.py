from pathlib import Path

SH = Path(__file__).resolve().parents[4] / "tutorials/manufacturedSolutions/bundle_paperI_data.sh"

def test_bundle_pulls_from_reference_not_results():
    text = SH.read_text()
    assert "/reference/" in text
    assert "*/setup/results/*_convergence.csv" not in text
