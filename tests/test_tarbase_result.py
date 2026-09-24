import json
from pathlib import Path
import pytest

P = Path(__file__).resolve().parents[1] / "results" / "tarbase_benchmark.json"


@pytest.mark.skipif(not P.exists(), reason="result not generated")
def test_tarbase_result_consistent():
    j = json.loads(P.read_text())
    assert j["n_rows"] > 0
    for lab in ("any", "low_yield", "clip"):
        r = j["labels"][lab]
        assert r["overlap_with_mirtarbase_pos"] + r["tarbase_pos_among_mirtarbase_neg"] == r["n_pos"]
        for v in r["pooled"].values():
            assert 0.0 <= v <= 1.0
