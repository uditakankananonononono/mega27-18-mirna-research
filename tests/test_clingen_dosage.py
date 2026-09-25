"""Hermetic checks on the committed ClinGen dosage test."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "clingen_dosage.json")))
P = J["per_mirna"]


def test_gates_and_panel():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"] and len(J["panel"]) == 43


def test_verdicts_follow_rule():
    for h in ["C1_cnn_vs_uniform", "C2_cnn_vs_expr", "C3_cnn_vs_utrlen"]:
        r = J[h]
        exp = "INCONCLUSIVE" if r["non_tied"] < 15 else ("CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
        assert r["verdict"] == exp


def test_pooled_counts():
    assert sum(v["cnn"] for v in P.values()) == J["C1_cnn_vs_uniform"]["pooled_cnn"]
    assert all(v["cnn"] <= 200 and v["universe_hi"] >= v["cnn"] for v in P.values())
