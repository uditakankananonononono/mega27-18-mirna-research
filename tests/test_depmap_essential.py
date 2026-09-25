"""Hermetic checks on the committed DepMap essential-gene test."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "depmap_essential.json")))
P = J["per_mirna"]


def test_gates_recorded():
    assert J["gates"]["G1_pass"] is False and J["gates"]["G2_pass"] and len(P) == 43


def test_hypotheses_falsified_and_counts():
    for name, k in (("E1_vs_uniform", "uniform_mean"), ("E2_vs_expr", "expr_mean"), ("E3_vs_utrlen", "len_mean")):
        assert J[name]["verdict"] == "FALSIFIED"
        assert sum(v["cnn"] < v[k] for v in P.values()) == J[name]["wins"]
