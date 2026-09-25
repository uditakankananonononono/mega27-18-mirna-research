"""Hermetic checks on the committed HPO disease-gene test."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "hpo_disease_genes.json")))
P = J["per_mirna"]


def test_gates():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"] and len(P) == 43


def test_counts_and_verdicts():
    for name, k in (("M1_vs_uniform", "uniform_mean"), ("M2_vs_expr", "expr_mean"), ("M3_vs_utrlen", "len_mean")):
        assert sum(v["cnn"] < v[k] for v in P.values()) == J[name]["wins"]
        assert sum(v["cnn"] for v in P.values()) == J[name]["pooled_cnn"]
        assert (J[name]["p_one_sided"] < 0.05) == (J[name]["verdict"] == "CONFIRMED")
