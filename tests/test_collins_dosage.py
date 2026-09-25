"""Hermetic checks on the committed Collins 2022 dosage test."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "collins_dosage.json")))
P = J["per_mirna"]


def test_gates():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"] and len(P) == 43


def test_counts_and_verdicts():
    for name, a, b in (("D1_phaplo_vs_uniform", "cnn", "uniform_mean"), ("D3_phaplo_vs_utrlen", "cnn", "len_mean"), ("T3_ptriplo_vs_utrlen", "cnn_tri", "len_tri_mean")):
        assert sum(v[a] < v[b] for v in P.values()) == J[name]["wins"]
        assert (J[name]["p_one_sided"] < 0.05) == (J[name]["verdict"] == "CONFIRMED")


def test_score_ranges():
    assert all(0 <= v["cnn"] <= 1 and 0 <= v["cnn_tri"] <= 1 for v in P.values())
