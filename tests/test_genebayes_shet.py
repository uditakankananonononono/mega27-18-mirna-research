"""Hermetic checks on the committed GeneBayes s_het test."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "genebayes_shet.json")))


def test_gates_and_verdicts():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"]
    for k in ["H1_vs_uniform", "H2_vs_expr", "H3_vs_utrlen", "H4_vs_pubmatched"]:
        assert (J[k]["p_one_sided"] < 0.05) == (J[k]["verdict"] == "CONFIRMED")


def test_counts_and_log_scale():
    P = J["per_mirna"]
    assert sum(v["cnn"] < v["len_mean"] for v in P.values()) == J["H3_vs_utrlen"]["wins"]
    assert all(v["cnn"] < 0 for v in P.values())  # log10 of a selection coefficient in (0, 1)
