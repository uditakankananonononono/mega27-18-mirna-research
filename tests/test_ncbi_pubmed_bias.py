"""Hermetic checks on the committed NCBI study-bias test and its post-hoc sensitivity."""
import json, os

R = os.path.join(os.path.dirname(__file__), "..", "results")
J = json.load(open(os.path.join(R, "ncbi_pubmed_bias.json")))
Q = json.load(open(os.path.join(R, "ncbi_pubmed_posthoc.json")))


def test_gates_recorded():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"] is False and J["gates"]["G2_map_rate"] < 0.95


def test_verdicts_and_posthoc_label():
    for k in ["P1_less_studied_vs_uniform", "P2_vs_expr", "P3_vs_utrlen", "P4_loeuf_vs_pubmatched"]:
        assert J[k]["verdict"] == "CONFIRMED"
    assert J["P5_hpo_vs_pubmatched"]["verdict"] == "FALSIFIED"
    assert "NOT pre-registered" in Q["note"] and Q["panel"] == J["panel"]


def test_counts_match_rows():
    P = J["per_mirna"]
    assert sum(v["cnn_loeuf"] > v["pub_loeuf_mean"] for v in P.values()) == J["P4_loeuf_vs_pubmatched"]["wins"]
