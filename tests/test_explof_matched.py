"""Hermetic checks on the committed exp_lof-matched control."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "explof_matched.json")))
P = J["per_mirna"]


def test_matching_worked():
    assert 0.95 <= J["matching_check_median_explof_ratio"] <= 1.05


def test_counts_and_verdicts():
    assert sum(v["cnn_loeuf"] > v["m_loeuf"] for v in P.values()) == J["L1_loeuf_vs_explof_matched"]["wins"]
    assert sum(v["cnn_explof"] < v["universe_explof"] for v in P.values()) == J["D1_cnn_shorter_cds"]["wins"]
    for k in ["L1_loeuf_vs_explof_matched", "L2_shet_vs_explof_matched"]:
        assert (J[k]["p_one_sided"] < 0.05) == (J[k]["verdict"] == "CONFIRMED")
