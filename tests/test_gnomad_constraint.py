"""Hermetic checks on the committed gnomAD constraint audit."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "gnomad_constraint.json")))
P = J["per_mirna"]


def test_gates():
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"] and J["gates"]["G2_coverage"] >= 0.85


def test_preregistered_falsified_and_reversal_consistent():
    for h in ["H1_cnn_vs_uniform", "H2_cnn_vs_sitecount", "H3_cnn_vs_exprmatched"]:
        assert J[h]["verdict"] == "FALSIFIED"
    assert sum(P[m]["cnn"] > P[m]["uniform_mean"] for m in P) == J["post_hoc_not_preregistered"]["cnn_higher_loeuf_than_uniform"]
    assert J["post_hoc_not_preregistered"]["p_two_sided_uniform"] <= 1


def test_positive_control_direction():
    assert sum(P[m]["ts_median"] < sum(P[m]["ts_rand"]) / 3 for m in P) == J["gates"]["G1_targetscan_vs_random"]["wins"]
