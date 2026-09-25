"""Hermetic checks on the committed g:Profiler functional-coherence audit."""
import json, os
import pytest

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "gprofiler_enrichment.json")))


def test_gates_pass():
    g = J["gates"]
    assert g["G1_cellcycle_GO0007049"] and g["G2_pass"] and g["G3_pass"]
    assert g["G2_median_random_nsig"] <= 2 and g["G3_median_unmapped"] < 0.10


def test_panel_shape():
    assert J["K"] == 200 and J["R"] == 3 and len(J["panel"]) == 13
    for m in J["panel"]:
        v = J["per_mirna"][m]
        assert len(v["random"]) == 3 and v["universe"] > J["K"]


def test_hypotheses_falsified_and_consistent():
    for h in ["H1_cnn_vs_random", "H2_termJaccard_cnn_vs_random", "H3_cnn_vs_expression"]:
        assert J[h]["verdict"] == "FALSIFIED" and J[h]["p_one_sided"] >= 0.05
    P = J["per_mirna"]
    w = sum(P[m]["cnn"]["n_sig"] > P[m]["rand_mean_nsig"] for m in P)
    assert w == J["H1_cnn_vs_random"]["wins"]
    assert J["H3_cnn_vs_expression"]["losses"] == sum(P[m]["expr"]["n_sig"] > P[m]["cnn"]["n_sig"] for m in P)


def test_expression_control_dominates():
    P = J["per_mirna"]
    assert min(P[m]["expr"]["n_sig"] for m in P) > max(P[m]["cnn"]["n_sig"] for m in P)


def test_jaccard_bounds_and_posthoc_labelled():
    for v in J["per_mirna"].values():
        assert 0.0 <= v["J_cnn_clip"] <= 1.0 and 0.0 <= v["J_rand_clip"] <= 1.0
    assert "not pre-registered" in J["posthoc_genome_background"]["note"]


def test_jaccard_helper():
    import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    pytest.importorskip("requests")
    from gprofiler_enrichment import jaccard, signtest
    assert jaccard([], []) == 0.0 and jaccard("ab", "bc") == pytest.approx(1 / 3)
    assert signtest(0, 0)["p_one_sided"] == 1.0 and signtest(10, 0)["p_one_sided"] < 0.01
