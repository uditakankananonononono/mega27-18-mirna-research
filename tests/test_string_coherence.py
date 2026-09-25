"""Hermetic checks on the committed STRING network-coherence audit."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "string_coherence.json")))
P = J["per_mirna"]


def test_gates():
    g = J["gates"]
    assert g["G1_pass"] and g["G1_cellcycle_edges"] >= 50
    assert g["G2_pass"] and g["G2_median_random_density"] < 0.01
    assert g["G3_pass"] and g["G3_median_mapped_frac"] >= 0.85


def test_density_definition():
    for v in P.values():
        for x in [v["cnn"], v["clip"], v["expr"]] + v["random"]:
            n = x["n_mapped"]
            assert abs(x["density"] - x["edges"] / (n * (n - 1) / 2)) < 1e-12


def test_preregistered_verdicts():
    assert J["H1_cnn_vs_random"]["verdict"] == "FALSIFIED"
    assert J["H3_cnn_vs_expression"]["verdict"] == "FALSIFIED" and J["H3_cnn_vs_expression"]["losses"] == 13
    assert J["H2_clip_vs_random"]["verdict"] == "CONFIRMED" and J["H2_clip_vs_random"]["wins"] == 13


def test_h2_explained_by_expression():
    ph = J["posthoc_H2_expression_matched_control"]
    assert "not pre-registered" in ph["note"]
    t = ph["clip_vs_matched"]
    assert t["verdict"] == "DOES NOT SURVIVE" and t["p_one_sided"] >= 0.05
    w = sum(v["clip_density"] > v["matched_mean_density"] for v in ph["per_mirna"].values())
    assert w == t["wins"]


def test_matched_sets_sized_like_clip():
    for v in J["posthoc_H2_expression_matched_control"]["per_mirna"].values():
        assert all(abs(s - v["clip_size"]) <= 2 for s in v["set_sizes"])
