"""Hermetic checks on the committed OmniPath independent-label audit."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "omnipath_independent.json")))
P = J["per_mirna"]


def test_gates():
    g = J["gates"]
    assert g["G1_pass"] and g["G1_n_eligible"] == len(P) >= 10
    assert g["G2_pass"] and 0.45 <= g["G2_median_perm_mean"] <= 0.55


def test_labels_exclude_mirtarbase():
    assert J["n_rows_non_mirtarbase"] < J["n_rows_all"]
    assert all(v["n_pos"] >= 5 for v in P.values())


def test_verdicts():
    assert J["H1_cnn_above_chance"]["verdict"] == "CONFIRMED"
    assert J["H2_cnn_vs_sitecount"]["verdict"] == "FALSIFIED" and J["H3_cnn_vs_expr"]["verdict"] == "FALSIFIED"
    assert sum(v["auroc_cnn"] > v["auroc_sitecount"] for v in P.values()) == J["H2_cnn_vs_sitecount"]["wins"]


def test_auroc_bounds():
    for v in P.values():
        for k in ["auroc_cnn", "auroc_sitecount", "auroc_expr", "perm_mean_cnn"]:
            assert 0.0 <= v[k] <= 1.0
