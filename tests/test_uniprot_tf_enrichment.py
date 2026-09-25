"""Hermetic checks on the committed UniProt TF-enrichment audit."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "uniprot_tf_enrichment.json")))
P = J["per_mirna"]


def test_gates_positive_control():
    g = J["gates"]
    assert g["G1_pass"] and g["G1_targetscan_vs_random"]["wins"] == 13 and g["G2_pass"]
    assert 0.05 <= g["G2_kw0805_prevalence"] <= 0.20


def test_cnn_arm_falsified():
    for h in ["H1_cnn_vs_uniform", "H2_cnn_vs_sitecount", "H3_cnn_vs_exprmatched"]:
        assert J[h]["verdict"] == "FALSIFIED"
    assert sum(P[m]["cnn"] > P[m]["sitecount"] for m in P) == J["H2_cnn_vs_sitecount"]["wins"]


def test_fraction_bounds_and_shapes():
    for v in P.values():
        assert 0 <= v["cnn"] <= 1 and 0 <= v["universe_tf_frac"] <= 1
        assert len(v["matched"]) == 3 and len(v["uniform"]) == 3 and v["ts_n"] == 200
