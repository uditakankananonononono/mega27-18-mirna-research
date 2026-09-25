"""Hermetic checks on the committed pre-registered gnomAD replication."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "gnomad_replication.json")))
P = J["per_mirna"]
DISCOVERY = {"hsa-let-7a-5p", "hsa-miR-1-3p", "hsa-miR-122-5p", "hsa-miR-125b-5p", "hsa-miR-145-5p", "hsa-miR-155-5p", "hsa-miR-16-5p",
             "hsa-miR-17-5p", "hsa-miR-19a-3p", "hsa-miR-200c-3p", "hsa-miR-21-5p", "hsa-miR-34a-5p", "hsa-miR-92a-3p"}


def test_heldout_disjoint_and_gates():
    assert not DISCOVERY & set(J["panel"]) and len(J["panel"]) == 30
    assert J["gates"]["G1_pass"] and J["gates"]["G2_pass"]


def test_sign_counts_match_rows():
    for k, h in [("uniform_mean", "R1_cnn_vs_uniform"), ("expr_mean", "R2_cnn_vs_expr"), ("len_mean", "R3_cnn_vs_utrlen")]:
        assert sum(P[m]["cnn"] > P[m][k] for m in P) == J[h]["wins"]
        assert (J[h]["p_one_sided"] < 0.05) == (J[h]["verdict"] == "CONFIRMED")
