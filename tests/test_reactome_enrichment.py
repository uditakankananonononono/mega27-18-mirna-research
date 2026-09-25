"""Hermetic checks on the committed Reactome audit (failed gates recorded honestly)."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "reactome_enrichment.json")))


def test_gate_outcomes_recorded():
    g = J["gates"]
    assert g["G1_pass"] is False and g["G1_cellcycle_nsig"] >= 20
    assert g["G2_pass"] is True and g["G3_pass"] is False and g["G3_median_not_found"] >= 0.10


def test_posthoc_review_labelled():
    r = J["posthoc_gate_review"]
    assert "post-hoc" in r["note"] and r["verdict_status"].startswith("INCONCLUSIVE")
    assert any("Cell Cycle" in n for n in r["G1_cellcycle_named_pathways_anywhere"])


def test_hypotheses_not_confirmed():
    for h in ["H1_cnn_vs_uniform", "H2_cnn_vs_sitecount", "H3_cnn_vs_exprmatched"]:
        assert J[h]["verdict"] == "FALSIFIED"
    P = J["per_mirna"]
    assert sum(P[m]["cnn"]["n_sig"] > P[m]["uniform_mean"] for m in P) == J["H1_cnn_vs_uniform"]["wins"]


def test_arm_shapes():
    for v in J["per_mirna"].values():
        assert len(v["matched"]) == 3 and len(v["uniform"]) == 3
        assert 0.0 <= v["cnn"]["not_found_frac"] <= 1.0
