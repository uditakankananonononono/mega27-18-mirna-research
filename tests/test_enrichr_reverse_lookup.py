"""Hermetic checks on the committed Enrichr reverse-lookup audit."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "enrichr_reverse_lookup.json")))
P = J["per_mirna"]


def test_gates():
    g = J["gates"]
    assert g["G1_pass"] and g["G1_miR21_rank"] == 1 and g["G2_pass"] and g["G3_pass"]
    assert J["library"] == "miRTarBase_2017" and J["library_terms"] == 3240


def test_preregistered_verdicts():
    assert J["H1_cnn_vs_random"]["verdict"] == "CONFIRMED" and J["H1_cnn_vs_random"]["wins"] == 13
    assert J["H2_clip_vs_random"]["verdict"] == "CONFIRMED"
    assert J["H3_cnn_vs_expression"]["verdict"] == "FALSIFIED"
    assert sum(P[m]["cnn"]["rank"] < P[m]["rand_mean_rank"] for m in P) == J["H1_cnn_vs_random"]["wins"]


def test_posthoc_controls_labelled_and_survive():
    sc, em = J["posthoc_sitecount_control"], J["posthoc_expression_matched_control"]
    assert "not pre-registered" in sc["note"] and "not pre-registered" in em["note"]
    assert sc["cnn_vs_sitecount"]["verdict"] == "SURVIVES" and sc["cnn_vs_sitecount"]["losses"] == 0
    assert em["cnn_vs_matched"]["verdict"] == "SURVIVES" and em["cnn_vs_matched"]["wins"] == 13


def test_rank_bounds():
    for v in P.values():
        for x in [v["cnn"], v["clip"], v["expr"]] + v["random"]:
            assert 1 <= x["rank"] <= x["n_terms"] + 1


def test_own_rank_helper():
    import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    import pytest; pytest.importorskip("requests")
    from enrichr_reverse_lookup import own_rank
    res = [[1, "a", 0.5, 0, 0, [], 0.9], [2, "b", 0.01, 0, 0, [], 0.1]]
    assert own_rank(res, "b")["rank"] == 1 and own_rank(res, "a")["rank"] == 2 and own_rank(res, "c")["rank"] == 3


def test_heldout_replication_recorded():
    R = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "enrichr_replication.json")))
    assert len(R["panel"]) == 30 and not set(R["panel"]) & set(J["panel"])
    assert R["R1_cnn_vs_uniform"]["verdict"] == "CONFIRMED" and R["R3_cnn_vs_exprmatched"]["verdict"] == "CONFIRMED"
    assert R["R2_cnn_vs_sitecount"]["verdict"] == "FALSIFIED" and R["discovery_replicates"] is False
