"""Hermetic checks on the committed GTEx target-avoidance audit."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "gtex_target_avoidance.json")))


def test_gates_and_design():
    g = J["gates"]
    assert g["G1_ALB_liver"] and g["G2_ACTA1_muscle"] and g["G3_pass"]
    assert len(J["pairs"]) == 8 and J["K"] == 200


def test_cnn_arm_falsified():
    assert J["H1_cnn_depletion"]["verdict"] == "FALSIFIED" and J["H1_cnn_depletion"]["wins"] == 1
    assert J["H3_cnn_tissue_specificity"]["verdict"] == "FALSIFIED"
    P = J["per_mirna"]
    assert sum(v["cnn_matched"] > 0.5 for v in P.values()) == J["H1_cnn_depletion"]["wins"]


def test_preregistered_h2_empty_by_construction():
    h2 = J["H2_targetscan_depletion"]
    assert h2["wins"] == 0 and h2["losses"] == 0 and h2["verdict"] == "FALSIFIED"
    assert all(v["targetscan_n"] == 0 for v in J["per_mirna"].values())


def test_posthoc_reference_recovers_avoidance():
    ph = J["posthoc_targetscan_reference"]
    assert "not pre-registered" in ph["note"]
    assert ph["depletion"]["wins"] == 8 and ph["depletion"]["p_one_sided"] < 0.01
    assert ph["specificity"]["wins"] == 7 and ph["specificity"]["p_one_sided"] < 0.05
    assert all(v["n_set"] == 200 for v in ph["per_mirna"].values())


def test_auc_bounds():
    for v in J["per_mirna"].values():
        for x in v["cnn"].values():
            assert 0.0 <= x["auc_dep"] <= 1.0 and x["n_set"] + x["n_rest"] == v["universe"]
