"""Hermetic checks on the committed ClinVar PLP-gene test."""
import json, os

J = json.load(open(os.path.join(os.path.dirname(__file__), "..", "results", "clinvar_pathogenic.json")))
P = J["per_mirna"]


def test_positive_control_gate_failed_and_recorded():
    g = J["gates"]
    assert g["G1_pass"] is False and g["G1_ts_vs_random"]["p_one_sided"] >= 0.05 and g["G2_pass"]


def test_counts_match_rows():
    for name, a, k in (("V1_vs_uniform", "cnn", "uniform_mean"), ("V3_vs_utrlen", "cnn", "len_mean"), ("S3_any_vs_utrlen", "cnn_any", "len_any_mean")):
        assert sum(v[a] < v[k] for v in P.values()) == J[name]["wins"]
    assert all(v["cnn"] <= v["cnn_any"] <= 200 for v in P.values())
